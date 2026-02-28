#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

// ===========================================================================
// Paetsch et al. (2012) Constitutive Model
// Replicating Figure 5: Passive to Active Transition
// ===========================================================================

class PaetschModel {
private:
    // --- Material Parameters (Table 1) ---
    // 単位: Pa (論文中のkPaは *1000 して統一)
    double mu;          // 755.5 Pa
    double mp;          // 86.13 (無次元)
    double l0;          // 0.95 (Active reference stretch)
    double c1;          // 470.0 (Active stiffness magnitude)
    double c2;          // 0.95 (Active overlap width)
    
    // Hysteresis Parameters (Table 1)
    // units for m are kPa -> converted to Pa
    double r_pass, r_act;
    double m_pass, m_act; 

public:
    PaetschModel() {
        mu = 755.5;           
        mp = 86.13;           
        l0 = 0.95;            
        c1 = 470.0;           
        c2 = 0.95;            

        // Unloading parameters (kPa -> Pa)
        r_pass = 2.5;   r_act = 30.0;
        m_pass = 2500.0; m_act = 30000.0; 
    }

    // --- Active Stiffness ma (Eq. 73) ---
    // 伸展 lambda に依存するActive剛性
    double calc_ma(double lambda) {
        // term = (lambda^2 / l0^2) - l0^2
        double term = (std::pow(lambda, 2.0) * std::pow(l0, -2.0)) - std::pow(l0, 2.0);
        return c1 * std::exp( - term / c2 );
    }

    // --- Strain Energy Density W0 (Eq. 72) ---
    // 仮想的な弾性エネルギー（ヒステリシス適用前）
    double calc_W0(double lambda, double alpha) {
        double I1 = std::pow(lambda, 2.0) + 2.0/lambda;
        double I4 = std::pow(lambda, 2.0);
        double hat_I4 = I4 * std::pow(l0, -2.0); // I4 / I40
        double I40 = std::pow(l0, 2.0);

        // 1. Matrix (Neo-Hookean)
        double W_iso = (mu / 2.0) * (I1 - 3.0);

        // 2. Passive Fiber
        double W_pass = (mu * mp / 2.0) * (1.0 - alpha) * std::pow(I4 - 1.0, 2.0);

        // 3. Active Fiber
        // Eq. 72: (mu * ma / 2) * alpha * (hat_I4 - I40)^2
        // hat_I4 - I40 = (lambda^2/l0^2) - l0^2
        double ma = calc_ma(lambda);
        double term_act = hat_I4 - I40;
        double W_act = (mu * ma / 2.0) * alpha * std::pow(term_act, 2.0);

        return W_iso + W_pass + W_act;
    }

    // --- Elastic Cauchy Stress sigma0 (Eq. 75) ---
    // ヒステリシス適用前の基本応力
    double calc_sigma0(double lambda, double alpha) {
        double ma = calc_ma(lambda);
        double I4 = std::pow(lambda, 2.0);
        double I40 = std::pow(l0, 2.0);
        double hat_I4 = I4 / I40;

        // Term 1: Matrix
        double s_iso = mu * (I4 - 1.0/lambda); 

        // Term 2: Passive Fiber
        // 2 * mp * mu * lambda^2 * (1-alpha) * (lambda^2 - 1)
        double s_pass = 2.0 * mp * mu * I4 * (1.0 - alpha) * (I4 - 1.0);

        // Term 3: Active Fiber
        // 微分計算が複雑なので Eq.75 をそのまま実装
        // term_bracket = (l^2 * l0^-2 - l0^2)
        double term_bracket = hat_I4 - I40;
        
        // Eq. 75 active part:
        // 2 * mu * ma * alpha * l0^-2 * lambda^2 * (term_bracket) * [1 - (term_bracket)/c2]
        // ※論文Eq.75の最後の項 [1 - .../c2] は指数の微分から来る項
        double s_act = 2.0 * mu * ma * alpha * std::pow(l0, -2.0) * I4 * term_bracket * (1.0 - term_bracket / c2);

        return s_iso + s_pass + s_act;
    }

    // --- Damage Parameter Z (Eq. 76, 82) ---
    // Unloading時のヒステリシス計算
    double calc_Z(double lambda, double alpha, double W0_max, double W0_current) {
        // Active/Passive混合パラメータ (Eq. 82)
        double r = (1.0 - alpha) * r_pass + alpha * r_act;
        double m = (1.0 - alpha) * m_pass + alpha * m_act;

        double energy_diff = W0_max - W0_current;
        if (energy_diff < 0) energy_diff = 0; // 数値誤差対策

        double arg = energy_diff / m;
        return 1.0 - (1.0/r) * std::tanh(arg);
    }
};

int main() {
    PaetschModel model;
    std::ofstream file("passive_active_transition.csv");
    file << "Phase,Stretch,Alpha,Stress,Energy" << std::endl;

    // --- Simulation Settings ---
    double lambda_start = 1.0;
    double lambda_max = 1.15;
    int steps = 100;
    
    double W0_max_at_turnaround = 0.0; // 折り返し地点でのエネルギー

    // ==========================================
    // Phase 1: Passive Loading (alpha = 0)
    // lambda: 1.0 -> 1.15
    // ==========================================
    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        double lambda = lambda_start + (lambda_max - lambda_start) * t;
        double alpha = 0.0;

        // Loading中は Z = 1.0 なので sigma = sigma0
        double stress = model.calc_sigma0(lambda, alpha);
        double energy = model.calc_W0(lambda, alpha);

        // 折り返し地点（最大伸展）でのエネルギーを保存
        if (i == steps) {
            W0_max_at_turnaround = energy;
        }

        file << "Loading," << lambda << "," << alpha << "," << stress << "," << energy << std::endl;
    }

    // ==========================================
    // Phase 2: Active Unloading (Transition)
    // lambda: 1.15 -> 1.0
    // alpha:  0.0  -> 1.0 (徐々にActive化)
    // ==========================================
    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        double lambda = lambda_max - (lambda_max - lambda_start) * t;
        
        // Transition: 単純な線形ランプ関数で遷移させる
        // Figure 5の挙動を再現するため、Unloadingの進行に合わせてalphaを上げる
        double alpha = t; 

        // 1. 現在の状態での仮想的な弾性エネルギー W0
        double W0_curr = model.calc_W0(lambda, alpha);

        // 2. 仮想的な弾性応力 sigma0
        double sigma0 = model.calc_sigma0(lambda, alpha);

        // 3. ダメージパラメータ Z (Eq. 76)
        // W0_max_at_turnaround は Phase 1 のピーク値を使う
        // （論文Eq.76の W0(l_max) は loading path 上の最大値を指すため）
        double Z = model.calc_Z(lambda, alpha, W0_max_at_turnaround, W0_curr);

        // 4. 実際の応力 (Eq. 77)
        double stress = Z * sigma0;

        file << "Unloading," << lambda << "," << alpha << "," << stress << "," << W0_curr << std::endl;
    }

    file.close();
    std::cout << "Data saved to passive_active_transition.csv" << std::endl;
    return 0;
}