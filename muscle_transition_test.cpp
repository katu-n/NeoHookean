#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>

// 論文のTable 1に基づくパラメータ 
struct MaterialParams {
    double mu;      // Matrix shear modulus (Pa)
    double mu_p;    // Passive anisotropy parameter
    double lambda_0;// Reference stretch for active fibers
    double c1;      // Active stiffness parameter
    double c2;      // Active stiffness parameter
    double r;       // Damage parameter (Hysteresis)
    double m;       // Hysteresis slope parameter (Pa)
    double a;       // Reloading parameter (Pa)
    double b;       // Secondary unloading parameter (Pa)
};

class MuscleModel {
private:
    MaterialParams p_pass; 
    MaterialParams p_act;  
    MaterialParams p_curr; 

    double alpha;      // Activation level (0.0 to 1.0)
    double lambda_max; // 負荷履歴の最大伸展比
    double W0_max;     // 負荷履歴の最大エネルギー

public:
    MuscleModel() {
        // Table 1 values
        p_pass = {755.5, 86.13, 0.95, 470.0, 0.95, 1.05, 2500.0, 2500.0, 2500.0}; 
        p_act  = {755.5, 86.13, 0.95, 470.0, 0.95, 1.55, 30000.0, 30000.0, 15000.0};
        
        alpha = 0.0;
        update_params();
        lambda_max = 1.0;
        W0_max = 0.0; 
    }

    // Alphaに基づいてパラメータを線形補間する Eq. 82
    void set_activation(double a_val) {
        alpha = std::max(0.0, std::min(1.0, a_val));
        update_params();
    }

    void update_params() {
        p_curr.mu       = p_pass.mu; 
        p_curr.mu_p     = p_pass.mu_p;
        p_curr.lambda_0 = p_pass.lambda_0;
        p_curr.c1       = p_pass.c1;
        p_curr.c2       = p_pass.c2;

        // ヒステリシス関連パラメータの補間
        p_curr.r = (1.0 - alpha) * p_pass.r + alpha * p_act.r;
        p_curr.m = (1.0 - alpha) * p_pass.m + alpha * p_act.m;
        p_curr.a = (1.0 - alpha) * p_pass.a + alpha * p_act.a;
        p_curr.b = (1.0 - alpha) * p_pass.b + alpha * p_act.b;
    }

    // Active stiffness function Eq. 43
    double get_mu_a(double lambda) const {
        double I4_hat = std::pow(lambda / p_curr.lambda_0, 2.0);
        double I4_0   = std::pow(p_curr.lambda_0, 2.0);
        double exponent = -(I4_hat - I4_0) / p_curr.c2;
        return p_curr.c1 * std::exp(exponent);
    }

    // Total Strain Energy W0 Eq. 44
    double get_W0(double lambda) const {
        // 1. Isotropic Matrix (Neo-Hookean base)
        double I1 = std::pow(lambda, 2.0) + 2.0 / lambda;
        double W_iso = (p_curr.mu / 2.0) * (I1 - 3.0);

        // 2. Passive Fiber
        double W_pass = 0.0;
        if (lambda > 1.0) {
            double I4 = std::pow(lambda, 2.0);
            W_pass = (p_curr.mu * p_curr.mu_p / 2.0) * (1.0 - alpha) * std::pow(I4 - 1.0, 2.0);
        }

        // 3. Active Fiber
        double W_act = 0.0;
        double I4_hat = std::pow(lambda / p_curr.lambda_0, 2.0);
        double I4_0   = std::pow(p_curr.lambda_0, 2.0);
        double mu_a_val = get_mu_a(lambda);
        W_act = (p_curr.mu * mu_a_val / 2.0) * alpha * std::pow(I4_hat - I4_0, 2.0);

        return W_iso + W_pass + W_act;
    }

    // Stress on the primary loading path (sigma_0) Eq. 75
    double get_sigma0(double lambda) const {
        // Isotropic part
        double sigma_iso = p_curr.mu * (std::pow(lambda, 2.0) - 1.0 / lambda);

        // Passive part
        double sigma_pass = 0.0;
        if (lambda > 1.0) {
            sigma_pass = 2.0 * p_curr.mu * p_curr.mu_p * (1.0 - alpha) * std::pow(lambda, 2.0) * (std::pow(lambda, 2.0) - 1.0);
        }

        // Active part
        double sigma_act = 0.0;
        double I4_hat = std::pow(lambda / p_curr.lambda_0, 2.0);
        double I4_0   = std::pow(p_curr.lambda_0, 2.0);
        double mu_a_val = get_mu_a(lambda);
        
        // 論文 Eq(75)の最後の項: [1 - (I4_hat - I4_0)/(2*c2)]
        double bracket = 1.0 - (I4_hat - I4_0) / (2.0 * p_curr.c2); 
        sigma_act = 2.0 * p_curr.mu * mu_a_val * alpha * I4_hat * (I4_hat - I4_0) * bracket;

        return sigma_iso + sigma_pass + sigma_act;
    }

    // ヒステリシスを含む現在の応力計算
    double compute_stress(double lambda) {
        double W0_curr = get_W0(lambda);
        double sigma0 = get_sigma0(lambda);

        // Loading process (Primary Path)
        if (lambda >= lambda_max) {
            lambda_max = lambda;
            W0_max = W0_curr;
            return sigma0; 
        } 
        
        // Unloading / Reloading process (Pseudo-elasticity)
        // 論文 Eq. 76: eta = 1 - (1/r) * tanh((W_max - W_curr)/m)
        // Active化中はパラメータが変わるため、現在のパラメータセットにおけるW0_maxを基準にする
        double W0_max_current_params = get_W0(lambda_max); 
        
        double arg = (W0_max_current_params - W0_curr) / p_curr.m;
        double eta = 1.0 - (1.0 / p_curr.r) * std::tanh(arg);

        return eta * sigma0;
    }
};

// 単純なNeo-Hookeanモデル（比較用）
// sigma = mu * (lambda^2 - 1/lambda)
double neo_hookean_stress(double lambda, double mu) {
    return mu * (std::pow(lambda, 2.0) - 1.0 / lambda);
}

int main() {
    std::ofstream file("muscle_comparison.csv");
    file << "Time,Lambda,Alpha,Stress_Model(MPa),Stress_NH(MPa),Phase\n";

    MuscleModel muscle;
    double mu_base = 755.5; // Neo-Hookean用（モデルのベース剛性と一致させる）

    // Simulation Parameters
    double lambda_start = 1.0;
    double lambda_peak = 1.2;
    double lambda_end = 1.0;

    int steps = 100;
    double time = 0.0;
    double dt = 0.1;

    // --- Phase 1: Passive Loading (1.0 -> 1.2) ---
    // Alpha = 0 固定
    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        double lam = lambda_start + t * (lambda_peak - lambda_start);
        
        double s_model = muscle.compute_stress(lam);
        double s_nh = neo_hookean_stress(lam, mu_base);
        
        file << time << "," << lam << "," << 0.0 << "," 
             << s_model / 1e6 << "," << s_nh / 1e6 << ",1_PassiveLoad\n";
        time += dt;
    }

    // --- Phase 2: Unloading with Gradual Activation (1.2 -> 1.0) ---
    // 除荷しながら、途中(例: 1.15)から活性化を開始し、1.05で完了させる
    double activation_start_lam = 1.15;
    double activation_end_lam   = 1.05;

    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        double lam = lambda_peak - t * (lambda_peak - lambda_end);
        
        // 活性化ロジック: Lambdaが特定範囲にある間にAlphaを0->1にする
        double alpha = 0.0;
        if (lam < activation_start_lam && lam > activation_end_lam) {
            // 線形補間: startで0, endで1
            alpha = (activation_start_lam - lam) / (activation_start_lam - activation_end_lam);
        } else if (lam <= activation_end_lam) {
            alpha = 1.0;
        }

        muscle.set_activation(alpha);
        
        double s_model = muscle.compute_stress(lam);
        double s_nh = neo_hookean_stress(lam, mu_base); // NHはalphaの影響を受けない

        file << time << "," << lam << "," << alpha << "," 
             << s_model / 1e6 << "," << s_nh / 1e6 << ",2_UnloadTrans\n";
        time += dt;
    }

    // --- Phase 3: Active Reloading (1.0 -> 1.2) ---
    // Alpha = 1 固定
    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        double lam = lambda_end + t * (lambda_peak - lambda_end);
        
        double s_model = muscle.compute_stress(lam);
        double s_nh = neo_hookean_stress(lam, mu_base);

        file << time << "," << lam << "," << 1.0 << "," 
             << s_model / 1e6 << "," << s_nh / 1e6 << ",3_ActiveReload\n";
        time += dt;
    }

    // --- Phase 4: Active Unloading (1.2 -> 1.0) ---
    // Alpha = 1 固定
    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        double lam = lambda_peak - t * (lambda_peak - lambda_end);
        
        double s_model = muscle.compute_stress(lam);
        double s_nh = neo_hookean_stress(lam, mu_base);

        file << time << "," << lam << "," << 1.0 << "," 
             << s_model / 1e6 << "," << s_nh / 1e6 << ",4_ActiveUnload\n";
        time += dt;
    }

    std::cout << "Simulation finished. Results in muscle_comparison.csv" << std::endl;
    return 0;
}