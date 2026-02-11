#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

// =========================================================
// 構成則モデル (Zhao vs Integrated)
// =========================================================
class MuscleModel {
private:
    // --- Zhao Passive Parameters (Outer Rectal) ---
    // ソース: newneoHookean.pdf (Corrigendum Table 4)
    double mu = 6.0 * 1000.0;     // Matrix (Pa) - 柔らかめ
    double k1 = 41.0 * 1000.0;    // Fiber stiffness (Pa)
    double k2 = 14.0;             // Fiber nonlinearity

    // --- Paetsch Active Parameters ---
    // ソース: passive.pdf, miosin.pdf (Cross-bridge stiffness)
    double lambda_0 = 0.7;        // Active自然長 (かなり短縮させて残留応力を出す)
    double S_act_max = 100.0 * 1000.0; // Active最大応力 (Pa)
    double c2 = 0.5;              // Active力の有効範囲 (Overlap shape)

public:
    // Zhaoモデル (Passiveのみ) の応力計算
    double get_zhao_stress(double lambda) {
        // 1. Matrix (Neo-Hookean)
        double sigma_iso = mu * (pow(lambda, 2.0) - 1.0 / lambda);

        // 2. Passive Fiber (Exponential)
        double sigma_fib = 0.0;
        if (lambda > 1.0) {
            double I4 = pow(lambda, 2.0);
            double term = k2 * pow(I4 - 1.0, 2.0);
            // 応力 = 2 * I4 * dW/dI4
            sigma_fib = 2.0 * I4 * k1 * (I4 - 1.0) * std::exp(term);
        }
        return sigma_iso + sigma_fib;
    }

    // 提案モデル (Integrated: Passive + Active) の応力計算
    double get_integrated_stress(double lambda, double alpha) {
        // ベースはZhao Passive
        double sigma_passive = get_zhao_stress(lambda);

        // Active成分 (Paetsch)
        double sigma_active = 0.0;
        if (alpha > 0.0) {
            // 中間配置理論: 現在の長さ lambda が active自然長 lambda_0 に対してどう変形しているか
            // lambda = 1.0 のとき、もし lambda_0 = 0.7 なら、筋肉は「伸ばされた状態」にあるため張力を発揮する
            
            // Overlap function (釣鐘型)
            double overlap = std::exp( -pow(lambda - lambda_0, 2.0) / c2 );
            
            // Active Stress = alpha * Stiffness * Strain_Factor * Overlap
            // ここでは簡易的に Paetsch Eq.75 の挙動を再現
            // (lambda/lambda_0)^2 * (I4_hat - I4_0) の項に相当する力
            double stretch_factor = (lambda / lambda_0); 
            sigma_active = alpha * S_act_max * overlap * (stretch_factor - 1.0); 
            
            // 補正: Active張力は負にはならない (縮みすぎたら力が出ない)
            if (sigma_active < 0) sigma_active = 0.0;
        }

        return sigma_passive + sigma_active;
    }
};

int main() {
    std::ofstream file("comparison_experiment.csv");
    file << "Step,Lambda,Alpha,Stress_Zhao(kPa),Stress_Integrated(kPa),Phase_Label\n";

    MuscleModel model;
    double lambda = 1.0;
    double alpha = 0.0;
    int step_count = 0;

    // --- 実験プロトコル ---
    
    // Phase 1: Passive Loading (1.0 -> 1.15)
    // ZhaoモデルとIntegratedモデル(alpha=0)は一致するはず
    for (int i = 0; i <= 50; ++i) {
        lambda = 1.0 + 0.15 * (double)i / 50.0;
        double s_zhao = model.get_zhao_stress(lambda);
        double s_int = model.get_integrated_stress(lambda, 0.0);
        file << step_count++ << "," << lambda << ",0," << s_zhao/1000.0 << "," << s_int/1000.0 << ",Passive_Load\n";
    }

    // Phase 2: Passive Unloading (1.15 -> 1.0)
    for (int i = 0; i <= 50; ++i) {
        lambda = 1.15 - 0.15 * (double)i / 50.0;
        double s_zhao = model.get_zhao_stress(lambda);
        double s_int = model.get_integrated_stress(lambda, 0.0);
        file << step_count++ << "," << lambda << ",0," << s_zhao/1000.0 << "," << s_int/1000.0 << ",Passive_Unload\n";
    }

    // Phase 3: Activation at Rest (Lambda=1.0, Alpha 0 -> 1)
    // ここで「残留応力」が発生する (垂直な立ち上がり)
    lambda = 1.0;
    for (int i = 0; i <= 20; ++i) {
        alpha = (double)i / 20.0;
        double s_zhao = model.get_zhao_stress(lambda); // 変化なし
        double s_int = model.get_integrated_stress(lambda, alpha); // 急上昇
        file << step_count++ << "," << lambda << "," << alpha << "," << s_zhao/1000.0 << "," << s_int/1000.0 << ",Activation\n";
    }

    // Phase 4: Active Loading (1.0 -> 1.15)
    // ここで「剛性増加」が見える (Zhaoよりも傾きが急になる)
    alpha = 1.0;
    for (int i = 0; i <= 50; ++i) {
        lambda = 1.0 + 0.15 * (double)i / 50.0;
        double s_zhao = model.get_zhao_stress(lambda);
        double s_int = model.get_integrated_stress(lambda, 1.0);
        file << step_count++ << "," << lambda << ",1," << s_zhao/1000.0 << "," << s_int/1000.0 << ",Active_Load\n";
    }

    // Phase 5: Active Unloading (1.15 -> 1.0)
    // 戻ってきても応力はゼロにならない
    for (int i = 0; i <= 50; ++i) {
        lambda = 1.15 - 0.15 * (double)i / 50.0;
        double s_zhao = model.get_zhao_stress(lambda);
        double s_int = model.get_integrated_stress(lambda, 1.0);
        file << step_count++ << "," << lambda << ",1," << s_zhao/1000.0 << "," << s_int/1000.0 << ",Active_Unload\n";
    }

    std::cout << "Data generated: comparison_experiment.csv" << std::endl;
    return 0;
}