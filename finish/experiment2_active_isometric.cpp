#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

class MuscleModel {
private:
    // --- Material Parameters ---
    double mu;      // Matrix shear modulus (Pa)
    
    // Passive Fiber (Zhao et al. Outer Colonic)
    double k1;      // Fiber stress parameter (Pa)
    double k2;      // Exponential coefficient (Reduced for visibility)

    // Active Contractile (Paetsch et al. modified)
    double lambda_0; // Optimal active length (relative to reference)
    double S_act_max; // Max active stress parameter (Pa)
    double c2;       // Width of the active force-length curve

public:
    MuscleModel() {
        // 1. Passive Parameters (Zhao et al. / Shokrani et al. 2024)
        mu = 10.0 * 1000.0;   // Matrixは柔らかく設定 (10 kPa)
        k1 = 50.0 * 1000.0;   // Passive係数 (50 kPa)
        k2 = 20.0;            // 指数係数 (少し抑えめに設定してカーブを見やすくする)

        // 2. Active Parameters (Physiological tuning)
        lambda_0 = 1.0;       // 自然長でピークが来るように設定（または0.95）
        S_act_max = 100.0 * 1000.0; // 最大Active応力を100kPaに設定
        c2 = 0.8;             // 釣鐘型の幅を広げる (小さいとすぐゼロになる)
    }

    double get_stress(double lambda, double alpha) {
        // 1. Matrix (Neo-Hookean)
        // 圧縮側(lambda < 1)でも計算する
        double sigma_iso = mu * (pow(lambda, 2.0) - 1.0/lambda);

        // 2. Passive Fiber (Zhao Exponential)
        double sigma_pass = 0.0;
        double I4 = pow(lambda, 2.0);
        // 繊維は引張(lambda > 1.0)のみ負担すると仮定
        if (lambda > 1.0) {
            double term = k2 * pow(I4 - 1.0, 2.0);
            if (term > 50.0) term = 50.0; // Overflow guard
            sigma_pass = 2.0 * I4 * k1 * (I4 - 1.0) * std::exp(term);
        }

        // 3. Active Fiber (Paetsch-like Gaussian/Bell shape)
        double sigma_act = 0.0;
        if (alpha > 0.0) {
            // 滑り説に基づく簡易的なActive張力曲線（正規分布型）
            // Peak at lambda_0, width controlled by c2
            // S_act = alpha * S_max * exp( - (lambda - lambda_0)^2 / (2 * width^2) )
            
            // Paetschの不変量形式を少し簡略化して挙動を安定させます
            double term = -pow(lambda - lambda_0, 2.0) / (2.0 * pow(c2/2.0, 2.0));
            double overlap_factor = std::exp(term);
            
            // 伸展比そのものによる幾何学的寄与も含める場合:
            sigma_act = alpha * S_act_max * overlap_factor * (lambda / lambda_0);
        }

        return sigma_iso + sigma_pass + sigma_act;
    }
};

int main() {
    std::ofstream file("experiment2_revised.csv");
    file << "Lambda,Passive_Stress(kPa),Total_Stress(kPa),Active_Stress(kPa)\n";

    MuscleModel model;

    // Simulation Range: 0.7 to 1.4
    double lam_start = 0.7;
    double lam_end = 1.4;
    int steps = 100;

    for (int i = 0; i <= steps; ++i) {
        double lam = lam_start + (double)i / steps * (lam_end - lam_start);

        double s_pass = model.get_stress(lam, 0.0); // alpha=0
        double s_total = model.get_stress(lam, 1.0); // alpha=1
        double s_active = s_total - s_pass;

        file << lam << "," 
             << s_pass / 1000.0 << "," 
             << s_total / 1000.0 << "," 
             << s_active / 1000.0 << "\n";
    }

    std::cout << "Done. Results in experiment2_revised.csv" << std::endl;
    return 0;
}