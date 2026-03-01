#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// --- Zhao et al. (2021) 構成則パラメータ構造体 ---
struct ZhaoParameters {
    double mu;     // kPa (マトリックスのせん断弾性率)
    double k1;     // kPa (繊維の初期剛性)
    double k2;     // 無次元 (繊維の非線形性)
    double alpha;  // 度 (繊維の配向角)
};

// --- 解析解（理論解）を計算するクラス ---
class ZhaoAnalyticalSolution {
private:
    ZhaoParameters params;
    double alpha_rad;

public:
    ZhaoAnalyticalSolution(ZhaoParameters p) : params(p) {
        // 角度をラジアンに変換
        alpha_rad = params.alpha * M_PI / 180.0;
    }

    // 等二軸伸展 (lambda_theta = lambda_z = lambda) 時のCauchy応力を計算
    // direction: 0 = circumferential (円周方向), 1 = axial (軸/縦方向)
    double calculate_cauchy_stress(double lambda, int direction) {
        // 非圧縮性条件より厚み方向の伸展比を計算
        double lambda_r = 1.0 / (lambda * lambda); 

        // 不変量 I4 の計算 (Eq. 6)
        double I4 = (lambda * lambda * std::cos(alpha_rad) * std::cos(alpha_rad)) + 
                    (lambda * lambda * std::sin(alpha_rad) * std::sin(alpha_rad));

        // 繊維のひずみエネルギーの微分 (Psi_4)
        double Psi4 = 0.0;
        if (I4 > 1.0) { // 繊維は引っ張られた時のみ抵抗を発揮する (Heaviside step function)
            Psi4 = params.k1 * (I4 - 1.0) * std::exp(params.k2 * std::pow(I4 - 1.0, 2.0));
        }

        // Cauchy応力の計算 (Eq. 7, 8 をベースに静水圧項を処理した標準形)
        double stress_matrix = params.mu * (lambda * lambda - lambda_r * lambda_r);
        
        if (direction == 0) {
            // 円周方向 (Circumferential: theta)
            double stress_fiber = 2.0 * lambda * lambda * std::cos(alpha_rad) * std::cos(alpha_rad) * Psi4;
            return stress_matrix + stress_fiber;
        } else {
            // 軸方向 (Axial: z)
            double stress_fiber = 2.0 * lambda * lambda * std::sin(alpha_rad) * std::sin(alpha_rad) * Psi4;
            return stress_matrix + stress_fiber;
        }
    }
};

// --- RMSE (二乗平均平方根誤差) を計算する関数 ---
double calculate_rmse(const std::vector<double>& sim_data, const std::vector<double>& analytical_data) {
    if (sim_data.size() != analytical_data.size() || sim_data.empty()) {
        std::cerr << "エラー: データサイズが一致しません。" << std::endl;
        return -1.0;
    }

    double sum_sq_error = 0.0;
    int n = sim_data.size();
    for (int i = 0; i < n; ++i) {
        double error = sim_data[i] - analytical_data[i];
        sum_sq_error += error * error;
    }

    return std::sqrt(sum_sq_error / n);
}

int main() {
    // 1. Zhao et al. (2021) Table 4に基づくパラメータ設定 (例: Outer Colonic)
    ZhaoParameters params = {83.10, 573.9, 69.75, 47.58};
    ZhaoAnalyticalSolution analytical_model(params);

    // 2. 比較する伸展比(Stretch ratio)のポイント
    std::vector<double> stretches = {1.0, 1.02, 1.04, 1.06, 1.08, 1.10};
    int n_points = stretches.size();

    // 3. シミュレーション結果 (ダミーデータ。実際のIGA出力値に置き換えてください)
    // ここでは、解析解から少しだけ誤差を含ませたダミー値を作成しています
    std::vector<double> sim_stress_theta = {0.0, 15.2, 45.1, 95.3, 175.5, 305.2};
    std::vector<double> sim_stress_z     = {0.0, 16.5, 48.0, 102.1, 185.0, 320.1};

    std::vector<double> analytical_stress_theta(n_points);
    std::vector<double> analytical_stress_z(n_points);

    // 4. 解析解の計算と結果の表示
    std::cout << "--- 応力比較 (kPa) ---" << std::endl;
    std::cout << std::setw(10) << "Stretch" 
              << std::setw(15) << "Sim_Theta" << std::setw(15) << "Analyt_Theta" 
              << std::setw(15) << "Sim_Axial" << std::setw(15) << "Analyt_Axial" << std::endl;

    for (int i = 0; i < n_points; ++i) {
        double lambda = stretches[i];
        analytical_stress_theta[i] = analytical_model.calculate_cauchy_stress(lambda, 0);
        analytical_stress_z[i]     = analytical_model.calculate_cauchy_stress(lambda, 1);

        std::cout << std::fixed << std::setprecision(3)
                  << std::setw(10) << lambda
                  << std::setw(15) << sim_stress_theta[i] << std::setw(15) << analytical_stress_theta[i]
                  << std::setw(15) << sim_stress_z[i]     << std::setw(15) << analytical_stress_z[i] << std::endl;
    }

    // 5. RMSEの計算
    double rmse_theta = calculate_rmse(sim_stress_theta, analytical_stress_theta);
    double rmse_z     = calculate_rmse(sim_stress_z, analytical_stress_z);

    std::cout << "\n--- 定量評価結果 (RMSE) ---" << std::endl;
    std::cout << "円周方向 (Circumferential) RMSE : " << rmse_theta << " kPa" << std::endl;
    std::cout << "軸方向   (Axial)           RMSE : " << rmse_z     << " kPa" << std::endl;

    // 総合的な評価（RMSEが最大応力の何%に相当するか等で議論可能）
    return 0;
}