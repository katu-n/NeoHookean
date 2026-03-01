#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

// --- 定数の定義 ---
const double PI = std::acos(-1.0);

const int voigt_map[3][3] = {
    {0, 3, 5},
    {3, 1, 4},
    {5, 4, 2}
};

// ==============================================================================
// IGAマテリアルルーチン: S と D を計算する関数
// ==============================================================================
void compute_S_and_D(const double C[3][3], double mu, double lambda_lame, 
                     double k1, double k2, const double M[3], 
                     double S[6], double D[6][6]) {
                     
    double detC = C[0][0]*(C[1][1]*C[2][2] - C[1][2]*C[2][1])
                - C[0][1]*(C[1][0]*C[2][2] - C[1][2]*C[2][0])
                + C[0][2]*(C[1][0]*C[2][1] - C[1][1]*C[2][0]);
                
    double J = std::sqrt(detC);
    double lnJ = std::log(J);
    
    double invC[3][3];
    invC[0][0] = (C[1][1]*C[2][2] - C[1][2]*C[2][1]) / detC;
    invC[0][1] = (C[0][2]*C[2][1] - C[0][1]*C[2][2]) / detC;
    invC[0][2] = (C[0][1]*C[1][2] - C[0][2]*C[1][1]) / detC;
    invC[1][0] = (C[1][2]*C[2][0] - C[1][0]*C[2][2]) / detC;
    invC[1][1] = (C[0][0]*C[2][2] - C[0][2]*C[2][0]) / detC;
    invC[1][2] = (C[0][2]*C[1][0] - C[0][0]*C[1][2]) / detC;
    invC[2][0] = (C[1][0]*C[2][1] - C[1][1]*C[2][0]) / detC;
    invC[2][1] = (C[0][1]*C[2][0] - C[0][0]*C[2][1]) / detC;
    invC[2][2] = (C[0][0]*C[1][1] - C[0][1]*C[1][0]) / detC;

    double I4 = 0.0;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            I4 += M[i] * C[i][j] * M[j];
        }
    }

    double fiber_stress_coeff = 0.0;
    double fiber_stiff_coeff  = 0.0;
    if (I4 > 1.0) { 
        double term = I4 - 1.0;
        double exp_term = std::exp(k2 * term * term);
        fiber_stress_coeff = 2.0 * k1 * term * exp_term;
        fiber_stiff_coeff  = 4.0 * k1 * (1.0 + 2.0 * k2 * term * term) * exp_term;
    }

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            int row = voigt_map[i][j];
            double I_ij = (i == j) ? 1.0 : 0.0;
            
            S[row] = mu * I_ij + (lambda_lame * lnJ - mu) * invC[i][j] + fiber_stress_coeff * (M[i] * M[j]);
                        
            for (int k = 0; k < 3; ++k) {
                for (int l = 0; l < 3; ++l) {
                    int col = voigt_map[k][l];
                    double I_Cinv = 0.5 * (invC[i][k] * invC[j][l] + invC[i][l] * invC[j][k]);
                    D[row][col] = lambda_lame * invC[i][j] * invC[k][l]
                                  + 2.0 * (mu - lambda_lame * lnJ) * I_Cinv
                                  + fiber_stiff_coeff * (M[i] * M[j] * M[k] * M[l]);
                }
            }
        }
    }
}

// ==============================================================================
// 論文(Zhao et al.)に基づく、完全非圧縮の理論解の計算
// ==============================================================================
double calculate_analytical_stress(double lambda, int direction, double mu, double k1, double k2, double alpha_deg) {
    double alpha_rad = alpha_deg * PI / 180.0;
    double lambda_r = 1.0 / (lambda * lambda); 
    double I4 = (lambda * lambda * std::cos(alpha_rad) * std::cos(alpha_rad)) + 
                (lambda * lambda * std::sin(alpha_rad) * std::sin(alpha_rad));
                
    double Psi4 = 0.0;
    if (I4 > 1.0) { 
        Psi4 = k1 * (I4 - 1.0) * std::exp(k2 * std::pow(I4 - 1.0, 2.0));
    }
    double stress_matrix = mu * (lambda * lambda - lambda_r * lambda_r);
    
    if (direction == 0) { // 円周方向 (theta)
        return stress_matrix + 2.0 * lambda * lambda * std::cos(alpha_rad) * std::cos(alpha_rad) * Psi4;
    } else {              // 軸方向 (z)
        return stress_matrix + 2.0 * lambda * lambda * std::sin(alpha_rad) * std::sin(alpha_rad) * Psi4;
    }
}

// ==============================================================================
// メイン処理
// ==============================================================================
// int main() {
//     double mu = 83.10;
//     double k1 = 573.9;
//     double k2 = 69.75;
//     double alpha_deg = 47.58;
//     double lambda_lame = mu * 10000.0; 
//     double alpha_rad = alpha_deg * PI / 180.0;
//     double M_vec[3] = {std::cos(alpha_rad), std::sin(alpha_rad), 0.0};

//     std::vector<double> stretches;
//     for (double s = 1.0; s <= 1.1001; s += 0.005) {
//         stretches.push_back(s);
//     }
    
//     std::ofstream ofs("stress_results.csv");
//     if (!ofs) {
//         std::cerr << "エラー: CSVファイルを開けませんでした。" << std::endl;
//         return 1;
//     }
//     // ヘッダーに NeoHookean を追加
//     ofs << "Stretch,Sim_Theta,Analyt_Theta,Sim_Axial,Analyt_Axial,NeoHookean\n";

//     for (double lambda : stretches) {
//         double lambda_3 = 1.0 / (lambda * lambda); 
//         double S[6];
//         double D[6][6];
        
//         for (int iter = 0; iter < 10; ++iter) {
//             double C[3][3] = {{0.0}};
//             C[0][0] = lambda * lambda;
//             C[1][1] = lambda * lambda;
//             C[2][2] = lambda_3 * lambda_3;
            
//             compute_S_and_D(C, mu, lambda_lame, k1, k2, M_vec, S, D);
            
//             double S33 = S[2]; 
//             if (std::abs(S33) < 1e-6) break; 
            
//             double dS33_dlambda3 = D[2][2] * lambda_3; 
//             lambda_3 -= S33 / dS33_dlambda3;
//         }
        
//         double J = lambda * lambda * lambda_3;
//         double sim_sigma_theta = (lambda * S[0] * lambda) / J;
//         double sim_sigma_z     = (lambda * S[1] * lambda) / J;

//         double analyt_sigma_theta = calculate_analytical_stress(lambda, 0, mu, k1, k2, alpha_deg);
//         double analyt_sigma_z     = calculate_analytical_stress(lambda, 1, mu, k1, k2, alpha_deg);

//         // Neo-Hookeanモデルの理論解（完全非圧縮）を計算
//         double lambda_r = 1.0 / (lambda * lambda);
//         double neo_hookean_stress = mu * (lambda * lambda - lambda_r * lambda_r);

//         ofs << std::fixed << std::setprecision(5)
//             << lambda << "," 
//             << sim_sigma_theta << "," << analyt_sigma_theta << ","
//             << sim_sigma_z << "," << analyt_sigma_z << ","
//             << neo_hookean_stress << "\n";
            
//         std::cout << "Stretch: " << lambda << " Done." << std::endl;
//     }

//     ofs.close();
//     std::cout << "CSVファイルの出力が完了しました: stress_results.csv" << std::endl;

//     return 0;
// }
// ==============================================================================
// メイン処理
// ==============================================================================
int main() {
    double mu = 83.10;
    double k1 = 573.9;
    double k2 = 69.75;
    double alpha_deg = 47.58;
    double lambda_lame = mu * 10000.0; 
    double alpha_rad = alpha_deg * PI / 180.0;
    double M_vec[3] = {std::cos(alpha_rad), std::sin(alpha_rad), 0.0};

    std::vector<double> stretches;
    for (double s = 1.0; s <= 1.1001; s += 0.005) {
        stretches.push_back(s);
    }
    
    // --- RMSE計算用の変数（追加） ---
    double sum_sq_error_theta = 0.0;
    double sum_sq_error_z = 0.0;

    std::ofstream ofs("stress_results.csv");
    if (!ofs) {
        std::cerr << "エラー: CSVファイルを開けませんでした。" << std::endl;
        return 1;
    }
    ofs << "Stretch,Sim_Theta,Analyt_Theta,Sim_Axial,Analyt_Axial,NeoHookean\n";

    for (double lambda : stretches) {
        double lambda_3 = 1.0 / (lambda * lambda); 
        double S[6];
        double D[6][6];
        
        for (int iter = 0; iter < 10; ++iter) {
            double C[3][3] = {{0.0}};
            C[0][0] = lambda * lambda;
            C[1][1] = lambda * lambda;
            C[2][2] = lambda_3 * lambda_3;
            
            compute_S_and_D(C, mu, lambda_lame, k1, k2, M_vec, S, D);
            
            double S33 = S[2]; 
            if (std::abs(S33) < 1e-6) break; 
            
            double dS33_dlambda3 = D[2][2] * lambda_3; 
            lambda_3 -= S33 / dS33_dlambda3;
        }
        
        double J = lambda * lambda * lambda_3;
        double sim_sigma_theta = (lambda * S[0] * lambda) / J;
        double sim_sigma_z     = (lambda * S[1] * lambda) / J;

        double analyt_sigma_theta = calculate_analytical_stress(lambda, 0, mu, k1, k2, alpha_deg);
        double analyt_sigma_z     = calculate_analytical_stress(lambda, 1, mu, k1, k2, alpha_deg);

        // RMSE用に二乗誤差を蓄積（追加）
        sum_sq_error_theta += std::pow(sim_sigma_theta - analyt_sigma_theta, 2.0);
        sum_sq_error_z     += std::pow(sim_sigma_z - analyt_sigma_z, 2.0);

        // Neo-Hookeanモデルの理論解
        double lambda_r = 1.0 / (lambda * lambda);
        double neo_hookean_stress = mu * (lambda * lambda - lambda_r * lambda_r);

        ofs << std::fixed << std::setprecision(5)
            << lambda << "," 
            << sim_sigma_theta << "," << analyt_sigma_theta << ","
            << sim_sigma_z << "," << analyt_sigma_z << ","
            << neo_hookean_stress << "\n";
            
        std::cout << "Stretch: " << lambda << " Done." << std::endl;
    }

    ofs.close();
    std::cout << "\nCSVファイルの出力が完了しました: stress_results.csv" << std::endl;

    // --- RMSEの計算とターミナル出力（追加） ---
    double rmse_theta = std::sqrt(sum_sq_error_theta / stretches.size());
    double rmse_z     = std::sqrt(sum_sq_error_z / stretches.size());

    std::cout << "\n--- 定量評価結果 (IGA vs 理論解) ---" << std::endl;
    std::cout << std::fixed << std::setprecision(6); // 誤差が小さい可能性があるので桁数を増やす
    std::cout << "円周方向 (Circumferential) RMSE : " << rmse_theta << " kPa" << std::endl;
    std::cout << "軸方向   (Axial)          RMSE : " << rmse_z     << " kPa" << std::endl;

    return 0;
}