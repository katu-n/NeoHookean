#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// ==============================================================================
// パラメータ設定 (Paetsch et al., 2012 / Table 1 に基づく)
// ==============================================================================
const double mu = 0.7555;      // kPa
const double m_p = 86.13;      // -
const double lambda_0 = 0.95;  // -
const double c1 = 470.0;       // -
const double c2 = 0.95;        // -

const double r_pass = 1.05;    // -
const double m_pass = 2.5;     // kPa
const double r_act = 1.55;     // -
const double m_act = 30.0;     // kPa

// ==============================================================================
// 論文の1次元理論式に基づく解析解(Analytical Solution)の計算関数 [Eq. 75, 77]
// ==============================================================================
double get_ma(double lambda) {
    double term = std::pow(lambda / lambda_0, 2.0) - std::pow(lambda_0, 2.0);
    return c1 * std::exp(-term / c2);
}

double get_W0(double lambda, double alpha) {
    double term_iso = lambda * lambda + 2.0 / lambda - 3.0;
    double term_pass = m_p * (1.0 - alpha) * std::pow(lambda * lambda - 1.0, 2.0);
    double term_act_base = std::pow(lambda / lambda_0, 2.0) - std::pow(lambda_0, 2.0);
    double term_act = alpha * get_ma(lambda) * std::pow(term_act_base, 2.0);
    return (mu / 2.0) * (term_iso + term_pass + term_act);
}

double get_sigma0(double lambda, double alpha) {
    double term_iso = mu * (lambda * lambda - 1.0 / lambda);
    double term_pass = 2.0 * mu * m_p * (1.0 - alpha) * lambda * lambda * (lambda * lambda - 1.0);
    double term_act_base = std::pow(lambda / lambda_0, 2.0) - std::pow(lambda_0, 2.0);
    double term_act = 2.0 * mu * get_ma(lambda) * alpha * std::pow(lambda / lambda_0, 2.0) * term_act_base * (1.0 - term_act_base / (2.0 * c2));
    return term_iso + term_pass + term_act;
}

double get_Z(double lambda, double alpha, double lambda_max) {
    double W0_max = get_W0(lambda_max, alpha);
    double W0_cur = get_W0(lambda, alpha);
    double r = (1.0 - alpha) * r_pass + alpha * r_act;
    double m = (1.0 - alpha) * m_pass + alpha * m_act;
    return 1.0 - (1.0 / r) * std::tanh((W0_max - W0_cur) / m);
}

// ==============================================================================
// RMSE計算関数
// ==============================================================================
double calculate_rmse(const std::vector<double>& sim, const std::vector<double>& analytical) {
    double sum_sq_error = 0.0;
    int n = sim.size();
    for (int i = 0; i < n; ++i) {
        double error = sim[i] - analytical[i];
        sum_sq_error += error * error;
    }
    return std::sqrt(sum_sq_error / n);
}

// ==============================================================================
// メイン処理: 解析解の算出とIGAシミュレーション結果との比較
// ==============================================================================
int main() {
    double lambda_max = 1.15;
    
    // 検証ポイント (Stretch: lambda)
    // 1. 負荷時(Loading): 受動状態 alpha = 0.0
    std::vector<double> load_lambdas = {1.0, 1.05, 1.10, 1.15};
    std::vector<double> load_alpha   = {0.0, 0.0,  0.0,  0.0};
    
    // 2. 除荷時(Unloading): 能動状態への遷移 alpha: 0.0 -> 1.0
    std::vector<double> unload_lambdas = {1.15, 1.10, 1.05, 1.00};
    std::vector<double> unload_alpha   = {0.0,  0.0,  0.5,  1.0};

    // --------------------------------------------------------------------------
    // あなたのIGAプログラムから得られたシミュレーション結果 (ダミー値)
    // ※ ここにIGAソルバーから出力されたCauchy応力を入力してください
    // --------------------------------------------------------------------------
    std::vector<double> iga_sim_load   = {0.0, 15.3, 40.2, 85.0}; // 負荷時のIGA応力(kPa)
    std::vector<double> iga_sim_unload = {85.0, 30.1, 80.5, 115.2}; // 除荷・活性化時のIGA応力(kPa)

    std::vector<double> analyt_load(load_lambdas.size());
    std::vector<double> analyt_unload(unload_lambdas.size());

    // 解析解の計算と出力
    std::cout << "--- 【負荷フェーズ】 受動状態 (Loading, Alpha=0) ---" << std::endl;
    std::cout << std::setw(10) << "Stretch" << std::setw(15) << "IGA_Sim(kPa)" << std::setw(15) << "Analytical" << std::endl;
    for (size_t i = 0; i < load_lambdas.size(); ++i) {
        analyt_load[i] = get_sigma0(load_lambdas[i], load_alpha[i]); // Z=1
        std::cout << std::fixed << std::setprecision(3)
                  << std::setw(10) << load_lambdas[i]
                  << std::setw(15) << iga_sim_load[i] 
                  << std::setw(15) << analyt_load[i] << std::endl;
    }

    std::cout << "\n--- 【除荷フェーズ】 能動状態への遷移 (Unloading, Alpha 0->1) ---" << std::endl;
    std::cout << std::setw(10) << "Stretch" << std::setw(10) << "Alpha" << std::setw(15) << "IGA_Sim(kPa)" << std::setw(15) << "Analytical" << std::endl;
    for (size_t i = 0; i < unload_lambdas.size(); ++i) {
        double sigma_0 = get_sigma0(unload_lambdas[i], unload_alpha[i]);
        double Z = get_Z(unload_lambdas[i], unload_alpha[i], lambda_max);
        analyt_unload[i] = Z * sigma_0;
        
        std::cout << std::fixed << std::setprecision(3)
                  << std::setw(10) << unload_lambdas[i] << std::setw(10) << unload_alpha[i]
                  << std::setw(15) << iga_sim_unload[i] 
                  << std::setw(15) << analyt_unload[i] << std::endl;
    }

    // RMSEの計算と評価
    double rmse_load = calculate_rmse(iga_sim_load, analyt_load);
    double rmse_unload = calculate_rmse(iga_sim_unload, analyt_unload);

    std::cout << "\n================ 定量評価結果 (RMSE) ================" << std::endl;
    std::cout << "負荷フェーズ (Loading) RMSE   : " << rmse_load << " kPa" << std::endl;
    std::cout << "除荷フェーズ (Unloading) RMSE : " << rmse_unload << " kPa" << std::endl;
    std::cout << "=====================================================" << std::endl;

    return 0;
}