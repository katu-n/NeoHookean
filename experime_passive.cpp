#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

// 数学定数
const double PI = 3.14159265358979323846;

// =========================================================
// 構成則モデルの基底クラス
// =========================================================
class ConstitutiveModel {
protected:
    double mu; // 基質剛性 (Pa)

public:
    ConstitutiveModel(double mu_val) : mu(mu_val) {}
    virtual ~ConstitutiveModel() {}

    // 応力計算 (Cauchy Stress)
    // lambda_c: 円周方向伸展比, lambda_a: 軸方向伸展比
    // 返り値: {sigma_circumferential, sigma_axial}
    virtual std::pair<double, double> get_stress(double lambda_c, double lambda_a) = 0;
};

// 1. Neo-Hookean モデル (等方性)
class NeoHookean : public ConstitutiveModel {
public:
    NeoHookean(double mu_val) : ConstitutiveModel(mu_val) {}

    std::pair<double, double> get_stress(double lambda_c, double lambda_a) override {
        // 非圧縮性条件: lambda_r = 1 / (lambda_c * lambda_a)
        // Sigma_i = mu * (lambda_i^2 - lambda_r^2)
        
        double lambda_r = 1.0 / (lambda_c * lambda_a);
        double s_c = mu * (pow(lambda_c, 2.0) - pow(lambda_r, 2.0));
        double s_a = mu * (pow(lambda_a, 2.0) - pow(lambda_r, 2.0));
        
        return {s_c, s_a};
    }
};

// 2. Paetsch et al. (多項式型異方性)
// W_fib = mu * mp / 2 * (I4 - 1)^2
class PaetschPassive : public ConstitutiveModel {
private:
    double mp;       // 繊維剛性係数 (無次元)
    double alpha_rad;// 繊維角度 (ラジアン)

public:
    PaetschPassive(double mu_val, double mp_val, double alpha_deg) 
        : ConstitutiveModel(mu_val), mp(mp_val) {
        alpha_rad = alpha_deg * PI / 180.0;
    }

    std::pair<double, double> get_stress(double lambda_c, double lambda_a) override {
        // I4 = lambda_c^2 * cos^2(alpha) + lambda_a^2 * sin^2(alpha)
        double cos_a = cos(alpha_rad);
        double sin_a = sin(alpha_rad);
        double I4 = pow(lambda_c * cos_a, 2.0) + pow(lambda_a * sin_a, 2.0);

        // 基質応力 (Neo-Hookean)
        double lambda_r = 1.0 / (lambda_c * lambda_a);
        double s_c_iso = mu * (pow(lambda_c, 2.0) - pow(lambda_r, 2.0));
        double s_a_iso = mu * (pow(lambda_a, 2.0) - pow(lambda_r, 2.0));

        // 繊維応力
        // dW/dI4 = mu * mp * (I4 - 1)
        // Sigma_fib_i = 2 * (dW/dI4) * (partial I4 / partial C_ii) * lambda_i^2
        //             = 2 * (dW/dI4) * M_i^2 * lambda_i^2
        double dWdI4 = 0.0;
        if (I4 > 1.0) {
            dWdI4 = mu * mp * (I4 - 1.0);
        }

        // 2つの繊維ファミリー (+alpha, -alpha) を考慮する場合、
        // sin^2, cos^2の項は対称なので、単純に2倍ではなく成分として加算される
        // ここでは有効応力として計算
        double s_c_fib = 2.0 * dWdI4 * pow(lambda_c, 2.0) * pow(cos_a, 2.0);
        double s_a_fib = 2.0 * dWdI4 * pow(lambda_a, 2.0) * pow(sin_a, 2.0);

        return {s_c_iso + s_c_fib, s_a_iso + s_a_fib};
    }
};

// 3. Zhao et al. (指数関数型異方性)
// W_fib = k1 / (2*k2) * [exp(k2 * (I4 - 1)^2) - 1]
class ZhaoPassive : public ConstitutiveModel {
private:
    double k1;       // 繊維応力パラメータ (Pa)
    double k2;       // 指数パラメータ (無次元)
    double alpha_rad;// 繊維角度 (ラジアン)

public:
    ZhaoPassive(double mu_val, double k1_val, double k2_val, double alpha_deg) 
        : ConstitutiveModel(mu_val), k1(k1_val), k2(k2_val) {
        alpha_rad = alpha_deg * PI / 180.0;
    }

    std::pair<double, double> get_stress(double lambda_c, double lambda_a) override {
        double cos_a = cos(alpha_rad);
        double sin_a = sin(alpha_rad);
        double I4 = pow(lambda_c * cos_a, 2.0) + pow(lambda_a * sin_a, 2.0);

        // 基質応力
        double lambda_r = 1.0 / (lambda_c * lambda_a);
        double s_c_iso = mu * (pow(lambda_c, 2.0) - pow(lambda_r, 2.0));
        double s_a_iso = mu * (pow(lambda_a, 2.0) - pow(lambda_r, 2.0));

        // 繊維応力
        // dW/dI4 = k1 * (I4 - 1) * exp(k2 * (I4 - 1)^2)
        double dWdI4 = 0.0;
        if (I4 > 1.0) {
            double exponent = k2 * pow(I4 - 1.0, 2.0);
            // オーバーフロー防止のキャップ（物理的にありえない高応力回避）
            if (exponent > 100.0) exponent = 100.0; 
            dWdI4 = k1 * (I4 - 1.0) * exp(exponent);
        }

        double s_c_fib = 2.0 * dWdI4 * pow(lambda_c, 2.0) * pow(cos_a, 2.0);
        double s_a_fib = 2.0 * dWdI4 * pow(lambda_a, 2.0) * pow(sin_a, 2.0);

        return {s_c_iso + s_c_fib, s_a_iso + s_a_fib};
    }
};

int main() {
    std::ofstream file("experiment1_passive.csv");
    file << "Lambda,Neo_Circ,Neo_Axial,Paetsch_Circ,Paetsch_Axial,Zhao_Circ,Zhao_Axial\n";

    // --- パラメータ設定 ---
    // ソース: Zhao et al. (2021) / Shokrani et al. (2024 Corrigendum) Table 4 (Outer Colonic)
    // 単位は kPa なので Pa に変換して入力 (x1000)
    double mu_zhao = 83.10 * 1000.0; 
    double k1_zhao = 573.9 * 1000.0;
    double k2_zhao = 69.75;
    double alpha_deg = 47.58; // ほぼ48度

    // Paetschモデルのパラメータ調整
    // 比較のため、基質剛性はZhaoモデルに合わせる
    // mp (繊維剛性) は Paetsch の元の比率 (mp ~ 86) を参考にしつつ、
    // k1のスケール感に合わせて調整（多項式と指数の比較のため）
    double mu_paetsch = mu_zhao;
    double mp_paetsch = 50.0; // 任意の係数

    // モデルのインスタンス化
    NeoHookean model_neo(mu_zhao);
    PaetschPassive model_paetsch(mu_paetsch, mp_paetsch, alpha_deg);
    ZhaoPassive model_zhao(mu_zhao, k1_zhao, k2_zhao, alpha_deg);

    // --- シミュレーションループ (等二軸伸展) ---
    double lambda_start = 1.0;
    double lambda_end = 1.15; // 15% 伸展
    int steps = 100;

    for (int i = 0; i <= steps; ++i) {
        double lam = lambda_start + (double)i / steps * (lambda_end - lambda_start);
        
        // 等二軸伸展: lambda_c = lambda_a = lam
        auto s_neo = model_neo.get_stress(lam, lam);
        auto s_paetsch = model_paetsch.get_stress(lam, lam);
        auto s_zhao = model_zhao.get_stress(lam, lam);

        // CSV出力 (単位をkPaに変換して保存)
        file << lam << "," 
             << s_neo.first / 1000.0 << "," << s_neo.second / 1000.0 << ","
             << s_paetsch.first / 1000.0 << "," << s_paetsch.second / 1000.0 << ","
             << s_zhao.first / 1000.0 << "," << s_zhao.second / 1000.0 << "\n";
    }

    std::cout << "Simulation complete. Results saved to experiment1_passive.csv" << std::endl;
    return 0;
}