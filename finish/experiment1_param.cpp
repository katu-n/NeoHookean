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

    virtual std::pair<double, double> get_stress(double lambda_c, double lambda_a) = 0;
};

// 1. Neo-Hookean モデル (等方性)
class NeoHookean : public ConstitutiveModel {
public:
    NeoHookean(double mu_val) : ConstitutiveModel(mu_val) {}

    std::pair<double, double> get_stress(double lambda_c, double lambda_a) override {
        double lambda_r = 1.0 / (lambda_c * lambda_a);
        double s_c = mu * (pow(lambda_c, 2.0) - pow(lambda_r, 2.0));
        double s_a = mu * (pow(lambda_a, 2.0) - pow(lambda_r, 2.0));
        return {s_c, s_a};
    }
};

// 2. Paetsch et al. (多項式型異方性)
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
        double cos_a = cos(alpha_rad);
        double sin_a = sin(alpha_rad);
        // I4は繊維方向の伸展の二乗
        double I4 = pow(lambda_c * cos_a, 2.0) + pow(lambda_a * sin_a, 2.0);

        // 基質
        double lambda_r = 1.0 / (lambda_c * lambda_a);
        double s_c_iso = mu * (pow(lambda_c, 2.0) - pow(lambda_r, 2.0));
        double s_a_iso = mu * (pow(lambda_a, 2.0) - pow(lambda_r, 2.0));

        // 繊維 (多項式型)
        double dWdI4 = 0.0;
        if (I4 > 1.0) {
            dWdI4 = mu * mp * (I4 - 1.0);
        }

        // 2倍するのを忘れない (2つの繊維群あるいは微分係数の2)
        double s_c_fib = 2.0 * dWdI4 * pow(lambda_c, 2.0) * pow(cos_a, 2.0);
        double s_a_fib = 2.0 * dWdI4 * pow(lambda_a, 2.0) * pow(sin_a, 2.0);

        return {s_c_iso + s_c_fib, s_a_iso + s_a_fib};
    }
};

// 3. Zhao et al. (指数関数型異方性)
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

        // 基質
        double lambda_r = 1.0 / (lambda_c * lambda_a);
        double s_c_iso = mu * (pow(lambda_c, 2.0) - pow(lambda_r, 2.0));
        double s_a_iso = mu * (pow(lambda_a, 2.0) - pow(lambda_r, 2.0));

        // 繊維 (指数関数型)
        double dWdI4 = 0.0;
        if (I4 > 1.0) {
            double exponent = k2 * pow(I4 - 1.0, 2.0);
            if (exponent > 50.0) exponent = 50.0; 
            dWdI4 = k1 * (I4 - 1.0) * exp(exponent);
        }

        double s_c_fib = 2.0 * dWdI4 * pow(lambda_c, 2.0) * pow(cos_a, 2.0);
        double s_a_fib = 2.0 * dWdI4 * pow(lambda_a, 2.0) * pow(sin_a, 2.0);

        return {s_c_iso + s_c_fib, s_a_iso + s_a_fib};
    }
};

int main() {
    std::ofstream file("experiment1_passive_adjusted.csv");
    file << "Lambda,Neo_Circ,Neo_Axial,Paetsch_Circ,Paetsch_Axial,Zhao_Circ,Zhao_Axial\n";

    // --- パラメータ調整 (比較のためにスケールを合わせる) ---
    
    // 1. 共通の基質剛性
    // 軟組織として一般的な値 (20kPa程度)
    double mu_val = 20.0 * 1000.0; 
    double alpha_deg = 47.58; 

    // 2. Zhaoパラメータの調整
    // k2 (指数部) が大きすぎると垂直な壁になるため、緩やかにする
    // k1 (係数部) も基質剛性と同程度に設定
    double k1_zhao = 10.0 * 1000.0; // 10 kPa
    double k2_zhao = 15.0;          // 69.75 -> 15.0 に変更（これで程よいJカーブになる）

    // 3. Paetschパラメータの調整
    // Zhaoの立ち上がり初期と重なる程度にする
    double mp_paetsch = 20.0;       // 無次元係数

    // モデル作成
    NeoHookean model_neo(mu_val);
    PaetschPassive model_paetsch(mu_val, mp_paetsch, alpha_deg);
    ZhaoPassive model_zhao(mu_val, k1_zhao, k2_zhao, alpha_deg);

    // --- シミュレーション (等二軸伸展) ---
    double lambda_start = 1.0;
    double lambda_end = 1.15; // 15% 伸展
    int steps = 100;

    for (int i = 0; i <= steps; ++i) {
        double lam = lambda_start + (double)i / steps * (lambda_end - lambda_start);
        
        // 等二軸伸展
        auto s_neo = model_neo.get_stress(lam, lam);
        auto s_paetsch = model_paetsch.get_stress(lam, lam);
        auto s_zhao = model_zhao.get_stress(lam, lam);

        // kPa 単位で出力
        file << lam << "," 
             << s_neo.first / 1000.0 << "," << s_neo.second / 1000.0 << ","
             << s_paetsch.first / 1000.0 << "," << s_paetsch.second / 1000.0 << ","
             << s_zhao.first / 1000.0 << "," << s_zhao.second / 1000.0 << "\n";
    }

    std::cout << "Simulation complete. Results saved to experiment1_passive_adjusted.csv" << std::endl;
    return 0;
}