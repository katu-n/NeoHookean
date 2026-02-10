#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

// 数学定数
const double E_NUM = 2.71828182845904523536;

class IntegratedMuscleModel {
public:
    enum PassiveType { POLY_PAETSCH, EXP_ZHAO };
    enum ActiveMechanism { PLASTICITY_ONLY, BIOPHYSICAL_STIMULATION };

private:
    double mu;      // 基質剛性 (Pa)
    // Paetsch (Passive) 用
    double mu_p;    // Passive繊維剛性係数
    // Zhao (Passive) 用
    double k1;      // 繊維剛性パラメータ (Pa)
    double k2;      // 指数パラメータ
    // Activeパラメータ
    double lambda_0; // 自然長の短縮率 (Active時)
    double mu_a_base;// Active剛性の基準値
    double c2;       // 剛性減衰パラメータ

    PassiveType p_type;
    ActiveMechanism a_mech;

public:
    IntegratedMuscleModel(PassiveType pt, ActiveMechanism am) : p_type(pt), a_mech(am) {
        // --- パラメータ設定 (論文値を参考) ---
        mu = 755.5;       // Paetsch Table 1

        // Paetsch Passive
        mu_p = 86.13;

        // Zhao Passive (Colonic Outer)
        // 比較のためPaetschの剛性感にスケールを合わせる調整が必要
        k1 = 2000.0;      
        k2 = 20.0;        

        // Active Parameters
        lambda_0 = 0.95;  // 活性化で5%
        mu_a_base = 470.0; // c1
        c2 = 0.95;
    }

    // -----------------------------------------------------
    // 応力計算 (Cauchy Stress Sigma_11)
    // alpha: 活性化度 (0.0=Passive, 1.0=Max Active)
    // lambda: 現在の伸展比
    // -----------------------------------------------------
    double calculate_stress(double lambda, double alpha) {
        double I4 = lambda * lambda;
        
        // iso-Neo-Hookean
        double sigma_iso = mu * (lambda * lambda - 1.0 / lambda);

        // 2. Passive繊維項 (モデル切り替え)
        double sigma_pass = 0.0;
        if (lambda > 1.0) {
            if (p_type == POLY_PAETSCH) {
                // Paetsch
                sigma_pass = 2.0 * mu * mu_p * I4 * (I4 - 1.0);
            } else if (p_type == EXP_ZHAO) {
                // Zhao
                double term = k2 * (I4 - 1.0) * (I4 - 1.0);
                sigma_pass = 2.0 * k1 * I4 * (I4 - 1.0) * std::exp(term);
            }
        }

        // 3. Active項 (メカニズム切り替え)
        double sigma_act = 0.0;
        
        // Active不変量 (中間配置基準)
        double I4_hat = (lambda * lambda) / (lambda_0 * lambda_0);
        double I4_0 = lambda_0 * lambda_0;

        if (alpha > 0.0) {
            double current_mu_a = 0.0;
            double bracket_term = 1.0;

            if (a_mech == PLASTICITY_ONLY) {
                // 自発的塑性変形:剛性(mu_a)は一定で、自然長(lambda_0)だけが変わるモデル
                current_mu_a = mu_a_base; 
                // 微分項(bracket)は1.0 (剛性がI4に依存しないため)
                bracket_term = 1.0; 
            } 
            else if (a_mech == BIOPHYSICAL_STIMULATION) {
                // 電気刺激・クロスブリッジとして捉える (Paetsch Original)
                // 伸展に応じて剛性が変化する (フィラメントの重なり減少)
                double exponent = -(I4_hat - I4_0) / c2;
                current_mu_a = mu_a_base * std::exp(exponent);
                // 剛性の変化率由来の項が含まれる
                bracket_term = 1.0 - (I4_hat - I4_0) / c2;
            }

            // Active応力計算 (共通形式)
            sigma_act = 2.0 * mu * current_mu_a * I4_hat * (I4_hat - I4_0) * bracket_term;
        }

        return sigma_iso + (1.0 - alpha) * sigma_pass + alpha * sigma_act;
    }
};

int main() {
    std::ofstream file("muscle_study.csv");
    file << "Lambda,Passive_Poly(Pa),Passive_Exp(Pa),Active_Plasticity(Pa),Active_Stimulation(Pa)\n";

    // 検討用モデルのインスタンス化
    // 1. Zhao vs Paetsch (Passive比較用)
    IntegratedMuscleModel model_poly(IntegratedMuscleModel::POLY_PAETSCH, IntegratedMuscleModel::BIOPHYSICAL_STIMULATION);
    IntegratedMuscleModel model_exp(IntegratedMuscleModel::EXP_ZHAO, IntegratedMuscleModel::BIOPHYSICAL_STIMULATION);

    // 2. Plasticity vs Stimulation (Activeメカニズム比較用 - Paetschベース)
    IntegratedMuscleModel model_plastic(IntegratedMuscleModel::POLY_PAETSCH, IntegratedMuscleModel::PLASTICITY_ONLY);
    IntegratedMuscleModel model_stim(IntegratedMuscleModel::POLY_PAETSCH, IntegratedMuscleModel::BIOPHYSICAL_STIMULATION);

    double lambda_start = 1.0;
    double lambda_end = 1.2;
    int steps = 100;

    for (int i = 0; i <= steps; ++i) {
        double lam = lambda_start + (double)i / steps * (lambda_end - lambda_start);

        // A. Passive特性の比較 (alpha = 0)
        double s_poly = model_poly.calculate_stress(lam, 0.0);
        double s_exp  = model_exp.calculate_stress(lam, 0.0);

        // B. Activeメカニズムの比較 (alpha = 1.0)
        // ※ここではPassiveモデルはPolyで固定してActiveの違いを見る
        double s_plas = model_plastic.calculate_stress(lam, 1.0);
        double s_stim = model_stim.calculate_stress(lam, 1.0);

        file << lam << "," << s_poly << "," << s_exp << "," << s_plas << "," << s_stim << "\n";
    }

    std::cout << "Data generated: muscle_study.csv" << std::endl;
    return 0;
}