#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>

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

// 状態管理用の列挙型
enum State {
    PRIMARY_LOADING,
    UNLOADING,
    RELOADING,
    SECONDARY_UNLOADING
};

class MuscleModel {
private:
    MaterialParams p_pass; 
    MaterialParams p_act;  
    MaterialParams p_curr; 

    double alpha;      // Activation level (0.0 to 1.0)
    
    // 状態管理
    State current_state;
    double lambda_prev; // 前回の伸展比（方向判定用）

    // 履歴変数
    double lambda_max; // 一次負荷経路の最大到達点
    double W0_max;     // 一次負荷経路の最大エネルギー

    // Reloading / Secondary Unloading 用の記憶変数
    double W0_reversal_reload;     // 再負荷開始時のエネルギー
    double eta_reversal_reload;    // 再負荷開始時の eta
    double W0_reversal_unload;     // 二次除荷開始時のエネルギー
    double eta_reversal_unload;    // 二次除荷開始時の eta
    double alpha_prev = 0.0; 

public:
    MuscleModel() {
        p_pass = {755.5, 86.13, 0.95, 470.0, 0.95, 1.05, 2500.0, 2500.0, 2500.0}; 
        p_act  = {755.5, 86.13, 0.95, 470.0, 0.95, 1.55, 30000.0, 30000.0, 15000.0};
        
        alpha = 0.0;
        update_params();
        
        // 初期状態
        current_state = PRIMARY_LOADING;
        lambda_prev = 1.0;
        lambda_max = 1.0;
        W0_max = 0.0; 
        
        // 内部変数の初期化
        W0_reversal_reload = 0.0;
        eta_reversal_reload = 1.0;
        W0_reversal_unload = 0.0;
        eta_reversal_unload = 1.0;
    }

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

        p_curr.r = (1.0 - alpha) * p_pass.r + alpha * p_act.r;
        p_curr.m = (1.0 - alpha) * p_pass.m + alpha * p_act.m;
        p_curr.a = (1.0 - alpha) * p_pass.a + alpha * p_act.a;
        p_curr.b = (1.0 - alpha) * p_pass.b + alpha * p_act.b;
    }

    double get_mu_a(double lambda) const {
        double I4_hat = std::pow(lambda / p_curr.lambda_0, 2.0);
        double I4_0   = std::pow(p_curr.lambda_0, 2.0);
        double exponent = -(I4_hat - I4_0) / p_curr.c2;
        return p_curr.c1 * std::exp(exponent);
    }

    double get_W0(double lambda) const {
        double I1 = std::pow(lambda, 2.0) + 2.0 / lambda;
        double W_iso = (p_curr.mu / 2.0) * (I1 - 3.0);

        double W_pass = 0.0;
        if (lambda > 1.0) {
            double I4 = std::pow(lambda, 2.0);
            W_pass = (p_curr.mu * p_curr.mu_p / 2.0) * (1.0 - alpha) * std::pow(I4 - 1.0, 2.0);
        }

        double W_act = 0.0;
        double I4_hat = std::pow(lambda / p_curr.lambda_0, 2.0);
        double I4_0   = std::pow(p_curr.lambda_0, 2.0);
        double mu_a_val = get_mu_a(lambda);
        W_act = (p_curr.mu * mu_a_val / 2.0) * alpha * std::pow(I4_hat - I4_0, 2.0);

        return W_iso + W_pass + W_act;
    }

    double get_sigma0(double lambda) const {
        double sigma_iso = p_curr.mu * (std::pow(lambda, 2.0) - 1.0 / lambda);

        double sigma_pass = 0.0;
        if (lambda > 1.0) {
            sigma_pass = 2.0 * p_curr.mu * p_curr.mu_p * (1.0 - alpha) * std::pow(lambda, 2.0) * (std::pow(lambda, 2.0) - 1.0);
        }

        double sigma_act = 0.0;
        double I4_hat = std::pow(lambda / p_curr.lambda_0, 2.0);
        double I4_0   = std::pow(p_curr.lambda_0, 2.0);
        double mu_a_val = get_mu_a(lambda);
        
        double bracket = 1.0 - (I4_hat - I4_0) / (2.0 * p_curr.c2); 
        sigma_act = 2.0 * p_curr.mu * mu_a_val * alpha * I4_hat * (I4_hat - I4_0) * bracket;

        return sigma_iso + sigma_pass + sigma_act;
    }

    // 現在のetaを計算するヘルパー関数（状態判定のために必要）
    double calculate_current_eta(double lambda, double W0_curr) {
        if (current_state == PRIMARY_LOADING) {
            return 1.0;
        } else if (current_state == UNLOADING) {
             double W0_max_curr = get_W0(lambda_max);
             double arg = (W0_max_curr - W0_curr) / p_curr.m;
             return 1.0 - (1.0 / p_curr.r) * std::tanh(arg);
        } else if (current_state == RELOADING) {
             double arg = (W0_curr - W0_reversal_reload) / p_curr.a;
             // eta = eta_r * (1 + ((1-eta_1)/eta_1) * tanh(...))
             double term = ((1.0 - eta_reversal_reload) / eta_reversal_reload) * std::tanh(arg);
             return eta_reversal_reload * (1.0 + term);
        } else if (current_state == SECONDARY_UNLOADING) {
             double arg = (W0_reversal_unload - W0_curr) / p_curr.b;
             return eta_reversal_unload * (1.0 - std::tanh(arg));
        }
        return 1.0;
    }
    double compute_stress(double lambda) {
        double d_lambda = lambda - lambda_prev;
        double W0_curr = get_W0(lambda);

        // --- 状態遷移ロジックの改善 ---
        // 伸展比の増減による基本的な状態遷移
        if (std::abs(d_lambda) > 1e-9) { 
            if (d_lambda < 0) { // 除荷方向
                 if (current_state == PRIMARY_LOADING) {
                    current_state = UNLOADING;
                 } else if (current_state == RELOADING) {
                    current_state = SECONDARY_UNLOADING;
                    W0_reversal_unload = get_W0(lambda_prev); // 反転直前の値を保存
                    eta_reversal_unload = calculate_current_eta(lambda_prev, W0_reversal_unload);
                 }
            } else { // 負荷方向
                 if (current_state == UNLOADING || current_state == SECONDARY_UNLOADING) {
                    current_state = RELOADING;
                    W0_reversal_reload = get_W0(lambda_prev);
                    eta_reversal_reload = calculate_current_eta(lambda_prev, W0_reversal_reload);
                 }
            }
        }
        lambda_prev = lambda;

        // --- 応力計算 ---
        double sigma0 = get_sigma0(lambda);
        double eta = 1.0;

        // Active化によりエネルギーが過去の最大を超えた場合の「自然なLoading復帰」判定
        // これは全ての状態でチェックすべき（特にUNLOADING中）
        if (W0_curr > W0_max) {
            current_state = PRIMARY_LOADING;
            W0_max = W0_curr;
            lambda_max = lambda;
        }

        if (current_state == PRIMARY_LOADING) {
            // Loading中は常に最大点を更新
            if (lambda > lambda_max) { 
                 lambda_max = lambda;
                 W0_max = W0_curr;
            }
            eta = 1.0; 
        } 
        else if (current_state == UNLOADING) {
            // ★重要修正: W0_max を再計算せず、履歴値(W0_max)をそのまま使う
            // Active化で W0_curr が増えると、(W0_max - W0_curr) が小さくなり、eta は 1.0 に近づく。
            // さらに W0_curr > W0_max になれば上のif文で PRIMARY_LOADING になる。
            // これにより滑らかな遷移が実現する。
            
            double arg = (W0_max - W0_curr) / p_curr.m;
            eta = 1.0 - (1.0 / p_curr.r) * std::tanh(arg);
        } 
        else if (current_state == RELOADING) {
            double arg = (W0_curr - W0_reversal_reload) / p_curr.a;
            double term = ((1.0 - eta_reversal_reload) / eta_reversal_reload) * std::tanh(arg);
            eta = eta_reversal_reload * (1.0 + term);
            
            // Primary Pathへの合流判定
            if (eta >= 1.0 || lambda >= lambda_max || W0_curr >= W0_max) {
                current_state = PRIMARY_LOADING;
                eta = 1.0;
                W0_max = W0_curr; // エネルギー準位を更新
                lambda_max = lambda;
            }
        } 
        else if (current_state == SECONDARY_UNLOADING) {
            double arg = (W0_reversal_unload - W0_curr) / p_curr.b;
            eta = eta_reversal_unload * (1.0 - std::tanh(arg));
        }

        return eta * sigma0;
    }
};

double neo_hookean_stress(double lambda, double mu) {
    return mu * (std::pow(lambda, 2.0) - 1.0 / lambda);
}

int main() {
    std::ofstream file("muscle_comparison_v2.csv");
    file << "Time,Lambda,Alpha,Stress_Model(MPa),Stress_NH(MPa),Phase\n";

    MuscleModel muscle;
    double mu_base = 755.5; 

    double lambda_start = 1.0;
    double lambda_peak = 1.2;
    double lambda_end = 1.0;

    int steps = 100;
    double time = 0.0;
    double dt = 0.1;

    // --- Phase 1: Passive Loading (1.0 -> 1.2) ---
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
    double activation_start_lam = 1.08;
    double activation_end_lam   = 1.00;

    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        double lam = lambda_peak - t * (lambda_peak - lambda_end);
        
        double alpha = 0.0;
        if (lam < activation_start_lam && lam > activation_end_lam) {
            alpha = (activation_start_lam - lam) / (activation_start_lam - activation_end_lam);
        } else if (lam <= activation_end_lam) {
            alpha = 1.0;
        }

        muscle.set_activation(alpha);
        
        double s_model = muscle.compute_stress(lam);
        double s_nh = neo_hookean_stress(lam, mu_base);

        file << time << "," << lam << "," << alpha << "," 
             << s_model / 1e6 << "," << s_nh / 1e6 << ",2_UnloadTrans\n";
        time += dt;
    }

    // --- Phase 3: Active Reloading (1.0 -> 1.2) ---
    // ここでActive状態での Reloading -> Primary Loading への復帰が発生する可能性がある
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
    // ここで Phase 3 の到達点によっては Secondary Unloading (param b) か Unloading (param m) かに分岐する
    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        double lam = lambda_peak - t * (lambda_peak - lambda_end);
        
        double s_model = muscle.compute_stress(lam);
        double s_nh = neo_hookean_stress(lam, mu_base);

        file << time << "," << lam << "," << 1.0 << "," 
             << s_model / 1e6 << "," << s_nh / 1e6 << ",4_ActiveUnload\n";
        time += dt;
    }

    std::cout << "Simulation finished. Results in muscle_comparison_v2.csv" << std::endl;
    return 0;
}