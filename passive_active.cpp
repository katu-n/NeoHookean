#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>

const double E_NUM = 2.718281;

class PassiveActive{
    private:
        double mu; //Matrix shear modulus (Pa)
        double mu_p; //DL
        double mu_a; //DL
        double lambda_0;
        double c1;
        double c2;

        double r; //Damage parameter
        double m_h; //Hystersis slope parameter (Pa)
        double a_r; //Reloading parameter (Pa)
        double b_u; //Secondary unloading parameter (Pa)

        double alpha; //Active level
        double lambda_max; //Max stretch reached on primary loading path
        double W0_max; //Max energy reached on primary loading path

        bool reloading;
        bool secondary_unloading;
        double lambda_reversal; // Stretch at load reversal point
        double Z_reversal; // Z value at reversal
        double W0_reversal; //Energy at reversal
    
    public:
        PassiveActive(bool active_mode){
            mu       = 755.5;
            mu_p     = 86.13;
            lambda_0 = 0.95;
            c1       = 470;
            c2       = 0.95;

            if(active_mode){
                alpha = 1.0;
                r     = 1.55;
                m_h   = 30000;
                a_r   = 30000;
                b_u   = 15000;
            } else {
                alpha = 0.0;
                r     = 1.05;
                m_h   = 2500;
                a_r   = 2500;
                b_u   = 2500;
            }

            lambda_max = 1.0;
            W0_max     = W0(1.0);
            reloading  = false;
            secondary_unloading = false;
        }

        void set_activetion(double a){
            alpha = a;
        }

        double ma(double lambda) const{
            double I4_hat = pow(lambda / lambda_0,2.0);
            double I4_0   = pow(lambda_0,2.0);
            double exponent = -(I4_hat - I4_0) /c2;
            return c1*std::exp(exponent);
        }

        double W0(double lambda) const{
            double term_iso = mu / 2.0 * (pow(lambda,2.0+2.0/lambda -3.0));

            double term_pass = 0.0;
            if(lambda > 1.0){
                term_pass = (mu*mu_p /2.0)*(1.0 - alpha)*pow(pow(lambda,2.0)-1.0,2.0);
            }

            double term_act = 0.0;
            double I4_hat   = pow(lambda / lambda_0, 2.0);
            double I4_0     = pow(lambda_0,2.0);
            double ma_val   = ma(lambda);

            term_act = (mu*ma_val /2.0)*alpha*pow(I4_hat-I4_0,2.0);

            return term_iso + term_pass + term_act;
        }

        double calc_sigma0(double lambda) const {
        double sigma_iso = mu * (pow(lambda, 2.0) - 1.0/lambda);

        double sigma_pass = 0.0;
        if (lambda > 1.0) {
            sigma_pass = 2.0 * mu * mu_p * (1.0 - alpha) * pow(lambda, 2.0) * (pow(lambda, 2.0) - 1.0);
        }

        double sigma_act = 0.0;
        double I4_hat = pow(lambda / lambda_0, 2.0);
        double I4_0 = pow(lambda_0, 2.0);
        double ma_val = ma(lambda);
        
        double bracket = 1.0 - (I4_hat - I4_0) / c2;
        
        sigma_act = 2.0 * mu * ma_val * alpha * I4_hat * (I4_hat - I4_0) * bracket;

        return sigma_iso + sigma_pass + sigma_act;
    }

    double update_stress(double lambda) {
            double W0_curr = W0(lambda);
            double sigma0 = calc_sigma0(lambda);
            double Z = 1.0;

            if (lambda >= lambda_max) {
                reloading = false;
                secondary_unloading = false;
                lambda_max = lambda;
                W0_max = W0_curr;
                Z = 1.0;
            } 
            else {
                double arg = (W0_max - W0_curr) / m_h;
                Z = 1.0 - (1.0 / r) * std::tanh(arg);
            }

            return Z * sigma0;
        }
};

int main() {
    std::ofstream file("muscle_response.csv");
    file << "Lambda,Sigma_Passive(Pa),Sigma_Active(Pa)\n";

    // Create two models to compare
    PassiveActive passiveModel(false); // Passive (alpha=0)
    PassiveActive activeModel(true);   // Active (alpha=1)

    // Simulation Protocol: Stretch from 1.0 to 1.15 then unload to 1.0
    double lambda_start = 1.0;
    double lambda_end = 1.15;
    int steps = 100;
    
    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        double lam = lambda_start + t * (lambda_end - lambda_start);
        
        double s_pass = passiveModel.update_stress(lam);
        double s_act = activeModel.update_stress(lam);
        
        file << lam << "," << s_pass << "," << s_act << "\n";
    }

    // 2. Unloading Phase
    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        double lam = lambda_end - t * (lambda_end - lambda_start);
        
        double s_pass = passiveModel.update_stress(lam);
        double s_act = activeModel.update_stress(lam);
        
        file << lam << "," << s_pass << "," << s_act << "\n";
    }

    std::cout << "Simulation finished. Data written to muscle_response.csv" << std::endl;
    return 0;
}