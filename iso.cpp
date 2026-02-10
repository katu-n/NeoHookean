#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>

class ConstitutiveModel {
private:
    double mu; // Shear modulus (Pa)
    double mu_p; // Dimensionless
    double k1; // Stress-like parameter (Pa)
    double k2; // Dimensionless exponential coefficient
    int model_type; // 0: Paetsch(Passive), 1: Zhao(Exponential)

public:
    ConstitutiveModel(int type) : model_type(type) {
        if (type == 0) {
            mu = 755.5; 
            mu_p = 86.13;
        } else if (type == 1) {
            mu = 154.0 * 1000.0; 
            k1 = 243.0 * 1000.0;
            k2 = 48.1; 
        }
    }

    double calculate_stress(double lambda) {
        double I4 = lambda * lambda; // 繊維方向の不変量 (線維角0度と仮定)
        double sigma_iso = mu * (pow(lambda, 2.0) - 1.0 / lambda);
        double sigma_fib = 0.0;

        if (I4 > 1.0) {
            if (model_type == 0) {
                sigma_fib = 2.0 * mu * mu_p * I4 * (I4 - 1.0);
            } 
            else if (model_type == 1) {
                double term = k2 * pow(I4 - 1.0, 2.0);
                sigma_fib = 2.0 * I4 * k1 * (I4 - 1.0) * std::exp(term);
            }
        }

        return sigma_iso + sigma_fib;
    }
};

int main() {
    std::ofstream file("comparison.csv");
    file << "Lambda,Stress_Paetsch(Pa),Stress_Zhao(Pa)\n";

    ConstitutiveModel modelPaetsch(0);
    ConstitutiveModel modelZhao(1);

    double lambda_start = 1.0;
    double lambda_end = 1.2; // 20% 伸展
    int steps = 100;

    for (int i = 0; i <= steps; ++i) {
        double lam = lambda_start + (double)i / steps * (lambda_end - lambda_start);
        
        double s_p = modelPaetsch.calculate_stress(lam);
        double s_z = modelZhao.calculate_stress(lam);

        file << lam << "," << s_p << "," << s_z << "\n";
    }
    
    std::cout << "Done. Check comparison.csv" << std::endl;
    return 0;
}