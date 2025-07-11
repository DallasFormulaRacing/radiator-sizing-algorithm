#include "../inc/radiator.hpp"
#include <cmath>

Radiator::Radiator(double numTubes, double tubeHeight, double radLength,
                                     double tubeWidth, double finDist,
                                     double finHeight, double finWidth)
    : NumberOfTubes(numTubes), TubeHeight(tubeHeight), RadiatorLength(radLength),
      TubeWidth(tubeWidth), FinDistance(finDist), FinHeight(finHeight), FinWidth(finWidth) {}

double Radiator::run() {
    const double rho_w = 1015.57, mu_w = 0.000344, Cp_w = 4181.92, k_w = 0.615;
    const double rho_a = 1.137, mu_a = 0.0000191, Cp_a = 1010.16;
    const double k_a = 0.029, V_air = 4.47, V_dot_cool = 0.000827;
    const double T_air = 310.9, T_cool = 394.2;

    double finRows = NumberOfTubes - 1;
    double finCols = RadiatorLength / FinDistance;
    double A_air = finRows * finCols * (2 * FinDistance * FinHeight + 2 * FinHeight * FinWidth);
    double A_cool = NumberOfTubes * (2 * TubeHeight * RadiatorLength + 2 * TubeWidth * RadiatorLength);
    double A_total = A_air + A_cool;

    double Dh_cool = 4 * TubeWidth * TubeHeight / (2 * (TubeWidth + TubeHeight));
    double V_cool = V_dot_cool / (NumberOfTubes * TubeWidth * TubeHeight);
    double Re_cool = Reynolds(rho_w, V_cool, Dh_cool, mu_w);
    double Pr_cool = Prandtl(Cp_w, mu_w, k_w);
    double Nu_cool = Nusselt(Re_cool, Pr_cool);
    double h_cool = HTC(Nu_cool, k_w, Dh_cool);

    double Dh_air = 4 * FinHeight * FinDistance / (2 * (FinDistance + 8 * FinDistance));
    double Re_air = Reynolds(rho_a, V_air, Dh_air, mu_a);

    double m_dot_cool = V_dot_cool * rho_w;
    double m_dot_air = 1.1086 * rho_a;
    double C_cool = Cp_w * m_dot_cool;
    double C_air = Cp_a * m_dot_air;
    double Cmin = std::min(C_cool, C_air);
    double Cmax = std::max(C_cool, C_air);
    double Cr = Cmin / Cmax;

    double nfha = 266.94;
    double U = 1 / (1 / (h_cool * A_cool) + 1 / (nfha * A_air));
    double Ntu = NTU(U, Cmin);
    double eff = Effectiveness(Cmax, Cmin, Cr, Ntu);

    double ITD = T_cool - T_air;
    return eff * Cmin * ITD;
}

double Radiator::getRadiatorHeight() const {
    return TubeHeight * NumberOfTubes + FinHeight + (NumberOfTubes + 1);
}

double Radiator::Reynolds(double rho, double v, double D, double mu) const {
    return (rho * v * D) / mu;
}

double Radiator::Prandtl(double Cp, double mu, double k) const {
    return (Cp * mu) / k;
}

double Radiator::Nusselt(double Re, double Pr) const {
    return 0.023 * std::pow(Re, 0.8) * std::pow(Pr, 1.0/3.0);
}

double Radiator::HTC(double Nu, double k, double D) const {
    return (Nu * k) / D;
}

double Radiator::NTU(double U, double Cmin) const {
    return U / Cmin;
}

double Radiator::Effectiveness(double Cmax, double Cmin, double Cr, double Ntu) const {
    return 1.0 - std::exp(-1.0 * (Cmax / Cmin) * (1.0 - std::exp(-Cr * Ntu)));
}
