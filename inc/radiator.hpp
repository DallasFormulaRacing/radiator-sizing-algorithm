#pragma once

class Radiator {

    public:
    Radiator(double numTubes, double tubeHeight, double radLength,
                      double tubeWidth, double finDist,
                      double finHeight, double finWidth);

    double run(); 
    double getRadiatorHeight() const;

private:
    
    double NumberOfTubes, TubeHeight, RadiatorLength, TubeWidth;
    double FinDistance, FinHeight, FinWidth;

    double Reynolds(double density, double velocity, double diameter, double viscosity) const;
    double Prandtl(double Cp, double mu, double k) const;
    double Nusselt(double Re, double Pr) const;
    double HTC(double Nu, double k, double D) const;
    double NTU(double U, double Cmin) const;
    double Effectiveness(double Cmax, double Cmin, double Cr, double Ntu) const;

};