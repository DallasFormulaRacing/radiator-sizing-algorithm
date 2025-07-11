#include <cmath>
// #include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

//////////////////////////////////
// FUNCTIONS DEFINING EQUATIONS //
//////////////////////////////////

double HeatTransEq(double e, double Cmin, double ITD) {
  double q = e * Cmin * ITD;
  return q;
}

double UHTE(double HeatTransferCoefficientCoolant, double AirSurfaceArea,
            double nfha, double CoolantSurfaceArea) {
  double rightside =
      (1 / (HeatTransferCoefficientCoolant * CoolantSurfaceArea)) +
      (1 / (nfha * AirSurfaceArea));
  double UHTE = 1 / rightside;
  return UHTE;
}

double ReynoldsEquation(double density, double velocity,
                        double HydraulicDiameter, double viscosity) {
  double ReynoldsNum = (density * velocity * HydraulicDiameter) / viscosity;
  return ReynoldsNum;
}

double HydraulicDiameter(double PipeArea, double Perimeter) {
  double HydraulicD = (4 * PipeArea) / (Perimeter);
  return HydraulicD;
}

double DittusBoelterEquation(double ReynoldsNum, double PrandtlNum) {
  double NusseltNum =
      .023 * (pow(ReynoldsNum, (0.8)) * pow(PrandtlNum, (0.333333)));
  return NusseltNum;
}

double PrandtlEquation(double SpecificHeat, double viscosity,
                       double ThermalConductivity) {
  double PrandtlNum = (SpecificHeat * viscosity) / ThermalConductivity;
  return PrandtlNum;
}

double NusseltEquation(double HeatTransferCoefficient, double HydraulicDiameter,
                       double ThermalConductivity) {
  double NusseltNum =
      (HeatTransferCoefficient * HydraulicDiameter) / ThermalConductivity;
  return NusseltNum;
}

double NtuEquation(double HeatTransferCoefficient, double Cmin) {
  double Ntu = (HeatTransferCoefficient / Cmin);
  return Ntu;
}

double ENtuEquation(double Cmax, double Cmin, double Cratio, double Ntu) {
  double e = 2.7183;
  double secondexponent = -1 * Cratio * Ntu;
  double exponent = -1 * (Cmax / Cmin) * (1 - pow(e, secondexponent));
  double ENtuNum = 1 - pow(e, exponent);
  return ENtuNum;
}

double ITDEquation(double CoolantTemperature, double AirTemperature) {
  double ITD = CoolantTemperature - AirTemperature;
  return ITD;
}

double HeatTransferCoefficient(double NusseltNum, double ThermalConductivity,
                               double HydraulicDiameter) {
  double HeatTransferCoefficient =
      (NusseltNum * ThermalConductivity) / HydraulicDiameter;
  return HeatTransferCoefficient;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Calculation(double NumberOfTubes, double TubeHeight,
                   double RadiatorLength, double TubeWidth, double FinDistance,
                   double FinHeight, double FinWidth) {
  const double DensityOfWater =
      1015.57; // KG/M^3 at 85 C //They have 1015.57 968.6
  const double ViscosityOfWater =
      .000344082; // Pa s //They have.000744082 .0003335
  const double SpecificHeatOfWater = 4181.92;        // THIS IS FOR GLYCOL
  const double ThermalConductivityOfWater = .615098; // W/mK at 85 C .6589
  const double DensityOfAir =
      1.13731; // KG/M^3 at 85 C NEED TO CHANGE TO 85DEG VALUE
  const double ViscosityOfAir = .00001912;  // Pa s at 85 C
  const double SpecificHeatOfAir = 1010.16; // J/g/K NEED TO REPLACE WITH 85DEG

  const double AirTemperature = 310.928;
  const double AirVolumetricFlow = 1.10860;
  const double CoolantTemperature = 394.261; // ADD VALUES HEREEEEEE
  const double VelocityOfAir = 4.4704;
  const double CoolantVolumetricFlow = .000827;

  double NumRowsOfFins = NumberOfTubes - 1;
  // printf("Numrowsfins: %g\n", NumRowsOfFins);

  double TotalNumberOfAirPassages =
      NumRowsOfFins * (RadiatorLength / FinDistance);
  //  printf("Total number Airpass: %g\n", TotalNumberOfAirPassages);

  double AirSurfaceArea =
      TotalNumberOfAirPassages *
      (2 * FinDistance * FinHeight + 2 * FinHeight * FinWidth);
  // printf("AirSurfaceArea: %g\n", AirSurfaceArea);
  double CoolantSurfaceArea =
      NumberOfTubes *
      (2 * TubeHeight * RadiatorLength +
       2 * TubeWidth *
           RadiatorLength); // These three equations are just to find surface
                            // area printf("CoolantSurfaceArea: %g\n",
                            // CoolantSurfaceArea);
  double AreaTotal = AirSurfaceArea + CoolantSurfaceArea;
  // printf("AreaTotal: %g\n", AreaTotal);

  double PipeAreaCoolant = TubeWidth * TubeHeight;
  // printf("PipeAreaCoolant: %g\n", PipeAreaCoolant);
  double Perimeter =
      2 *
      (TubeWidth +
       TubeHeight); // These three equations are just to find hydraulic diameter
  // printf("perimeter: %g\n", Perimeter);
  double HydraulicD = HydraulicDiameter(PipeAreaCoolant, Perimeter);
  // printf("hydraulicD: %g\n", HydraulicD);

  double VelocityOfCoolant =
      CoolantVolumetricFlow / (NumberOfTubes * PipeAreaCoolant);
  // printf("velocity of coolant: %g\n", VelocityOfCoolant);
  double ReynoldsNumCoolant = ReynoldsEquation(
      DensityOfWater, VelocityOfCoolant, HydraulicD, ViscosityOfWater);
  // printf("reynoldsnumcoolant: %g\n", ReynoldsNumCoolant);
  double PrandtlNumCoolant = PrandtlEquation(
      SpecificHeatOfWater, ViscosityOfWater,
      ThermalConductivityOfWater); // These use the Reynolds equation and
                                   // Prandtl Equation to find the heat transfer
                                   // coefficient
  // printf("prandtlnumcool: %g\n", PrandtlNumCoolant);
  double NusseltNumCoolant =
      DittusBoelterEquation(ReynoldsNumCoolant, PrandtlNumCoolant);
  // printf("nusseltnumcoolant: %g\n", NusseltNumCoolant);
  double HeatTransferCoefficientCoolant = HeatTransferCoefficient(
      NusseltNumCoolant, ThermalConductivityOfWater, HydraulicD);
  // printf("HTCcoolant: %g\n", HeatTransferCoefficientCoolant);

  double PipeAreaAir = FinHeight * FinDistance;
  // printf("Pipeareaair: %g\n", PipeAreaAir);
  double PerimeterAir = 2 * (FinDistance + 8 * FinDistance);
  // printf("perimeterair: %g\n", PerimeterAir);
  double HydraulicDiameterAir = HydraulicDiameter(
      PipeAreaAir,
      PerimeterAir); // FOUND REYNOLDS NUM FOR AIR
                     // printf("hydrodiamair: %g\n", HydraulicDiameterAir);
  double ReynoldsNumAir = ReynoldsEquation(
      DensityOfAir, VelocityOfAir, HydraulicDiameterAir, ViscosityOfAir);
  // printf("reynoldsnumair: %g\n", ReynoldsNumAir);

  double MassFlowRateCoolant = CoolantVolumetricFlow * DensityOfWater;
  // printf("massflowratecool: %g\n", MassFlowRateCoolant);
  double MassFlowRateAir = AirVolumetricFlow * DensityOfAir;
  // printf("massflowrateair: %g\n", MassFlowRateAir);
  double ThermalCapacityCoolant = SpecificHeatOfWater * MassFlowRateCoolant;
  // printf("thermalcapacitycool: %g\n", ThermalCapacityCoolant);
  double ThermalCapacityAir = SpecificHeatOfAir * MassFlowRateAir;
  // printf("thermalcapacityair: %g\n", ThermalCapacityAir);
  double CapacityRatio = ThermalCapacityAir / ThermalCapacityCoolant;
  // printf("capcityratio: %g\n", CapacityRatio);

  double ITDValue = ITDEquation(CoolantTemperature, AirTemperature);
  // printf("ITD: %g\n", ITDValue);
  double nfha = 266.9425414086; // Converted Number given by paper to SI units
  // printf("nfha: %g\n", nfha);
  double OverallHTC = UHTE(HeatTransferCoefficientCoolant, AirSurfaceArea, nfha,
                           CoolantSurfaceArea);
  // printf("overallhtc: %g\n", OverallHTC);
  double NTU = NtuEquation(
      OverallHTC,
      ThermalCapacityAir); // I believe Cmin to be ThermalCapacityAir
  // printf("ntu: %g\n", NTU);
  double E = ENtuEquation(ThermalCapacityCoolant, ThermalCapacityAir,
                          CapacityRatio, NTU);
  // printf("E: %g\n", E);

  double q = HeatTransEq(E, ThermalCapacityAir, ITDValue);

  // printf("The value for q is: %g\n", q);
  // printf("The value for q should be: 60674.152");

  return q;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void sendtofile(double NumberOfTubes[], double TubeHeight[],
                double RadiatorLength[], double TubeWidth[],
                double FinDistance[], double FinHeight[], double FinWidth[],
                double RadiatorHeight[], double qvalue[], int ArraySize) {
  FILE *resultsfile;
  resultsfile = fopen("onevar.txt", "w");
  if (resultsfile == NULL) {
    printf("Error opening file...");
    return;
  }

  int i = 0;
  fprintf(resultsfile, " Number Of Tubes \n\n");
  while (i < ArraySize) {
    fprintf(resultsfile, "%g,", NumberOfTubes[i]);
    i++;
  }
  i = 0;
  fprintf(resultsfile, "\n\n Tube Height \n\n");
  while (i < ArraySize) {
    fprintf(resultsfile, "%g, ", TubeHeight[i]);
    i++;
  }
  i = 0;
  fprintf(resultsfile, "\n\n Radiator Length \n\n");
  while (i < ArraySize) {
    printf("RADLENGTH %g\n\n", RadiatorLength[i]); // for debugging
    printf("INDEX %d\n\n", i);                     // for debugging
    printf("ARRAYSIZE %d\n\n", ArraySize);         // for debugging

    fprintf(resultsfile, "%g, ", RadiatorLength[i]);
    i++;
  }
  i = 0;
  fprintf(resultsfile, "\n\n Tube Width \n\n");
  while (i < ArraySize) {
    fprintf(resultsfile, "%g, ", TubeWidth[i]);
    i++;
  }
  i = 0;
  fprintf(resultsfile, "\n\n Fin Distance \n\n");
  while (i < ArraySize) {
    fprintf(resultsfile, "%g, ", FinDistance[i]);
    i++;
  }
  i = 0;
  fprintf(resultsfile, "\n\n Fin Height \n\n");
  while (i < ArraySize) {
    fprintf(resultsfile, "%g, ", FinHeight[i]);
    i++;
  }
  i = 0;
  fprintf(resultsfile, "\n\n Fin Width \n\n");
  while (i < ArraySize) {
    fprintf(resultsfile, "%g, ", FinWidth[i]);
    i++;
  }
  i = 0;
  fprintf(resultsfile, "\n\n Radiator Height \n\n");
  while (i < ArraySize) {
    fprintf(resultsfile, "%g, ", RadiatorHeight[i]);
    i++;
  }
  i = 0;
  fprintf(resultsfile, "\n\n Q Value \n\n");
  while (i < ArraySize) {
    fprintf(resultsfile, "%g, ", qvalue[i]);
    i++;
  }
  fclose(resultsfile);
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
  int storagecounter = 0;
  int run = 1;
  int NumberOfSuccess = 0;

  double iterations = 100;

  double RadiatorLengthmin = 0;
  double RadiatorLengthmax = 1;
  double RadiatorLengthInc =
      (double)(RadiatorLengthmax - RadiatorLengthmin) / iterations;
  double RadiatorLength = RadiatorLengthmin;

  double TubeHeightmin = 0.00156261;

  double TubeHeight = TubeHeightmin;

  double TubeWidthmin = 0.0246063;

  double TubeWidth = TubeWidthmin;

  double FinDistancemin = 0.0015875;

  double FinDistance = FinDistancemin;

  double FinHeightmin = 0.118813;

  double FinHeight = FinHeightmin;

  double FinWidthmin = 0.0246063;

  double FinWidth = FinWidthmin;

  double NumberOfTubesmin = 60;

  double NumberOfTubes = NumberOfTubesmin;

  double qvalue;
  double RadiatorHeight;

  while (run == 1) {

    printf("1  %g  %g  %g  %g  %g  %g  %g\n\n", NumberOfTubes, TubeHeight,
           RadiatorLength, TubeWidth, FinDistance, FinHeight,
           FinWidth); // for debugging

    qvalue = Calculation(NumberOfTubes, TubeHeight, RadiatorLength, TubeWidth,
                         FinDistance, FinHeight, FinWidth);
    NumberOfSuccess++;

    RadiatorLength = RadiatorLength + RadiatorLengthInc;
    if (RadiatorLength > RadiatorLengthmax) {
      run++;
    }
  }

  RadiatorLength = RadiatorLengthmin;

  double NumberOfTubesArray[NumberOfSuccess];
  double TubeHeightArray[NumberOfSuccess];
  double RadiatorLengthArray[NumberOfSuccess];
  double TubeWidthArray[NumberOfSuccess];
  double FinDistanceArray[NumberOfSuccess];
  double FinHeightArray[NumberOfSuccess];
  double FinWidthArray[NumberOfSuccess];
  double RadiatorHeightArray[NumberOfSuccess];
  double QValueArray[NumberOfSuccess];

  while (run == 2) {

    RadiatorHeight =
        TubeHeight * NumberOfTubes + FinHeight + (NumberOfTubes + 1);

    qvalue = Calculation(NumberOfTubes, TubeHeight, RadiatorLength, TubeWidth,
                         FinDistance, FinHeight, FinWidth);
    printf("2  %g  %g  %g  %g  %g  %g  %g  %g\n\n", NumberOfTubes, TubeHeight,
           RadiatorLength, TubeWidth, FinDistance, FinHeight, FinWidth,
           qvalue); // for debugging
    RadiatorLength = RadiatorLength + RadiatorLengthInc;

    NumberOfTubesArray[storagecounter] = NumberOfTubes;
    TubeHeightArray[storagecounter] = TubeHeight;
    RadiatorLengthArray[storagecounter] = RadiatorLength;
    TubeWidthArray[storagecounter] = TubeWidth;
    FinDistanceArray[storagecounter] = FinDistance;
    FinHeightArray[storagecounter] = FinHeight;
    FinWidthArray[storagecounter] = FinWidth;
    RadiatorHeightArray[storagecounter] = RadiatorHeight;
    QValueArray[storagecounter] = qvalue;
    storagecounter++;
    if (RadiatorLength > RadiatorLengthmax) {
      run++;
    }
  }

  sendtofile(NumberOfTubesArray, TubeHeightArray, RadiatorLengthArray,
             TubeWidthArray, FinDistanceArray, FinHeightArray, FinWidthArray,
             RadiatorHeightArray, QValueArray, NumberOfSuccess);

  return 0;
}
