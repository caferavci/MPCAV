#pragma once
#include "gurobi_c++.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm> 
#include "DynamicArray.h"
#include <string>

class MPC
{
private:


	double	simulationTimeStep = 0.1, simulationStartTime = 0, CACC_k = 1.0, CACC_kA = 1.0, CACC_kB = 0.58, CACC_kC = 0.1, CACC_H = 0.5, minimunDistance = 2.0, minimumDisturbanceDistance = 1;
	double	ACCLimit = 2.0, DECCLimit = -3.0, x0Acc = 0.0, xMin = 0.0, xMax = 1000.0, uMin = 0.0, uMax = 33.0, epsilonMin = 0.0, epsilonMax = 10.0, exponentialAlpha = 0.95, exponentialBeta = 0.99, avgHOVSpeed = 0.0;
	double	alfaobj = 100, gama = 100;
	double finalVelocity = 0.0;
	double displacement = 0.0;
	const int vehicleCounts, tHor, simHor;
	char comaType = ',';
	const int sizeRefG = 4;
	double refG[4] = { CACC_kC, CACC_kB, -CACC_kC,  -CACC_kB - (CACC_kC * CACC_H) };

	double* vehicles_X_V;
	double* vehiclesLength;

	double leaderAcc;
	double* optimizedCAVs;
	double* optimizedCAVsACC;
	double** trjHOV;

	void fillB(double* matrix, const int* vectorLength, double timeStep);
	void fillB_2(double** matrix, const int* vectorLength, double timeStep);
	void getSystemMatrix(double** matrix, double timeStep, const int* matrixSize);
	void fillRef5vehicles(double** matrix, double timeStep);
	void getSystemMatrix_2(double** matrix, double timeStep, const int* matrixSize);
	void fillRef5vehicles_2(double** matrix, double timeStep);
	double getR();
	void writeTrajectories(GRBVar* CAVs, int matrixLength, int vectorLength, std::string file_name);// CSV file
	void fillOutput2(GRBVar* CAVs, int matrixLength, int vectorLength);

public:
	MPC();
	MPC(double* vehicles_X_V, double* vehiclesLength, double leaderAcc, const int vehicleCounts, const int tHor, double simulationTimeStep);
	bool optimizePlatoonwithDisturbanceMinACC(double* HDV_X_V, const int HOVsize, const int initialTime, double* speedLimits);
	void setDisturbanceTrajectories(double** disturbanceTrj, double* HDV_X_V, const int HOVsize, const int initialTime, int disturbanceSize);
	void getOptimizedTrajectories(double* Matrix);

	//vehicleCounts
	int getvehicleCounts();
	//simulationTimeStep
	void setfinalVelocity(double s);
	double getfinalVelocity();
	//simulationTimeStep
	void setdisplacement(double s);
	double getdisplacement();
	//simulationTimeStep
	void setsimulationTimeStep(double s);
	double getsimulationTimeStep();
	//simulationStartTime
	void setsimulationStartTime(double s);
	double getsimulationStartTime();
	//minimunDistance
	void setminimunDistance(double s);
	double getminimunDistance();
	//minimunDisturbanceDistance
	void setminimunDisturbanceDistance(double s);
	double getminimunDisturbanceDistance();
	//CACC_k
	void setCACC_k(double s);
	double getCACC_k();
	//CACC_kA
	void setCACC_kA(double s);
	double getCACC_kA();
	//CACC_kB
	void setCACC_kB(double s);
	double getCACC_kB();
	//CACC_kC
	void setCACC_kC(double s);
	double getCACC_kC();
	//CACC_H
	void setCACC_H(double s);
	double getCACC_H();
	//ACCLimit
	void setACCLimit(double s);
	double getACCLimit();
	//DECCLimit
	void setDECCLimit(double s);
	double getDECCLimit();
	//x0Acc
	void setx0Acc(double s);
	double getx0Acc();
	//xMin
	void setxMin(double s);
	double getxMin();
	//xMax
	void setxMax(double s);
	double getxMax();
	//uMin
	void setuMin(double s);
	double getuMin();
	//uMax
	void setuMax(double s);
	double getuMax();
	//epsilonMin
	void setepsilonMin(double s);
	double getepsilonMin();
	//epsilonMax
	void setepsilonMax(double s);
	double getepsilonMax();
	//tHor
	int gettHor();
	//alfaobj
	void setalfaobj(double s);
	double getalfaobj();
	//gama
	void setgama(double s);
	double getgama();
	//leaderACC
	void setleaderACC(double s);
	double getleaderACC();
	void getDisturbanceTrajectories(double** disturbanceTrj);
	void getOptimizedACCs(double* Matrix);
};