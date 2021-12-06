#include "MPC.h"
#define MAX_SIZE 16382
using namespace std;
MPC::MPC() :vehicleCounts(1), tHor(50), simHor(200)
{


}
MPC::MPC(double* vehicles_X_V, double* vehiclesLength, double leaderAcc, const int vehicleCounts, const int tHor, double simulationTimeStep) : vehicleCounts(vehicleCounts), tHor(tHor), simHor(simHor)
{
    this->vehiclesLength = vehiclesLength;
    this->vehicles_X_V = vehicles_X_V;
    this->leaderAcc = leaderAcc;
    this->simulationTimeStep = simulationTimeStep;
}

void MPC::getDisturbanceTrajectories(double** disturbanceTrj)
{
    {
        ofstream FILE;
        //writing trajectories
        std::string fileName_ = "Predicted_Trajectoriess.csv";
        FILE.open(fileName_);
        int t_ = 0;
        int loopSize = 0;
        int countTimer = 0;
        for (int i = 0; i < tHor; i++)
        {
            FILE << trjHOV[i][0] << ',' << trjHOV[i][1] << ',' << trjHOV[i][2] << endl;
        }
        FILE.close();
    }

}

bool MPC::optimizePlatoonwithDisturbanceMinACC(double* HDV_X_V, const int HOVsize, const int initialTime, double* speedLimits)
{

    GRBEnv* CAVenv = new GRBEnv();
    GRBModel* CAVmodel = new GRBModel(*CAVenv);
    //CAVmodel->set(GRB_IntParam_Method, 0);
    int disturbanceSize = tHor;
    trjHOV = AllocateDynamicArray<double>(tHor, 3);
    setDisturbanceTrajectories(trjHOV, HDV_X_V, HOVsize, initialTime, disturbanceSize);

    const int gVectorLength = vehicleCounts * 2;
    const int gVectorCount = vehicleCounts - 1;
    const int vectorLength = vehicleCounts * 2;
    const int QhorLength = vectorLength * tHor;

    double** matrixA = AllocateDynamicArray<double>(vectorLength, vectorLength);
    double* matrixB = AllocateDynamicVector<double>(vectorLength);



    GRBVar* xN = AllocateDynamicVector<GRBVar>(QhorLength);
    GRBVar* epsilonN = AllocateDynamicVector<GRBVar>(QhorLength);
    GRBVar* MuN = AllocateDynamicVector<GRBVar>(QhorLength);
    GRBVar* aN = AllocateDynamicVector<GRBVar>(QhorLength);
    GRBVar* thetaN = AllocateDynamicVector<GRBVar>(QhorLength / vectorLength);
    GRBVar* lambdaN = AllocateDynamicVector<GRBVar>(QhorLength / vectorLength);
    for (int i = 0, inc = 1; i < QhorLength; i += inc)
    {
        xN[i] = CAVmodel->addVar(uMin, 1 * GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x(" + std::to_string(i) + "," + std::to_string(i % vectorLength) + ")");
        epsilonN[i] = CAVmodel->addVar(0, 1 * GRB_INFINITY, 0.0, GRB_CONTINUOUS, "Epsilon=" + std::to_string(i));
        MuN[i] = CAVmodel->addVar(-1 * GRB_INFINITY, 0, 0.0, GRB_CONTINUOUS, "MuN=" + std::to_string(i));
        aN[i] = CAVmodel->addVar(DECCLimit, ACCLimit, 0.0, GRB_CONTINUOUS, "u(" + std::to_string(i / vectorLength) + ")");
    }
    for (int i = 0; i < QhorLength / vectorLength; i++)
    {
        thetaN[i] = CAVmodel->addVar(-1 * GRB_INFINITY, 0, 0.0, GRB_CONTINUOUS, "Theta=" + std::to_string(i));
        lambdaN[i] = CAVmodel->addVar(0, 1 * GRB_INFINITY, 0.0, GRB_CONTINUOUS, "Lambda=" + std::to_string(i));
    }
    for (int i = vectorLength + 1, inc = 2; i < QhorLength; i += inc)
    {
        xN[i].set(GRB_StringAttr_VarName, "vAll(" + std::to_string(i) + "," + std::to_string(i % vectorLength) + ")");
        xN[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
        xN[i].set(GRB_DoubleAttr_LB, uMin);
        xN[i].set(GRB_DoubleAttr_UB, speedLimits[((i % vectorLength) - 1) / 2]);
    }
    for (int i = vectorLength + 1, inc = vectorLength; i < QhorLength; i += inc)
    {
        xN[i].set(GRB_StringAttr_VarName, "vLeader(" + std::to_string(i) + ")");
        xN[i].set(GRB_DoubleAttr_LB, uMin);
        xN[i].set(GRB_DoubleAttr_UB, 1 * GRB_INFINITY);
    }

    if (vehicles_X_V[1] > finalVelocity)
    {
        double uni_t = (2 * displacement) / (finalVelocity + vehicles_X_V[1]);
        if (uni_t > 0 && uni_t < 10.0)
        {
            double uni_a = (finalVelocity - vehicles_X_V[1]) / uni_t;
            for (int m = 0; m < QhorLength; m += vectorLength)
                aN[m].set(GRB_DoubleAttr_UB, uni_a);
        }
    }



    for (int i = 1; i < disturbanceSize; i++)
        xN[(i)*vectorLength].set(GRB_DoubleAttr_UB, trjHOV[i][1] - minimumDisturbanceDistance);
    int index = 0;
    if (vehicleCounts > 1)
    {
        getSystemMatrix(matrixA, simulationTimeStep, &vectorLength);
        fillB(matrixB, &vectorLength, simulationTimeStep);
        //set initial conditions for time = 0
        for (int i = 0, inc = 2; i < vectorLength; i += inc)
        {
            CAVmodel->addConstr(xN[i], '=', vehicles_X_V[i], "LTI-Position(index=" + std::to_string(index) + ",Time=" + std::to_string(0) + ")"); //new location and speed
            CAVmodel->addConstr(xN[i + 1], '=', vehicles_X_V[i + 1], "LTI-Speed(index=" + std::to_string(index) + ",Time=" + std::to_string(0) + ")"); //new location and speed
        }
        //set system constraints when time > 0
        for (int m = 0, inc = vectorLength; m < QhorLength - vectorLength; m += inc)
        {
            {
                GRBLinExpr Aexpr = 0;
                index = m % vectorLength;
                for (int j = 0; j < vectorLength; j++)
                    Aexpr += xN[m + j] * matrixA[index][j];
                Aexpr += aN[m] * 0.5 * simulationTimeStep * simulationTimeStep;
                CAVmodel->addConstr(xN[m + vectorLength], '=', Aexpr, "LTI-LeaderPos(index=" + std::to_string(index) + ",Time=" + std::to_string((m + vectorLength) / vectorLength) + ")"); //new location and speed
            }
            for (int i = m + 1, inc = 2; i < m + vectorLength; i += inc)
            {
                GRBLinExpr Aexpr = 0;
                index = i % vectorLength;
                for (int j = 0; j < vectorLength; j++)
                    Aexpr += xN[m + j] * matrixA[index][j];
                Aexpr += aN[m] * matrixB[index];
                CAVmodel->addConstr(xN[i + vectorLength], '=', Aexpr, "LTI-Speed(index=" + std::to_string(index) + ",Time=" + std::to_string((m + vectorLength) / vectorLength) + ")"); //new location and speed
            }
            for (int i = m + 2, inc = 2; i < m + vectorLength; i += inc)
            {
                GRBLinExpr Aexpr = 0;
                index = i % vectorLength;
                for (int j = 0; j < vectorLength; j++)
                    Aexpr += xN[m + j] * matrixA[index][j];
                CAVmodel->addConstr(xN[i + vectorLength], '=', Aexpr, "LTI-FollowersPos(index=" + std::to_string(index) + ",Time=" + std::to_string((m + vectorLength) / vectorLength) + ")"); //new location and speed
            }
        }
        //Theta constraint
        for (int m = vectorLength, inc = vectorLength; m < QhorLength; m += inc)
        {
            CAVmodel->addConstr((trjHOV[m / vectorLength][1] - xN[m]), GRB_GREATER_EQUAL, (xN[m + 1] * 1.2) + thetaN[m / vectorLength], "timeHeadwayHOVandCAV" + std::to_string(m));//CAV-HOV-S-Limit
            CAVmodel->addConstr(xN[m + 1], GRB_LESS_EQUAL, uMax + lambdaN[m / vectorLength], "speedLimitSoftConstraint" + std::to_string(m));//CAV-HOV-S-Limit

        }
        //Epsilon and Mu constraints
        for (int m = vectorLength, inc = vectorLength; m < QhorLength; m += inc)
            for (int i = m + 3, inc2 = 2; i < m + vectorLength; i += inc2)
            {
                CAVmodel->addConstr(xN[i] - xN[i - vectorLength], GRB_LESS_EQUAL, epsilonN[i] + (ACCLimit)*simulationTimeStep, "acc0");//acc limit
                CAVmodel->addConstr(xN[i] - xN[i - vectorLength], GRB_GREATER_EQUAL, MuN[i] + (DECCLimit)*simulationTimeStep, "decc0");//decc limit
            }
    }
    //Objective function
    GRBQuadExpr qObjFunc = 0;
    //t-tHor:V(0)-Vmax
    for (int tT = vectorLength + 1; tT < QhorLength; tT += vectorLength)
        qObjFunc += (xN[tT] - uMax) * (xN[tT] - uMax);
    //t-tHor:gama*theta(t)^2
    for (int tT = 1; tT < QhorLength / vectorLength; tT++)
    {
        qObjFunc += gama * (thetaN[tT] * thetaN[tT]);
        qObjFunc += gama * (lambdaN[tT - 1] * lambdaN[tT - 1]);
    }

    //qObjFunc += gama * (thetaN[tT] * thetaN[tT]);
//t-(tHor - 1):gama*a_i(t)^2+gama*epsilon_i(t)^2+gama*Mu_i(t)^2
    for (int tT = vectorLength + 3; tT < QhorLength - vectorLength; tT += vectorLength)
    {
        for (int tTi = tT; tTi < tT + vectorLength; tTi += 2)
        {
            qObjFunc += alfaobj * (xN[tTi] - xN[tTi - vectorLength]) * (xN[tTi] - xN[tTi - vectorLength]);
            qObjFunc += (epsilonN[tTi] * epsilonN[tTi]) * gama;
            qObjFunc += (MuN[tTi] * MuN[tTi]) * gama;
        }
    }
    CAVmodel->setObjective(qObjFunc, GRB_MINIMIZE);
    CAVmodel->optimize();

    if (CAVmodel->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        fillOutput2(xN, QhorLength, vectorLength);
        //cleanMemory(CAVmodel);
        DeallocateDynamicArray<double>(matrixA, vectorLength, vectorLength);
        DeallocateDynamicVector<double>(matrixB, vectorLength);
        DeallocateDynamicVector<GRBVar>(xN, QhorLength);
        DeallocateDynamicVector<GRBVar>(epsilonN, QhorLength);
        DeallocateDynamicVector<GRBVar>(MuN, QhorLength);
        DeallocateDynamicVector<GRBVar>(aN, QhorLength);
        DeallocateDynamicVector<GRBVar>(thetaN, QhorLength / vectorLength);
        DeallocateDynamicVector<GRBVar>(lambdaN, QhorLength / vectorLength);
        DeallocateDynamicArray<double>(trjHOV,tHor, 3);


        int counters = CAVmodel->get(GRB_IntAttr_NumVars);
        for (int i = 0; i < counters; i++)
            CAVmodel->remove(CAVmodel->getVar(0));
        counters = CAVmodel->get(GRB_IntAttr_NumConstrs);
        for (int i = 0; i < counters; i++)
            CAVmodel->remove(CAVmodel->getConstr(0));
        counters = CAVmodel->get(GRB_IntAttr_NumQConstrs);
        for (int i = 0; i < counters; i++)
            CAVmodel->remove(CAVmodel->getQConstrs()[0]);
        CAVmodel->reset(0);
        delete CAVmodel;
        delete CAVenv;
        return true;
    }
    else
    {
        //cleanMemory(CAVmodel);
        DeallocateDynamicArray<double>(matrixA, vectorLength, vectorLength);
        DeallocateDynamicVector<double>(matrixB, vectorLength);
        DeallocateDynamicVector<GRBVar>(xN, QhorLength);
        DeallocateDynamicVector<GRBVar>(epsilonN, QhorLength);
        DeallocateDynamicVector<GRBVar>(MuN, QhorLength);
        DeallocateDynamicVector<GRBVar>(aN, QhorLength);
        DeallocateDynamicVector<GRBVar>(thetaN, QhorLength / vectorLength);
        DeallocateDynamicVector<GRBVar>(lambdaN, QhorLength / vectorLength);
        DeallocateDynamicArray<double>(trjHOV, tHor, 3);

        int counters = CAVmodel->get(GRB_IntAttr_NumVars);
        for (int i = 0; i < counters; i++)
            CAVmodel->remove(CAVmodel->getVar(0));
        counters = CAVmodel->get(GRB_IntAttr_NumConstrs);
        for (int i = 0; i < counters; i++)
            CAVmodel->remove(CAVmodel->getConstr(0));
        counters = CAVmodel->get(GRB_IntAttr_NumQConstrs);
        for (int i = 0; i < counters; i++)
            CAVmodel->remove(CAVmodel->getQConstrs()[0]);
        CAVmodel->reset(0);
        delete CAVmodel;
        delete CAVenv;
        return false;
    }
}



void MPC::fillOutput2(GRBVar* CAVs, int matrixLength, int vectorLength)
{
    optimizedCAVs = AllocateDynamicVector<double>(vectorLength);
    optimizedCAVsACC = AllocateDynamicVector<double>(vectorLength / 2);
    double acc = 0;
    int j = 1;
    acc = CAVs[j + vectorLength].get(GRB_DoubleAttr_X) - CAVs[j].get(GRB_DoubleAttr_X);
    optimizedCAVsACC[(j - 1) / 2] = acc;
    if (acc > 0.3 && acc < -0.4)
        bool error = true;

    for (j = 3; j < vectorLength; j += 2)
    {
        acc = CAVs[j + vectorLength].get(GRB_DoubleAttr_X) - CAVs[j].get(GRB_DoubleAttr_X);
        acc = max(min(ACCLimit * 1.05 * simulationTimeStep, acc), DECCLimit * 1.05 * simulationTimeStep);
        optimizedCAVsACC[(j - 1) / 2] = acc;
        if (acc > 0.3 && acc < -0.4)
            bool error = true;
    }
    for (int j = 0; j < vectorLength; j++)
        optimizedCAVs[j] = CAVs[j + vectorLength].get(GRB_DoubleAttr_X);
}

void MPC::getOptimizedTrajectories(double* Matrix)
{
    for (int j = 0; j < vehicleCounts * 2; j++)
        Matrix[j] = optimizedCAVs[j];
}
void MPC::getOptimizedACCs(double* Matrix)
{
    for (int j = 0; j < vehicleCounts; j++)
        Matrix[j] = optimizedCAVsACC[j];
}

void MPC::setDisturbanceTrajectories(double** disturbanceTrj, double* HDV_X_V, const int HOVsize, const int initialTime, int disturbanceSize)
{

    double s1, b1, st, bt, ft, prevst, speed;
    s1 = HDV_X_V[1];
    b1 = HDV_X_V[3] - HDV_X_V[1];
    st = (exponentialAlpha * HDV_X_V[3]) + ((1 - exponentialAlpha) * (s1 + b1));
    bt = (exponentialBeta * (st - s1)) + ((1 - exponentialBeta) * b1);

    for (int i = 5, inc = 2; i < HOVsize * 2; i += inc)
    {
        prevst = st;
        st = (exponentialAlpha * HDV_X_V[i]) + ((1 - exponentialAlpha) * (st + bt));
        bt = (exponentialBeta * (st - prevst)) + ((1 - exponentialBeta) * bt);
    }

    speed = max(0.0, st + bt);
    disturbanceTrj[0][0] = initialTime;
    disturbanceTrj[0][1] = HDV_X_V[(HOVsize * 2) - 2]; //x(t+1) = x(t) + Tv(t); //x(0)
    disturbanceTrj[0][2] = HDV_X_V[(HOVsize * 2) - 1]; //v(0)

    for (int i = 1; i < disturbanceSize; i++)
    {
        disturbanceTrj[i][0] = i;;
        disturbanceTrj[i][1] = disturbanceTrj[i - 1][1] + simulationTimeStep * speed; //x(t+1) = x(t) + Tv(t); //x(0)
        disturbanceTrj[i][2] = speed;
    }
}
void MPC::fillB(double* matrix, const int* vectorLength, double timeStep)
{
    for (int i = 0, inc = 2; i < *vectorLength; i += inc)
    {
        matrix[i] = 0;
    }
    if (*vectorLength < 3)
    {
        matrix[1] = timeStep;
    }
    if (*vectorLength > 3)
    {
        double square = 1;
        for (int i = 1, inc = 2; i < *vectorLength; i += inc)
        {
            matrix[i] = timeStep * pow(CACC_kA, square);
            square++;
        }
    }
}
void MPC::fillB_2(double** matrix, const int* vectorLength, double timeStep)
{
    //filling leader ACC part
    for (int i = 0, inc = 2; i < *vectorLength; i += inc)
    {
        matrix[0][i] = 0;
    }
    if (*vectorLength < 3)
    {
        matrix[0][1] = timeStep;
    }
    if (*vectorLength > 3)
    {
        double square = 1;
        for (int i = 1, inc = 2; i < *vectorLength; i += inc)
        {
            matrix[0][i] = timeStep * pow(CACC_kA, square);
            square++;
        }
    }
    //filling dn(t) part
    for (int i = 1; i <= *vectorLength; i++)
        for (int j = 0; j < *vectorLength; j++)
            matrix[i][j] = 0;
    for (int i = 3, inc = 2; i < *vectorLength; i += inc)
    {
        matrix[i + 1][i] = -1 * timeStep * CACC_kC;
        for (int j = i + 2, inc = 2; j < *vectorLength; j += inc)
            matrix[j + 1][i] = matrix[j + 1 - 2][i] * CACC_kA;
    }
}

void MPC::getSystemMatrix(double** matrix, double timeStep, const int* matrixSize)
{
    double** matrixA = AllocateDynamicArray<double>(10, 10);
    fillRef5vehicles(matrixA, simulationTimeStep);
    if (*matrixSize <= 6)
    {
        for (int i = 0; i < *matrixSize; i++)
        {
            for (int j = 0; j < *matrixSize; j++)
            {
                matrix[i][j] = matrixA[i][j];
            }
        }
    }
    if (*matrixSize > 6)
    {
        for (int i = 0; i < 6; i++)//fill by reference table in new table
        {
            for (int j = 0; j < 6; j++)
            {
                matrix[i][j] = matrixA[i][j];
            }
        }
        for (int i = 0; i < 6; i++)//zeros for rest of columns of reference table
        {
            for (int j = 6; j < *matrixSize; j++)
            {
                matrix[i][j] = 0;
            }
        }
        for (int i = 6; i < *matrixSize; i++)
        {
            if ((i % 2) == 0)//fill the row of x_i(t)
            {
                for (int j = 0; j < *matrixSize; j++)
                    matrix[i][j] = 0;
                matrix[i][i] = 1;
                matrix[i][i + 1] = simulationTimeStep;
            }
            else//fill the row of v_i(t)
            {
                for (int j = i + 1; j < *matrixSize; j++)//zero for rest columns
                {
                    matrix[i][j] = 0;
                }
                for (int j = 0; j <= i - 4; j++)
                {
                    matrix[i][j] = matrix[i - 2][j] * CACC_kA;
                }
                matrix[i][i - 3] = (matrix[i - 2][i - 3] * CACC_kA) + (timeStep * CACC_kC);
                matrix[i][i - 2] = ((matrix[i - 2][i - 2] - 1) * CACC_kA) + (timeStep * CACC_kB);
                matrix[i][i - 1] = -1 * timeStep * CACC_kC;
                matrix[i][i] = 1 - timeStep * (CACC_kB + (CACC_kC * CACC_H));
                //for (int j = i; j > i - 4; j--)
                //{
                //    matrix[i][j] = matrix[i - 2][j - 2];
                //}
            }
        }
    }
    DeallocateDynamicArray<double>(matrixA,10, 10);
}
void MPC::fillRef5vehicles(double** matrix, double timeStep)
{

    matrix[0][0] = 1;
    matrix[0][1] = timeStep;
    matrix[0][2] = 0;
    matrix[0][3] = 0;
    matrix[0][4] = 0;
    matrix[0][5] = 0;
    matrix[0][6] = 0;
    matrix[0][7] = 0;
    matrix[0][8] = 0;
    matrix[0][9] = 0;

    matrix[1][0] = 0;
    matrix[1][1] = 1;
    matrix[1][2] = 0;
    matrix[1][3] = 0;
    matrix[1][4] = 0;
    matrix[1][5] = 0;
    matrix[1][6] = 0;
    matrix[1][7] = 0;
    matrix[1][8] = 0;
    matrix[1][9] = 0;

    matrix[2][0] = 0;
    matrix[2][1] = 0;
    matrix[2][2] = 1;
    matrix[2][3] = timeStep;
    matrix[2][4] = 0;
    matrix[2][5] = 0;
    matrix[2][6] = 0;
    matrix[2][7] = 0;
    matrix[2][8] = 0;
    matrix[2][9] = 0;

    matrix[3][0] = timeStep * CACC_kC;
    matrix[3][1] = timeStep * CACC_kB;
    matrix[3][2] = -1 * timeStep * CACC_kC;
    matrix[3][3] = 1 - (timeStep * (CACC_kB + CACC_kC * CACC_H));
    matrix[3][4] = 0;
    matrix[3][5] = 0;
    matrix[3][6] = 0;
    matrix[3][7] = 0;
    matrix[3][8] = 0;
    matrix[3][9] = 0;

    matrix[4][0] = 0;
    matrix[4][1] = 0;
    matrix[4][2] = 0;
    matrix[4][3] = 0;
    matrix[4][4] = 1;
    matrix[4][5] = timeStep;
    matrix[4][6] = 0;
    matrix[4][7] = 0;
    matrix[4][8] = 0;
    matrix[4][9] = 0;

    matrix[5][0] = (timeStep * CACC_kA * CACC_kC);
    matrix[5][1] = (timeStep * CACC_kA * CACC_kB);
    matrix[5][2] = timeStep * (CACC_kC - (CACC_kA * CACC_kC));
    matrix[5][3] = timeStep * (CACC_kB - (CACC_kA * CACC_kB) - (CACC_kA * CACC_kC * CACC_H));
    matrix[5][4] = -1 * (timeStep * CACC_kC);
    matrix[5][5] = 1 - timeStep * (CACC_kB + (CACC_kC * CACC_H));
    matrix[5][6] = 0;
    matrix[5][7] = 0;
    matrix[5][8] = 0;
    matrix[5][9] = 0;

    matrix[6][0] = 0;
    matrix[6][1] = 0;
    matrix[6][2] = 0;
    matrix[6][3] = 0;
    matrix[6][4] = 0;
    matrix[6][5] = 0;
    matrix[6][6] = 1;
    matrix[6][7] = timeStep;
    matrix[6][8] = 0;
    matrix[6][9] = 0;

    matrix[7][0] = timeStep * CACC_kA * CACC_kA * CACC_kC;//x1
    matrix[7][1] = timeStep * CACC_kA * CACC_kA * CACC_kB;//v1
    matrix[7][2] = timeStep * CACC_kA * (CACC_kC - (CACC_kA * CACC_kC));//x2
    matrix[7][3] = timeStep * CACC_kA * (CACC_kB - (CACC_kA * CACC_kB) - (CACC_kA * CACC_kC * CACC_H));//v2
    matrix[7][4] = timeStep * (CACC_kC - (CACC_kA * CACC_kC));//x3
    matrix[7][5] = timeStep * (CACC_kB - (CACC_kA * CACC_kB) - (CACC_kA * CACC_kC * CACC_H));//v3
    matrix[7][6] = -1 * timeStep * CACC_kC;//x4
    matrix[7][7] = 1 - timeStep * (CACC_kB + CACC_kC * CACC_H);//v4
    matrix[7][8] = 0;
    matrix[7][9] = 0;

    matrix[8][0] = 0;
    matrix[8][1] = 0;
    matrix[8][2] = 0;
    matrix[8][3] = 0;
    matrix[8][4] = 0;
    matrix[8][5] = 0;
    matrix[8][6] = 0;
    matrix[8][7] = 0;
    matrix[8][8] = 1;
    matrix[8][9] = timeStep;

    matrix[9][0] = timeStep * CACC_kA * CACC_kA * CACC_kA * CACC_kC;
    matrix[9][1] = timeStep * CACC_kA * CACC_kA * CACC_kA * CACC_kB;
    matrix[9][2] = timeStep * CACC_kA * CACC_kA * (CACC_kC - (CACC_kA * CACC_kC));
    matrix[9][3] = timeStep * CACC_kA * CACC_kA * (CACC_kB - (CACC_kA * CACC_kB) - (CACC_kA * CACC_kC * CACC_H));
    matrix[9][4] = timeStep * CACC_kA * (CACC_kC - (CACC_kA * CACC_kC));
    matrix[9][5] = timeStep * CACC_kA * (CACC_kB - (CACC_kA * CACC_kB) - (CACC_kA * CACC_kC * CACC_H));
    matrix[9][6] = timeStep * (CACC_kC - (CACC_kA * CACC_kC));
    matrix[9][7] = timeStep * (CACC_kB - (CACC_kA * CACC_kB) - (CACC_kA * CACC_kC * CACC_H));
    matrix[9][8] = -1 * timeStep * CACC_kC;
    matrix[9][9] = 1 - timeStep * (CACC_kB + (CACC_kC * CACC_H));
    //cout << "t0 = " << timeStep << " q[" << 9 << "][" << 9<< "]" << matrix[9][9] << endl;
}
void MPC::getSystemMatrix_2(double** matrix, double timeStep, const int* matrixSize)
{
    double** matrixA = AllocateDynamicArray<double>(10, 10);
    fillRef5vehicles_2(matrixA, simulationTimeStep);
    for (int i = 0; i < *matrixSize; i++)
    {
        for (int j = 0; j < *matrixSize; j++)
        {
            matrix[i][j] = matrixA[i][j];
        }
    }
}
void MPC::fillRef5vehicles_2(double** matrix, double timeStep)
{

    matrix[0][0] = 1;
    matrix[0][1] = timeStep;
    matrix[0][2] = 0;
    matrix[0][3] = 0;
    matrix[0][4] = 0;
    matrix[0][5] = 0;
    matrix[0][6] = 0;
    matrix[0][7] = 0;
    matrix[0][8] = 0;
    matrix[0][9] = 0;

    matrix[1][0] = 0;
    matrix[1][1] = 1;
    matrix[1][2] = 0;
    matrix[1][3] = 0;
    matrix[1][4] = 0;
    matrix[1][5] = 0;
    matrix[1][6] = 0;
    matrix[1][7] = 0;
    matrix[1][8] = 0;
    matrix[1][9] = 0;

    matrix[2][0] = 0;
    matrix[2][1] = 0;
    matrix[2][2] = 1;
    matrix[2][3] = timeStep;
    matrix[2][4] = 0;
    matrix[2][5] = 0;
    matrix[2][6] = 0;
    matrix[2][7] = 0;
    matrix[2][8] = 0;
    matrix[2][9] = 0;

    matrix[3][0] = timeStep * CACC_kC;
    matrix[3][1] = timeStep * CACC_kB;
    matrix[3][2] = -1 * timeStep * CACC_kC;
    matrix[3][3] = 1 - (timeStep * (CACC_kB));
    matrix[3][4] = 0;
    matrix[3][5] = 0;
    matrix[3][6] = 0;
    matrix[3][7] = 0;
    matrix[3][8] = 0;
    matrix[3][9] = 0;

    matrix[4][0] = 0;
    matrix[4][1] = 0;
    matrix[4][2] = 0;
    matrix[4][3] = 0;
    matrix[4][4] = 1;
    matrix[4][5] = timeStep;
    matrix[4][6] = 0;
    matrix[4][7] = 0;
    matrix[4][8] = 0;
    matrix[4][9] = 0;

    matrix[5][0] = (timeStep * CACC_kA * CACC_kC);
    matrix[5][1] = (timeStep * CACC_kA * CACC_kB);
    matrix[5][2] = timeStep * (CACC_kC - (CACC_kA * CACC_kC));
    matrix[5][3] = timeStep * (CACC_kB - (CACC_kA * CACC_kB));
    matrix[5][4] = -1 * (timeStep * CACC_kC);
    matrix[5][5] = 1 - timeStep * (CACC_kB);
    matrix[5][6] = 0;
    matrix[5][7] = 0;
    matrix[5][8] = 0;
    matrix[5][9] = 0;

    matrix[6][0] = 0;
    matrix[6][1] = 0;
    matrix[6][2] = 0;
    matrix[6][3] = 0;
    matrix[6][4] = 0;
    matrix[6][5] = 0;
    matrix[6][6] = 1;
    matrix[6][7] = timeStep;
    matrix[6][8] = 0;
    matrix[6][9] = 0;

    matrix[7][0] = timeStep * CACC_kA * CACC_kA * CACC_kC;
    matrix[7][1] = timeStep * CACC_kA * CACC_kA * CACC_kB;
    matrix[7][2] = timeStep * CACC_kA * (CACC_kC - (CACC_kA * CACC_kC));
    matrix[7][3] = timeStep * CACC_kA * (CACC_kB - (CACC_kA * CACC_kB));
    matrix[7][4] = timeStep * (CACC_kC - (CACC_kA * CACC_kC));
    matrix[7][5] = timeStep * (CACC_kB - (CACC_kA * CACC_kB));
    matrix[7][6] = -1 * timeStep * CACC_kC;
    matrix[7][7] = 1 - timeStep * (CACC_kB);
    matrix[7][8] = 0;
    matrix[7][9] = 0;

    matrix[8][0] = 0;
    matrix[8][1] = 0;
    matrix[8][2] = 0;
    matrix[8][3] = 0;
    matrix[8][4] = 0;
    matrix[8][5] = 0;
    matrix[8][6] = 0;
    matrix[8][7] = 0;
    matrix[8][8] = 1;
    matrix[8][9] = timeStep;

    matrix[9][0] = timeStep * CACC_kA * CACC_kA * CACC_kA * CACC_kC;
    matrix[9][1] = timeStep * CACC_kA * CACC_kA * CACC_kA * CACC_kB;
    matrix[9][2] = timeStep * CACC_kA * CACC_kA * (CACC_kC - (CACC_kA * CACC_kC));
    matrix[9][3] = timeStep * CACC_kA * CACC_kA * (CACC_kB - (CACC_kA * CACC_kB));
    matrix[9][4] = timeStep * CACC_kA * (CACC_kC - (CACC_kA * CACC_kC));
    matrix[9][5] = timeStep * CACC_kA * (CACC_kB - (CACC_kA * CACC_kB));
    matrix[9][6] = timeStep * (CACC_kC - (CACC_kA * CACC_kC));
    matrix[9][7] = timeStep * (CACC_kB - (CACC_kA * CACC_kB));
    matrix[9][8] = -1 * timeStep * CACC_kC;
    matrix[9][9] = 1 - timeStep * (CACC_kB);
    //cout << "t0 = " << timeStep << " q[" << 9 << "][" << 9<< "]" << matrix[9][9] << endl;
}
double MPC::getR()
{
    double tmpTerm = 0.0;
    for (int i = 1; i <= vehicleCounts; i++)
        tmpTerm += pow(CACC_kA, (i - 1.0) * 2.0);
    tmpTerm = tmpTerm * 2.0 * (double)alfaobj;
    return tmpTerm;
}
void MPC::writeTrajectories(GRBVar* CAVs, int matrixLength, int vectorLength, std::string file_name)// CSV file
{
    cout << "now writing files" << endl;
    ofstream FILE;
    //writing trajectories
    {
        std::string fileName_ = file_name + "_Trajectories.csv";
        FILE.open(fileName_);
        int t_ = 0;
        int loopSize = 0;
        int countTimer = 0;
        for (int i = 0, inc = vectorLength; i < matrixLength; i += vectorLength)
        {
            loopSize = i + vectorLength;

            FILE << countTimer * simulationTimeStep << comaType;
            countTimer++;
            //location
            for (int m = i, inc2 = 2; m < loopSize; m += inc2)
            {
                FILE << CAVs[m].get(GRB_DoubleAttr_X) << comaType;
            }
            //speed
            for (int m = i + 1, inc2 = 2; m < loopSize; m += inc2)
            {
                FILE << CAVs[m].get(GRB_DoubleAttr_X) << comaType;
            }
            for (int m = i + 1, inc2 = 2; m < loopSize; m += inc2)
            {
                if (i > 0)
                    FILE << (CAVs[m].get(GRB_DoubleAttr_X) - CAVs[m - vectorLength].get(GRB_DoubleAttr_X)) / (double)simulationTimeStep << comaType;
                else
                    FILE << 0 << comaType;

            }
            FILE << endl;
        }
        FILE << endl;
        FILE.close();
    }
}

//vehicleCounts
int MPC::getvehicleCounts() {
    return vehicleCounts;
}

//simulationTimeStep
void MPC::setfinalVelocity(double s) {
    finalVelocity = s;
}
double MPC::getfinalVelocity() {
    return finalVelocity;
}
//simulationTimeStep
void MPC::setdisplacement(double s) {
    displacement = s;
}
double MPC::getdisplacement() {
    return displacement;
}
//simulationTimeStep
void MPC::setsimulationTimeStep(double s) {
    simulationTimeStep = s;
}
double MPC::getsimulationTimeStep() {
    return simulationTimeStep;
}
//simulationStartTime
void MPC::setsimulationStartTime(double s) {
    simulationStartTime = s;
}
double MPC::getsimulationStartTime() {
    return simulationStartTime;
}
//minimunDistance
void MPC::setminimunDistance(double s) {
    minimunDistance = s;
}
double MPC::getminimunDistance() {
    return minimunDistance;
}
//minimunDisturbanceDistance
void MPC::setminimunDisturbanceDistance(double s)
{
    minimumDisturbanceDistance = s;
}
double MPC::getminimunDisturbanceDistance()
{
    return minimumDisturbanceDistance;
}
//CACC_k
void MPC::setCACC_k(double s) {
    CACC_k = s;
}
double MPC::getCACC_k() {
    return CACC_k;
}
//CACC_kA
void MPC::setCACC_kA(double s) {
    CACC_kA = s;
}
double MPC::getCACC_kA() {
    return CACC_kA;
}
//CACC_kB
void MPC::setCACC_kB(double s) {
    CACC_kB = s;
}
double MPC::getCACC_kB() {
    return CACC_kB;
}
//CACC_kC
void MPC::setCACC_kC(double s) {
    CACC_kC = s;
}
double MPC::getCACC_kC() {
    return CACC_kC;
}
//CACC_H
void MPC::setCACC_H(double s) {
    CACC_H = s;
}
double MPC::getCACC_H() {
    return CACC_H;
}
//ACCLimit
void MPC::setACCLimit(double s) {
    ACCLimit = s;
}
double MPC::getACCLimit() {
    return ACCLimit;
}
//DECCLimit
void MPC::setDECCLimit(double s) {
    DECCLimit = s;
}
double MPC::getDECCLimit() {
    return DECCLimit;
}
//x0Acc
void MPC::setx0Acc(double s) {
    x0Acc = s;
}
double MPC::getx0Acc() {
    return x0Acc;
}
//xMin
void MPC::setxMin(double s) {
    xMin = s;
}
double MPC::getxMin() {
    return xMin;
}
//xMax
void MPC::setxMax(double s) {
    xMax = s;
}
double MPC::getxMax() {
    return xMax;
}
//uMin
void MPC::setuMin(double s) {
    uMin = s;
}
double MPC::getuMin() {
    return uMin;
}
//uMax
void MPC::setuMax(double s) {
    uMax = s;
}
double MPC::getuMax() {
    return uMax;
}
//epsilonMin
void MPC::setepsilonMin(double s) {
    epsilonMin = s;
}
double MPC::getepsilonMin() {
    return epsilonMin;
}
//epsilonMax
void MPC::setepsilonMax(double s) {
    epsilonMax = s;
}
double MPC::getepsilonMax() {
    return epsilonMax;
}
//tHor
int MPC::gettHor() {
    return tHor;
}
//alfaobj
void MPC::setalfaobj(double s) {
    alfaobj = s;
}
double MPC::getalfaobj() {
    return alfaobj;
}
//gama
void MPC::setgama(double s) {
    gama = s;
}
double MPC::getgama() {
    return gama;
}
//leaderACC
void MPC::setleaderACC(double s) {
    leaderAcc = s;
}
double MPC::getleaderACC() {
    return leaderAcc;
}