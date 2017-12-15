#define KALMANFILTER_CPP

#include <iostream>
#include <map>

#ifdef USING_OMP_PRAGMAS
    #include <omp.h>
#endif

#include <limits>
#include <set>
#include <iterator>


#include "KalmanFilter.h"

using std::cerr;
using std::endl;


Kalmanfilter::Kalmanfilter()
{
 
}

Kalmanfilter::~Kalmanfilter()
{
}

void Kalmanfilter::init(vpColVector& f,vpColVector& p, vpHomogeneousMatrix &_cMo, vpMatrix& stiffnessMatrix)
{

H.resize(2*f.getRows(),2*f.getRows());
H.setIdentity();

H_0.resize(2*f.getRows(),2*f.getRows());

J.resize(2*f.getRows(),2*f.getRows());
J.setIdentity();

Q.resize(2*f.getRows(),2*f.getRows());
R.resize(2*f.getRows(),2*f.getRows());
PEst.resize(2*f.getRows(),2*f.getRows());
I.resize(2*f.getRows(),2*f.getRows());
I.setIdentity();

KPose.resize(6,6);
KPose.setIdentity();

qT = 0.02;
qR = 0.001;
pT = 0.001;
pR = 0.0001;

std::cout << " init 1 " << std::endl;

estimatedVelPose.resize(6);
QPose.resize(12,12);
RPose.resize(6,6);
JPose.resize(12,12);
HPose.resize(6,12);
H_0Pose.resize(6,12);
PEstPose.resize(12,12);
JPose.setIdentity();
QPose.setIdentity();	
LPose.resize(6,12);
QPose =qT*qT*QPose;
PEstPose.setIdentity();
PEstPose = pT*pT*PEstPose;

	KPose.resize(12,6);

	QVel.resize(6,6);
	RVel.resize(6,6);
	JVel.resize(6,6);
	KVel.resize(6,6);
	KVel.setIdentity();
	PEstVel.resize(6,6);
	JVel.setIdentity();
	QVel.setIdentity();
	RVel.setIdentity();
	QVel =qR*qR*QVel;
	PEstVel.setIdentity();
	PEstVel = pR*pR*PEstVel;


	for(int j = 0; j<12; j++)
		for(int i = 0; i<12; i++)
		{
			//R[i][j] = 0.0001;
			if(i>=6 && j<6)
				{
				if(j==i-6)
					{
					KPose[i][j] = 1;
					}
				}
			if(i<6 && j>=6)
				{
				if(i==j-6)
					{JPose[i][j] = 1;
					H_0Pose[i][j] = 1;
					}

				LPose[i][j] = 1;
				}
			else if(j<6 && i<6)
				{
				QPose[i][j] = 0;
				if(i == j)
				{HPose[i][j] = 1;
				RPose[i][j] = 0.02;
				//RVel[i][j] = 0.00000001;
				RVel[i][j] = 0.001;
				}
				if (j<3 && i<3){
					if(i == j)
						{QVel[i][j] = qT*qT;
						PEstVel[i][j] = pT*pT;
						//R[i][j] = 3e-09;
						RPose[i][j] = 0.001;
						//RVel[i][j] = 0.0000001;
						RVel[i][j] = 0.0001;
						}
				}
				}
		}
	QPose[9][9] = qR*qR;
	QPose[10][10] = qR*qR;
	QPose[11][11] = qR*qR;
	estimatedcMo = _cMo;
	estimatedcMo_0 = _cMo;
	estimatedVelPose[0] = 0;
	estimatedVelPose[1] = 0;
	estimatedVelPose[2] = 0;
	estimatedVelPose[3] = 0;
	estimatedVelPose[4] = 0;
	estimatedVelPose[5] = 0;
	measuredVelPose = estimatedVelPose;

std::cout << " init 3 " << std::endl;

//PEst = R;

qF = 0.001;
rF = 0.0001;

qX = 0.03;
rX = 0.001;

for (int k = 0; k < f.getRows(); k++)
	{
		Q[k][k] = (double)qX;
		Q[k+f.getRows()][k+f.getRows()] = (double)qF;
		R[k][k] = (double)rX;
		R[k+f.getRows()][k+f.getRows()] = (double)rF;
	}	
	
	estimatedForces = f;
	estimatedPositions = p;
	
	estimatedState = estimatedPositions;
	estimatedState.stack(estimatedForces);	
	
    predictedForces = estimatedForces;
	predictedPositions = estimatedPositions;
	}

void Kalmanfilter::predictPose(vpMatrix& stiffnessMatrix)
{
	/*for (int k = 0; k < stiffnessMatrix.getCols(); k++)
	for (int l = 0; l < stiffnessMatrix.getCols(); l++)
	std::cout <<k<<" "<<l<<"   "<< stiffnessMatrix[k][l] << std::endl;*/
	
	predictedVelPose = estimatedVelPose;
	vpMatrix Iv(6,6);
	Iv.setIdentity();
	vpMatrix I(12,12);
	I.setIdentity();
	predictedcMo = vpExponentialMap::direct(predictedVelPose).inverse() * estimatedcMo;
	predictedcMo_0 = vpExponentialMap::direct((Iv)*predictedVelPose).inverse() * estimatedcMo;
	//cMoPred_0 = vpExponentialMap::direct(vMes-vPred).inverse() * cMoEst;
	//cMoPred_0 = cMoEst;
	PPredPose = JPose*PEstPose*JPose.transpose() + QPose;
    PPredVel = PEstVel + QVel;
    estimatedcMo_0 = estimatedcMo;
	
	predictedForces = estimatedForces;
	predictedPositions = estimatedPositions + estimatedForces; //stiffnessMatrix.pseudoInverse()*estimatedForces;
		
	predictedState = predictedPositions;
	predictedState.stack(predictedForces);
	
	for (int k = 0; k < predictedForces.getRows(); k++)
	for (int l = predictedForces.getRows(); l < 2*predictedForces.getRows(); l++)
		J[k][l] = 1;//stiffnessMatrix[k][l-predictedForces.getRows()];
		
	//PPred = J*PEst*J.transpose() + Q;
    //PfPred = PEst + Q;
		

}

void Kalmanfilter::predictPoints(vpColVector& p)
{
	
	predictedForces = estimatedForces;
	estimatedPositions = p;
	predictedPositions = estimatedPositions + estimatedForces; //stiffnessMatrix.pseudoInverse()*estimatedForces;
		
	predictedState = predictedPositions;
	predictedState.stack(predictedForces);
	
	for (int k = 0; k < predictedForces.getRows(); k++)
	for (int l = predictedForces.getRows(); l < 2*predictedForces.getRows(); l++)
		J[k][l] = 1;//stiffnessMatrix[k][l-predictedForces.getRows()];
		
	//PPred = J*PEst*J.transpose() + Q;
    //PfPred = PEst + Q;
		

}


void Kalmanfilter::estimatePose(vpHomogeneousMatrix &_cMo)
{
	
innovationcMo = _cMo;//measuredcMo*(estimatedcMo.inverse());
//K = PPred*((PPred + R).pseudoInverse());
//KPose = PPredPose*HPose.transpose()*((HPose*PPredPose*HPose.transpose() + RPose).pseudoInverse());
//Kv = PvPred*((PvPred + R).pseudoInverse());
KVel = PPredVel*((PPredVel + RVel).pseudoInverse());
vpMatrix Iv(6,6);
Iv.setIdentity();
measuredVelPose = vpExponentialMap::inverse((innovationcMo).inverse());

std::cout << " Kv " << KVel << std::endl;
KVel.setIdentity();
KVel *= 0.7;
estimatedVelPose = predictedVelPose + (KVel)*(measuredVelPose-predictedVelPose);
//vEst = vPred + (H_0)*K*(vMes-vPred);
estimatedcMo = vpExponentialMap::direct(KVel*(measuredVelPose-predictedVelPose)).inverse()*predictedcMo;
//cMoEst = vpExponentialMap::direct(Kv*vpExponentialMap::inverse(cMoInnov.inverse())).inverse()*cMoPred;
//cMoEst = cMoMes;
	
}

void Kalmanfilter::estimatePoints(vpColVector& f, vpColVector& p)
{
	
//K = PPred*((PPred + R).pseudoInverse());

//PEst = (I-K*H)*PPred;
//PEstVel = (Iv-KVel)*PPredVel;	

measuredState = p;
measuredState.stack(f);

//K.setIdentity();
//K *= 0.7;

estimatedState = predictedState + 1*(measuredState - predictedState);
//PEst = (I-K)*PPred;

for (int l = 0; l < predictedForces.getRows(); l++)
{
estimatedPositions[l] = estimatedState[l];
estimatedForces[l] = estimatedState[l+predictedForces.getRows()];
}
//estimatedForces = estimatedState.rows(f.getRows(), 2*f.getRows()-1);

for (int l = 0; l < predictedForces.getRows(); l++){
	//std::cout <<l<<"   "<< predictedForces[l] << " " << estimatedForces[l] << std::endl;
}

//PvEst = (Iv-Kv)*PvPred;	
	
}


