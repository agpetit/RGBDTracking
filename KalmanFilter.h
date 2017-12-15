#ifndef KALMANFILTER_H
#define KALMANFILTER_H


#include <set>

#ifdef WIN32
#include <process.h>
#else
#include <pthread.h>
#endif

#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>

#include <cuda_runtime.h>
#include <npp.h>
#include <nppi.h>


#include <visp/vpIoTools.h>
#include <visp/vpImageIo.h>
#include <visp/vpParseArgv.h>
#include <visp/vpMatrix.h>
#include <visp/vpPose.h>
#include <visp/vpExponentialMap.h>
#include <visp/vpMatrixException.h>
#include <visp/vpImagePoint.h>

#include <iostream>
#include <string>
#include <map>
#include <boost/thread.hpp>
#include <sys/times.h>

using namespace std;

class Kalmanfilter
{
	
public:
    Kalmanfilter();
    virtual ~Kalmanfilter();
	
	void init(vpColVector& f,vpColVector& p, vpHomogeneousMatrix &_cMo, vpMatrix& stiffnessMatrix);
    void predictPose(vpMatrix& stiffnessMatrix);
    void estimatePose(vpHomogeneousMatrix &_cMo);
	void estimatePoints(vpColVector& f,vpColVector& p);
	void predictPoints(vpColVector& p);
	void convert(vpHomogeneousMatrix &_cMo, vpColVector &poseQuat_);
	void convert(vpColVector &poseQuat_, vpHomogeneousMatrix &_cMo);
	void computeJacobian(double dt);
	vpMatrix getTheta(vpColVector &w);
	vpMatrix computeExp(vpMatrix &Theta, vpColVector &w);
	vpMatrix computePhi(vpColVector &w);
	
	double qX;
    double qF;
	double rX;
    double rF;
	
	double qT;
    double qR;
	double rT;
    double rR;
	double pT;
    double pR;
	
    // kalman vectors
	vpColVector estimatedState;
	vpColVector predictedState;
	vpColVector measuredState;
	vpColVector estimatedPositions;
	vpColVector predictedPositions;
	vpColVector estimatedForces;
	vpColVector predictedForces;
	vpColVector predictedVelPose;
	vpColVector estimatedVelPose;
	vpColVector measuredVelPose;
	vpColVector estimatedPose;
	vpColVector predictedPose;
	vpColVector measuredPose;
	vpColVector estimatedStatePose;
	vpColVector predictedStatePose;
	vpColVector measuredStatePose;
	
	vpHomogeneousMatrix predictedcMo;
	vpHomogeneousMatrix predictedcMo_0;
	vpHomogeneousMatrix estimatedcMo;
	vpHomogeneousMatrix estimatedcMo_0;
	vpHomogeneousMatrix innovationcMo;
	vpHomogeneousMatrix measuredcMo;

	vpMatrix PEst;
	vpMatrix PPred;
	vpMatrix PEstVel;
	vpMatrix PPredVel;
	vpMatrix PfEst;
	vpMatrix PfPred;
	vpMatrix Q;
	vpMatrix R;
	vpMatrix J;
	vpMatrix K;
	vpMatrix H;
    vpMatrix I;
	vpMatrix H_0;
	vpMatrix LPose;
	
	vpMatrix PEstPose;
	vpMatrix PPredPose;
	vpMatrix PfEstPose;
	vpMatrix PfPredPose;
	vpMatrix QPose;
	vpMatrix QVel;
	vpMatrix RPose;
	vpMatrix RVel;
	vpMatrix JPose;
	vpMatrix JVel;
	vpMatrix KPose;
	vpMatrix KVel;
	vpMatrix HPose;
    vpMatrix IPose;
	vpMatrix H_0Pose;
	
};

#endif
