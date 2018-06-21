/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#ifndef SOFA_RGBDTRACKING_REGISTRATIONRIGID_H
#define SOFA_RGBDTRACKING_REGISTRATIONRIGID_H

#include <RGBDTracking/config.h>
#include <image/ImageTypes.h>
#include <sofa/core/core.h>

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>

#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>

#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <visp/vpKltOpencv.h>
#include <SofaGeneralEngine/NormalsFromPoints.h>
//#include <sofa/helper/kdTree.inl>
#include "KalmanFilter.h"
#include <visp/vpDisplayX.h>
#include <algorithm>    
#ifdef WIN32
#include <process.h>
#else
#include <pthread.h>
#endif

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <sys/times.h>

#include <opencv/cv.h>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <visp/vpIoTools.h>
#include <visp/vpImageIo.h>
#include <visp/vpParseArgv.h>
#include <visp/vpMatrix.h>

#include <string>
#include <boost/thread.hpp>
#include "ccd.h"
#include "RGBDDataProcessing.h"
#include "MeshProcessing.h"
#include "ImageConverter.h"


using namespace std;
using namespace cv;

typedef struct {
  // for softassign
  sofa::defaulttype::Vector3 coef; 
  int triangle;   // parameter for outliers (see the original Softassign paper). default: 3.0
   
} mapping;

typedef struct point_struct{

  double x;
  double y; 
  double z;
}point_struct;

namespace sofa
{

namespace core
{


using helper::vector;
using namespace sofa::defaulttype;

template<class DataTypes>
class RegistrationRigidInternalData
{
public:
};

template<class DataTypes>
class RegistrationRigid : public virtual objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RegistrationRigid,DataTypes),sofa::core::objectmodel::BaseObject);

    typedef sofa::core::objectmodel::BaseObject Inherit;
	Kalmanfilter kalman;
    typedef defaulttype::ImageF DepthTypes;
	
	int npoints;

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecReal VecReal;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;
	typedef sofa::defaulttype::Vector4 Vector4;
	
	typedef core::topology::BaseMeshTopology::Edge Edge;

    typedef core::behavior::MechanicalState<DataTypes> MechanicalState;
    enum { N=DataTypes::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;

    typedef helper::fixed_array <unsigned int,3> tri;
    //typedef helper::kdTree<Coord> KDT;
    //typedef typename KDT::distanceSet distanceSet;
	
	double timef;
	
	vpHomogeneousMatrix cMo;
	Data< VecReal > translation;
	Data< VecReal > rotation;
	
	Data< double > errorfunction;
	
	typename core::behavior::MechanicalState<DataTypes> *mstate;

public:
    RegistrationRigid();
    virtual ~RegistrationRigid();
	
	static std::string templateName(const RegistrationRigid<DataTypes>* = NULL) { return DataTypes::Name();    }
    virtual std::string getTemplateName() const    { return templateName(this);    }

    // -- Base object interface
    void reinit();
    void init();
	void handleEvent(sofa::core::objectmodel::Event *event);
	void RegisterRigid();
	
	protected :
	
	VecCoord displ;
    Real min,max;
	
	typename sofa::core::objectmodel::RGBDDataProcessing<DataTypes>::SPtr rgbddataprocessing;
	typename sofa::core::objectmodel::MeshProcessing<DataTypes>::SPtr meshprocessing;
			
	Data<Vector4> barycenter;

    // source mesh data
    Data< helper::vector< tri > > sourceTriangles;
    Data< VecCoord > sourceNormals;
	Data< VecCoord > sourceSurfacePositions;
    Data< VecCoord > sourceSurfaceNormals;
	Data< VecCoord > sourceContourPositions;

	vector< bool > sourceVisible;  // flag ignored vertices
	vector< bool > sourceSurface;
	vector< bool > targetBackground;  // flag ignored vertices

	std::vector<int> indices;
	std::vector<int> indicesTarget;
	std::vector<int> indicesVisible;
	std::vector<int> sourceSurfaceMapping;
	
	VecCoord f_ ;       //WDataRefVecDeriv f(_f);
    VecCoord  x_ ;			//RDataRefVecCoord x(_x);
    VecCoord v_;	

    // target point cloud data
	Data< VecCoord > sourcePositions;
    Data< VecCoord > targetPositions;
	Data< VecCoord > targetGtPositions;
    Data< VecCoord > targetNormals;
    Data< helper::vector< tri > > targetTriangles;
    Data< VecCoord > rigidForces;
	
	vector< bool > targetBorder;
	Data< VecCoord > targetContourPositions;
    vector < double > sourceWeights;
	vector < double > combinedWeights;
	
	int ind;
	Data< VecCoord > sourceVisiblePositions;
	
	cv::Mat depth,depth_1, depthrend, depth00, depth01;	
	cv::Mat color, ir, ig, ib, gray;
	cv::Mat color_1,color_2, color_3, color_4, color_5, color_init;
	cv::Mat depthMap;
	cv::Mat silhouetteMap;

	// Number of iterations
	Data<int> niterations;
	Data<int> nimages;
	int npasses;
	
	// Paths
	Data<std::string> inputPath;
	Data<std::string> outputPath;
	Data<std::string> dataPath;

    cv::Mat foreground;
	Data<bool> useContour;
	Data<bool> useVisible;
	Data<bool> useRealData;
	Data<bool> useGroundTruth;
	Data<bool> generateSynthData;
	Data<bool> useSensor;
	Data<int> sensorType;
	Data<Vector4> cameraIntrinsicParameters;
	Eigen::Matrix3f rgbIntrinsicMatrix;
	int ntargetcontours;
	
	int iter_im;
	double timeOverall;
	double timeRigid;
	double timei,timeii;
	int timer;
	
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr source;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr source0;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr target;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetP;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetGt;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr source_registered;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr source_registered0;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetContour;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetPointCloud;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr sourceSurfacePointCloud;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr sourceSurfacePointCloud_registered;

	Eigen::Matrix4f transformation_matrix;

    std::vector<int> source2target_;
    std::vector<int> target2source_;
	std::vector<int> source2target_distances_;
    std::vector<int> target2source_distances_;
    pcl::CorrespondencesPtr correspondences_;
	std::vector<int> distances_;

    sofa::core::behavior::MechanicalState< DataTypes > *mstateRigid;
    Data< std::string > rigidState;

	void determineRigidTransformation ();
	void determineRigidTransformationVisible ();
	double determineErrorICP();
	};


/*#if defined(SOFA_EXTERN_TEMPLATE) && !defined(RegistrationRigid_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API RegistrationRigid<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API RegistrationRigid<defaulttype::Vec3fTypes>;
#endif
#endif*/


//

} //

} // namespace sofa

#endif
