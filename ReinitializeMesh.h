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

#ifndef SOFA_RGBDTRACKING_ReinitializeMesh_H
#define SOFA_RGBDTRACKING_ReinitializeMesh_H


#include <RGBDTracking/config.h>
#include <ImageTypes.h>
#include <sofa/core/core.h>
#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/component/interactionforcefield/SpringForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/component/component.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Quat.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/component/topology/TopologyData.h>
#include <visp/vpKltOpencv.h>
#include <sofa/helper/kdTree.inl>

#include <set>
#include "KalmanFilter.h"
#include <visp/vpDisplayX.h>
#include <algorithm>    // std::max

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
#include <sys/times.h>

#include <cuda_runtime.h>
#include <npp.h>
#include <nppi.h>

#define GL_GLEXT_PROTOTYPES 1
#define GL4_PROTOTYPES 1
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/glext.h>
#include <GL/glu.h>
#include <cuda_gl_interop.h>
#include <helper_cuda.h>
#include <helper_string.h>

#include <iostream>
#include <string>
#include <map>
//#include <XnCppWrapper.h>
#include <opencv2/opencv.hpp>
#include <boost/thread.hpp>

using namespace std;
using namespace cv;


namespace sofa
{

namespace component
{

namespace forcefield
{

using helper::vector;
using namespace sofa::defaulttype;
using namespace sofa::component::topology;

template<class DataTypes>
class ReinitializeMeshInternalData
{
public:
};

template<class DataTypes>
class ReinitializeMesh : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ReinitializeMesh,DataTypes),SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));
	
	int npoints;
    typedef core::behavior::ForceField<DataTypes> Inherit;

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;
	typedef sofa::defaulttype::Vector4 Vector4;
	
	enum { N=Vec3dTypes::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;
	typedef helper::fixed_array <unsigned int,3> tri;
	
	cv::Rect rectRtt;
    Real min,max;

    Data<Real> outlierThreshold;
    Data<Real> normalThreshold;
    Data<bool> rejectBorders;
    Data<bool> rejectOutsideBbox;
    defaulttype::BoundingBox targetBbox;
	int ind;
	//SoftKinetic softk;
	cv::Mat depth,depth_1, depthrend;
	cv::Mat color, ir, ig, ib, gray;
	//cv::Mat color_1,color_2, color_3, color_4, color_5, color_init;
	cv::Mat depthMap;

	// Number of iterations
	Data<int> nimages;
	Data<Real> sigmaWeight;
	int npasses;
	Data<int> sensorType;
	
	Data<Vector4> cameraIntrinsicParameters;
	Eigen::Matrix3f rgbIntrinsicMatrix;

	std::vector<int> indicesVisible;	
	
	int iter_im;
	
	    // source mesh data
	Data<Real> visibilityThreshold;
	Data< VecCoord > sourcePositions;
    Data< helper::vector< tri > > sourceTriangles;
    Data< VecCoord > sourceNormals;
	Data< VecCoord > sourceSurfacePositions;
    Data< VecCoord > sourceSurfaceNormals;
	Data< VecCoord > sourceSurfaceNormalsM;
	Data< VecCoord > sourceContourPositions;
	Data< VecCoord > sourceVisiblePositions;
	
	Data<int> borderThdSource;
	
	vector< bool > sourceBorder;
    vector< bool > sourceIgnored;  // flag ignored vertices
	vector< bool > sourceVisible;  // flag ignored vertices
	vector< bool > sourceSurface;
	vector < double > sourceWeights;
	
	cv::Mat rtd;

    ReinitializeMesh(core::behavior::MechanicalState<DataTypes> *mm = NULL);
    virtual ~ReinitializeMesh();

    void init(){};
	
	VecCoord getSourcePositions(){return sourcePositions.getValue();}
	VecCoord getSourceVisiblePositions(){return sourceVisiblePositions.getValue();}
	VecCoord getSourceContourPositions(){return sourceContourPositions.getValue();}

	void updateSourceSurface(); // built k-d tree and identify border vertices
	
	void extractSourceContour();
	void extractSourceVisibleContour();
	void getSourceVisible(double znear, double zfar);
	void updateSourceVisible();
	void updateSourceVisibleContour();
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(ReinitializeMesh_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API ReinitializeMesh<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API ReinitializeMesh<defaulttype::Vec3fTypes>;
#endif
#endif


} //

} //

} // namespace sofa

#endif
