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

#ifndef SOFA_RGBDTRACKING_RENDERTEXTUREAR_H
#define SOFA_RGBDTRACKING_RENDERTEXTUREAR_H

#include <RGBDTracking/config.h>
#include <image/ImageTypes.h>
#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>

#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Quat.h>
#include <sofa/defaulttype/Vec.h>
#include <SofaBaseTopology/TopologyData.h>
#include <visp/vpKltOpencv.h>

#include <set>
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

#include <opencv/cv.h>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <boost/thread.hpp>

#include <iostream>
#include <string>
#include <map>
//#include <XnCppWrapper.h>

#include "luaconfig.h"
using namespace std;
using namespace cv;


namespace sofa
{

namespace core
{

namespace objectmodel
{

using helper::vector;
using namespace sofa::defaulttype;
using namespace sofa::component::topology;

template<class DataTypes>
class RenderTextureARInternalData
{
public:
};

template<class DataTypes>
class RenderTextureAR : public core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RenderTextureAR,DataTypes), sofa::core::objectmodel::BaseObject);
	
	typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;
	typedef sofa::defaulttype::Vector4 Vector4;
	
	int npoints;
    typedef sofa::core::objectmodel::BaseObject Inherit;

	cv::Rect rectRtt;
    Real min,max;
	int ind;
	Data<Vector4> cameraIntrinsicParameters;
	Eigen::Matrix3f rgbIntrinsicMatrix;

    RenderTextureAR();
    virtual ~RenderTextureAR();

    void init(){};
	void renderToTexture(cv::Mat &_rtt);
    void renderToTextureD(cv::Mat &_rtt,cv::Mat &color_1);
	void renderToTextureDepth(cv::Mat &_rtt, cv::Mat &_rttdepth);

};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(RenderTextureAR_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API RenderTextureAR<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API RenderTextureAR<defaulttype::Vec3fTypes>;
#endif
#endif


} //

} //

} // namespace sofa

#endif
