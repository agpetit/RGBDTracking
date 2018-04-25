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

#ifndef SOFA_RGBDTRACKING_IMAGECONVERTER_H
#define SOFA_RGBDTRACKING_IMAGECONVERTER_H

#include <RGBDTracking/config.h>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

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
#include <visp/vpKltOpencv.h>
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

#include <cuda_runtime.h>
#include <npp.h>
#include <nppi.h>

#define GL_GLEXT_PROTOTYPES 1
#define GL4_PROTOTYPES 1
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/glext.h>
#include <GL/glu.h>

#include <string>
#include <boost/thread.hpp>


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
using cimg_library::CImg;

template<class DataTypes, class _ImageTypes>
class ImageConverterInternalData
{
public:
};

template<class DataTypes, class _ImageTypes>
class ImageConverter : public sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(ImageConverter,DataTypes,_ImageTypes), sofa::core::objectmodel::BaseObject);
    typedef sofa::core::objectmodel::BaseObject Inherit;
    typedef _ImageTypes DepthTypes;
    typedef typename DepthTypes::T dT;
    typedef helper::WriteAccessor<Data< DepthTypes > > waDepth;
    typedef helper::ReadAccessor<Data< DepthTypes > > raDepth;
    Data< DepthTypes > depthImage;
	
    typedef defaulttype::ImageUC ImageTypes;
	typedef typename ImageTypes::T T;
	typedef helper::WriteAccessor<Data< ImageTypes > > waImage;
    typedef helper::ReadAccessor<Data< ImageTypes > > raImage;
    Data< ImageTypes > image;
	
	typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecReal VecReal;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;

	int npoints;

    typedef core::behavior::MechanicalState<DataTypes> MechanicalState;
	typename core::behavior::MechanicalState<DataTypes> *mstate;

	double timef;
	cv::Rect rectRtt;

public:
    ImageConverter();
    virtual ~ImageConverter();

    core::behavior::MechanicalState<DataTypes>* getObject() { return mstate; }
	
	static std::string templateName(const ImageConverter<DataTypes,DepthTypes >* = NULL) { return DataTypes::Name()+ std::string(",")+DepthTypes::Name();    }
    virtual std::string getTemplateName() const    { return templateName(this);    }

    // -- ForceField interface
    void init();
	void handleEvent(sofa::core::objectmodel::Event *event);

    void draw(const core::visual::VisualParams* vparams);

	void getImages();
		
	cv::Mat depth,depth_1, depthrend, depth00, depth01;	
	cv::Mat color, ir, ig, ib, gray;
	cv::Mat color_1,color_2, color_3, color_4, color_5, color_init;
	cv::Mat depthMap;
	cv::Mat silhouetteMap;
	
	cv::Mat getColor(){return color;}


	// Number of iterations
	Data<int> niterations;
	Data<int> nimages;
	int npasses;

	Data<bool> useRealData;
	Data<bool> useSensor;
	Data<int> sensorType;

	Data<Vector4> cameraIntrinsicParameters;
	Eigen::Matrix3f rgbIntrinsicMatrix;
	
	int iter_im;
		
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(ImageConverter_CPP)
    #ifndef SOFA_FLOAT
     extern template class SOFA_RGBDTRACKING_API ImageConverter<Vec3dTypes,ImageUC>;
     extern template class SOFA_RGBDTRACKING_API ImageConverter<Vec3dTypes,ImageUS>;
     extern template class SOFA_RGBDTRACKING_API ImageConverter<Vec3dTypes,ImageF>;
    #endif
    #ifndef SOFA_DOUBLE
    extern template class SOFA_RGBDTRACKING_API ImageConverter<Vec3fTypes,ImageUC>;
    extern template class SOFA_RGBDTRACKING_API ImageConverter<Vec3fTypes,ImageUS>;
    extern template class SOFA_RGBDTRACKING_API ImageConverter<Vec3fTypes,ImageF>;
    #endif

#endif


} //

} //

} // namespace sofa

#endif
