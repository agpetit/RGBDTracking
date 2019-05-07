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

#ifndef SOFA_RGBDTRACKING_DATAIO_H
#define SOFA_RGBDTRACKING_DATAIO_H

#include <pcl/common/common_headers.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <RGBDTracking/config.h>
#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>

#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/defaulttype/VecTypes.h>

#include <opencv/cv.h>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <set>


#ifdef WIN32
#include <process.h>
#else
#include <pthread.h>
#endif

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>

#include <sys/times.h>

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

template<class DataTypes>
class DataIOInternalData
{
public:
};

template<class DataTypes>
class DataIO : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(DataIO,DataTypes), sofa::core::objectmodel::BaseObject);
		
    typedef sofa::core::objectmodel::BaseObject Inherit;
	
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecReal VecReal;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;

    // Number of iterations
    Data<int> nimages;
    Data<int> startimage;
    int npasses;
	
    // Paths
    Data<std::string> inputPath;
    Data<std::string> outputPath;
    Data<std::string> dataPath;
    Data<std::string> ipad;
	
    bool pcl;
    bool disp;
    Data<bool> useRealData;
    Data<bool> useGroundTruth;
    Data<bool> useKLTPoints;
    Data<int> sensorType;
    Data<bool> useSensor;
    Data<int> niterations;
    Data<bool> newImages;
	
    int ntargetcontours;
	
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr target;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetP;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetGt;

    int ind;
    cv::Mat depth,depth_1, depthrend, depth00, depth01;
    cv::Mat color, ir, ig, ib, gray;
    cv::Mat color_1,color_2, color_3, color_4, color_5, color_init;
    cv::Mat depthMap;

    int hght,wdth;
	
    Data< VecCoord > targetPositions;
    Data< VecCoord > targetGtPositions;
 
    std::vector<cv::Mat*> listimg;
    std::vector<cv::Mat*> listimgklt;
    std::vector<cv::Mat*> listimgseg;
    std::vector<cv::Mat*> listdepth;
    std::vector<cv::Mat*> listrtt;
    std::vector<cv::Mat*> listrttstress;
    std::vector<cv::Mat*> listrttstressplast;
    std::vector<std::vector<Vec3d>*> listpcd;
    std::vector<std::vector<bool>*> listvisible;
    std::vector<std::vector<double>*> listvm;
    std::vector<std::vector<double>*> listps;
    std::vector<std::vector<double>*> listes;
    std::vector<std::vector<double>*> listts;
    std::vector<std::vector<double>*> listesnode;
    std::vector<std::vector<double>*> listpsnode;
    std::vector<std::vector<double>*> listtsnode;
	
    std::vector<double> vonmisesstressGt;
    std::vector<double> elasticstrainsGt;
    std::vector<double> plasticstrainsGt;
    std::vector<double> totalstrainsGt;
    std::vector<double> elasticstrainsnodeGt;
    std::vector<double> plasticstrainsnodeGt;
    std::vector<double> totalstrainsnodeGt;

    cv::Mat* imgl;
    cv::Mat* imgklt;
    cv::Mat* imglsg;
    cv::Mat* depthl;
    cv::Mat* rtt;
	
    int iter_im;
    cv::Mat rtd;
	
    DataIO();
    virtual ~DataIO();

    void init();
    void handleEvent(sofa::core::objectmodel::Event *event);
	
    void readData();
    void readImages();
	
    void writeImages();
    void writeImagesSynth();
    void writeData();
		
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(DataIO_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API DataIO<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API DataIO<defaulttype::Vec3fTypes>;
#endif
#endif

}

} //

} // namespace sofa

#endif
