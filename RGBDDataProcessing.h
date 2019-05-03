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

#ifndef SOFA_RGBDTRACKING_RGBDDATAPROCESSING_H
#define SOFA_RGBDTRACKING_RGBDDATAPROCESSING_H

#include <RGBDTracking/config.h>
#include <boost/thread.hpp>
#include "DataIO.h"

#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <SofaBaseTopology/TopologyData.h>
#include <sofa/gui/BaseGUI.h>
#include <sofa/gui/BaseViewer.h>
#include <sofa/gui/GUIManager.h>

#define GL_GLEXT_PROTOTYPES 1
#define GL4_PROTOTYPES 1
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/glext.h>
#include <GL/glu.h>

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

#include <visp/vpIoTools.h>
#include <visp/vpImageIo.h>
#include <visp/vpParseArgv.h>
#include <visp/vpMatrix.h>
#include <visp/vpKltOpencv.h>

#include "segmentation.h"

//#include "ImageConverter.h"

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
class RGBDDataProcessingInternalData
{
public:
};

template<class DataTypes>
class RGBDDataProcessing : public sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RGBDDataProcessing,DataTypes),sofa::core::objectmodel::BaseObject);
	
    typedef sofa::core::objectmodel::BaseObject Inherit;

    int npoints;

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;
    typedef sofa::defaulttype::Vector4 Vector4;
    typedef sofa::defaulttype::Vector3 Vec3;
	
    typedef defaulttype::ImageF DepthTypes;
	
    enum { N=Vec3dTypes::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;
	
    double timef;
    cv::Rect rectRtt;

    VecCoord IntensityCurrent, IntensityPrev;
    Real min,max;

    Data<Real> outlierThreshold;
    Data<Real> normalThreshold;
    Data<bool> rejectBorders;
    Data<bool> rejectOutsideBbox;
    defaulttype::BoundingBox targetBbox;
			
    typedef helper::fixed_array <unsigned int,3> tri;

    // target point cloud data
    Data< VecCoord > targetPositions;
    Data< VecCoord > targetContourPositions;
    Data< VecCoord > targetGtPositions;
    Data< VecCoord > targetNormals;
    Data< VecCoord > targetKLTPositions;

    Data< helper::vector< tri > > targetTriangles;
    Data< helper::vector< bool > > targetBorder;
    Data< helper::vector< double > > targetWeights;

    int ind;

    cv::Mat depth,depth_1, depth00;
    cv::Mat color, ir, ig, ib, gray;
    cv::Mat color_1;
    cv::Mat depthMap;
    cv::Mat silhouetteMap;
    cv::Mat distimage, dotimage;

    Data<int> imagewidth;
    Data<int> imageheight;

    // Number of iterations
    Data<int> niterations;
    int npasses;
    Data<int> nimages;
    Data<Real> sigmaWeight;

    Data<int> samplePCD;
    Data<int> offsetX, offsetY;
    Data<int> borderThdPCD;
    Data<int> windowKLT;
    Data<bool> useDistContourNormal;

    // Paths
    Data<std::string> inputPath;
    Data<std::string> outputPath;
    Data<std::string> dataPath;
    Data<std::string> ipad;

    segmentation seg;
    Data<int> segNghb;
    Data<int> segImpl;
    Data<int> segMsk;
	
    cv::Mat foreground, foregroundS;
    bool pcl;
    bool disp;

    Data<bool> useContour;
    Data<bool> useRealData;
    Data<bool> useGroundTruth;
    Data<bool> useKLTPoints;
    Data<int> sensorType;
    Data<bool> useSensor;
    Data<bool> cameraChanged;

    Data< int > scaleImages;
    Data< bool > displayImages;
    Data< int > displayDownScale;
    Data< bool > saveImages;
    Data< bool > displaySegmentation;
    Data< int > scaleSegmentation;
    Data<bool> drawPointCloud;
    Data<bool> displayBackgroundImage;

    Data<bool> useCurvature;
    Data< helper::vector< double > > curvatures;
    Data<bool> useSIFT3D;

    Data<Vector4> cameraIntrinsicParameters;
    Eigen::Matrix3f rgbIntrinsicMatrix;

    Data<Vec3> cameraPosition;
    Data<Quat> cameraOrientation;

    Data<bool> stopatinit;
    Data<bool> safeModeSeg;
    Data<double> segTolerance;

	
    int ntargetcontours;
    int iter_im;
    int sizeinit;

    double timeSegmentation;
    double timePCD;

    bool initsegmentation;
	
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr target;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetP;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetGt;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetContour;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetPointCloud;
		
    RGBDDataProcessing();
    virtual ~RGBDDataProcessing();

    void init();
	void handleEvent(sofa::core::objectmodel::Event *event);
	void setRGBDData(cv::Mat &_color, cv::Mat &_depth)
	{
		color = _color;
		depth = _depth;
	}
	
	
    void computeTargetNormals();
    Eigen::Matrix<float,1,3>& computeJacobian(int i);
    Eigen::Matrix<float,3,3>& computeJacobianColor(int i);

    void initTarget();  // built k-d tree and identify border vertices
    void initTargetContour();  // built k-d tree and identify border vertices
    vector < double > combinedWeights;
    VecCoord getTargetPositions(){return targetPositions.getValue();}
    VecCoord getTargetContourPositions(){return targetContourPositions.getValue();}

    void detectBorder(vector<bool> &border,const helper::vector< tri > &triangles);
    void computeCenter(vpImage<unsigned char> &Itemp, vpImagePoint &cog,double &angle, int &surface);
    void extractTargetPCD();
    void extractTargetPCDContour();
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr PCDFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr PCDContourFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage, cv::Mat& distImage, cv::Mat& dotImage);
    void setCameraPose();

    void initSegmentation();
    void segment();
    void segmentSynth();
    void ContourFromRGBSynth(cv::Mat& rgbImage, cv::Mat& distImage, cv::Mat& dotImage);
    void draw(const core::visual::VisualParams* vparams) ;
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(RGBDDataProcessing_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API RGBDDataProcessing<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API RGBDDataProcessing<defaulttype::Vec3fTypes>;
#endif
#endif


} //

} //

} // namespace sofa

#endif
