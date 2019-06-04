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

#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/core/objectmodel/Link.h>
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
//#include <GL/glext.h>
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

#include <pcl/search/impl/search.hpp>

#include "DataIO.h"
#include "ImageConverter.h"
#include <RGBDTracking/src/Segmentor/segmentation.h>

using namespace std;
using namespace cv;


namespace sofa {

namespace rgbdtracking {

using helper::vector;
using namespace sofa::defaulttype;
using namespace sofa::component::topology;

class MouseEventHandler {
public :

    static bool destroy;
    static cv::Rect box;
    static bool drawing_box ;

    static void init () {
        destroy=false ;
        drawing_box = false ;
    }

    // used in mouse callback function
    static void draw_box(cv::Mat _img, cv::Rect rect) {
        cv::rectangle(
            _img,
            cvPoint(rect.x, rect.y),
            cvPoint(rect.x+rect.width, rect.y+rect.height),
            cvScalar(0,0,255),
            2);
        //cv::Rect rect2=cv::Rect(box.x,box.y,box.width,box.height);
        //cvSetImageROI(image, rect2);   //here I wanted to set the drawn rect as ROI
    }

    // Implement mouse callback
    static void my_mouse_callback( int event, int x, int y, int /*flags*/, void* param ) {
      cv::Mat* frame = (cv::Mat*) param;

      switch( event ) {
          case CV_EVENT_MOUSEMOVE:
              if( drawing_box ) {
                  box.width = x-box.x;
                  box.height = y-box.y;
              }
          break;

          case CV_EVENT_LBUTTONDOWN:
              drawing_box = true;
              box = cvRect( x, y, 0, 0 );
          break;

          case CV_EVENT_LBUTTONUP:
              drawing_box = false;
              if( box.width < 0 ) {
                  box.x += box.width;
                  box.width *= -1;
              }

              if( box.height < 0 ) {
                  box.y += box.height;
                  box.height *= -1;
              }

              MouseEventHandler::draw_box(*frame, box);
          break;

          case CV_EVENT_RBUTTONUP:
              destroy=true;
          break;

          default:
          break;
       }

    }

} ;

template<class DataTypes>
class RGBDDataProcessing : public sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RGBDDataProcessing,DataTypes),sofa::core::objectmodel::BaseObject);
	
    typedef sofa::core::objectmodel::BaseObject Inherit;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;
    typedef sofa::defaulttype::Vector4 Vector4;
    typedef sofa::defaulttype::Vector3 Vec3;
    typedef helper::fixed_array <unsigned int,3> tri;
//    typedef defaulttype::ImageF DepthTypes;
	
    typedef defaulttype::Mat<
        Vec3dTypes::spatial_dimensions,
        Vec3dTypes::spatial_dimensions,
        Real
    > Mat;

    core::objectmodel::SingleLink<
        RGBDDataProcessing<DataTypes>,
        ImageConverter<DataTypes, ImageF>,
        BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_imconv ;

    core::objectmodel::SingleLink<
        RGBDDataProcessing<DataTypes>,
        DataIO<DataTypes>,
        BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_dataio ;

    Data<Vector4> cameraIntrinsicParameters;

    Data<bool> useContour;
    Data<bool> useSensor;
    Data<int> sensorType;

    Data<int> niterations;

    Data<int> samplePCD;
    Data<Real> sigmaWeight;
    Data<int> borderThdPCD;

    //segmentation param
    segmentation seg;
    Data<int> segNghb;
    Data<int> segImpl;
    Data<int> segMsk;

    Data< int > scaleSegmentation;

    //display param
    Data< int > scaleImages;
    Data< bool > displayImages;
    Data< int > displayDownScale;
    Data< bool > displaySegmentation;
    Data<bool> drawPointCloud;
    Data<bool> displayBackgroundImage;

    Data< bool > saveImages;
    Data<bool> useCurvature;

    //outputs
    Data< VecCoord > targetPositions;
    Data< VecCoord > targetNormals;
    Data< VecCoord > targetContourPositions;
    Data< VecCoord > targetGtPositions;
    Data< helper::vector< bool > > targetBorder;
    Data< helper::vector< double > > targetWeights;
    Data< helper::vector< double > > curvatures;

    Data<Vec3> cameraPosition;
    Data<Quat> cameraOrientation;
    Data<bool> cameraChanged;

    Data<bool> stopatinit;
    Data<bool> safeModeSeg;
    Data<double> segTolerance;

    Eigen::Matrix3f rgbIntrinsicMatrix;

    cv::Mat foreground, foregroundS;
    cv::Mat depth, color ;
    cv::Mat distimage, dotimage;

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr target ;

    int ntargetcontours;
    int sizeinit;

    RGBDDataProcessing();
    virtual ~RGBDDataProcessing();

    void init();
    void handleEvent(sofa::core::objectmodel::Event *event);
    void setRGBDData(cv::Mat &_color, cv::Mat &_depth) {
        color = _color;
        depth = _depth;
    }
	
    vector < double > combinedWeights;
    VecCoord getTargetPositions(){return targetPositions.getValue();}
    VecCoord getTargetContourPositions(){return targetContourPositions.getValue();}

    void computeCenter(vpImage<unsigned char> &Itemp, vpImagePoint &cog,double &angle, int &surface);
    void extractTargetPCD();
    void extractTargetPCDContour();
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr PCDFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr PCDContourFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage, cv::Mat& distImage, cv::Mat& dotImage);
    void setCameraPose();

    void initSegmentation();
    void segment();
    void draw(const core::visual::VisualParams* vparams) ;

private :
    bool loadFromImageConverter () ;

    bool loadFromDataIO () ;
    void saveToDataIO () ;
    void saveSegmentationToDataIO () ;

    void displayDownScaledImage () ;
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(RGBDDataProcessing_CPP)
#ifndef SOFA_FLOAT
    extern template class SOFA_RGBDTRACKING_API RGBDDataProcessing<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
    extern template class SOFA_RGBDTRACKING_API RGBDDataProcessing<defaulttype::Vec3fTypes>;
#endif
#endif


} // rgbdtracking

} // namespace sofa

#endif
