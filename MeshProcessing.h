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

#ifndef SOFA_RGBDTRACKING_MESHPROCESSING_H
#define SOFA_RGBDTRACKING_MESHPROCESSING_H

#include <opencv/cv.h>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/helper/gl/FrameBufferObject.h>

#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <SofaBaseTopology/TopologyData.h>
#include <RGBDTracking/config.h>
#include <algorithm>    // std::max

#include <stdlib.h>
#include <vector>
#include <math.h>
#include <sys/times.h>

#define GL_GLEXT_PROTOTYPES 1
#define GL4_PROTOTYPES 1
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/glext.h>
#include <GL/glu.h>

#include <boost/thread.hpp>

#include <image/ImageTypes.h>
#include "RenderingManager.h"

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
class MeshProcessingInternalData
{
public:
};

template<class DataTypes>
class MeshProcessing : public sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(MeshProcessing,DataTypes),sofa::core::objectmodel::BaseObject);
	
	int npoints;
    typedef sofa::core::objectmodel::BaseObject Inherit;

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
	
    typename core::behavior::MechanicalState<DataTypes> *mstate;
    typename sofa::component::visualmodel::RenderingManager::SPtr renderingmanager;

    cv::Rect rectRtt;
    Data<Vector4> BBox;
    Real min,max;
    int hght;
    int wdth;

    int ind;
    //SoftKinetic softk;
    cv::Mat depth,depth_1, depthrend;
    cv::Mat color, ir, ig, ib, gray;
    //cv::Mat color_1,color_2, color_3, color_4, color_5, color_init;
    cv::Mat depthMap;

    // Number of iterations
    Data<int> niterations;
    Data<int> nimages;
    Data<Real> sigmaWeight;
    int npasses;
	
    int iter_im;
	
    // source mesh data
		
    Data<bool> useContour;
    Data<bool> useVisible;
    Data<bool> useRealData;
    Data<Vector4> cameraIntrinsicParameters;
    Eigen::Matrix3f rgbIntrinsicMatrix;
    Data<Real> visibilityThreshold;
    Data< helper::vector<int> > indicesVisible;

    Data< VecCoord > sourcePositions;
    Data< helper::vector< tri > > sourceTriangles;
    Data< VecCoord > sourceNormals;
    Data< VecCoord > sourceSurfacePositions;
    Data< VecCoord > sourceSurfaceNormals;
    Data< VecCoord > sourceSurfaceNormalsM;
    Data< VecCoord > sourceContourPositions;
    Data< VecCoord > sourceVisiblePositions;
	
    Data<int> borderThdSource;

    Data< helper::vector< bool > > sourceBorder;
    vector< bool > sourceIgnored;  // flag ignored vertices
    Data<helper::vector< bool > > sourceVisible;  // flag ignored vertices
    vector< bool > sourceSurface;
    vector < double > sourceWeights;

    double timeMeshProcessing;

    MeshProcessing();
    virtual ~MeshProcessing();

    void init();
    void handleEvent(sofa::core::objectmodel::Event *event);

    VecCoord getSourcePositions(){return sourcePositions.getValue();}
    VecCoord getSourceVisiblePositions(){return sourceVisiblePositions.getValue();}
    VecCoord getSourceContourPositions(){return sourceContourPositions.getValue();}
    void setViewPoint();

    void updateSourceSurface(); // built k-d tree and identify border vertices
	
    void extractSourceContour();
    void extractSourceVisibleContour();
    void getSourceVisible(double znear, double zfar);
    void updateSourceVisible();
    void updateSourceVisibleContour();
    void draw(const core::visual::VisualParams* vparams);
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(MeshProcessing_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API MeshProcessing<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API MeshProcessing<defaulttype::Vec3fTypes>;
#endif
#endif


} //

} //

} // namespace sofa

#endif
