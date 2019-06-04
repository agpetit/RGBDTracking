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

#ifndef SOFA_RGBDTRACKING_CLOSESTPOINTFORCEFIELD_H
#define SOFA_RGBDTRACKING_CLOSESTPOINTFORCEFIELD_H

#include <opencv/cv.h>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <visp/vpKltOpencv.h>
#include <visp/vpDisplayX.h>

#include "DataIO.h"

#include <image/ImageTypes.h>
#include <sofa/core/core.h>
#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/core/behavior/ForceField.h>
#include <SofaDeformable/SpringForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <SofaGeneralEngine/NormalsFromPoints.h>
//#include <sofa/helper/kdTree.inl>
#include <RGBDTracking/config.h>
#include "KalmanFilter.h"
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

#define GL_GLEXT_PROTOTYPES 1
#define GL4_PROTOTYPES 1
//#include <GL/glew.h>
//#include <GL/freeglut.h>
//#include <GL/glext.h>
//#include <GL/glu.h>

#include <string>
#include <boost/thread.hpp>
#include "ClosestPoint.h"
#include "RGBDDataProcessing.h"


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

template<class DataTypes>
class ClosestPointForceFieldInternalData
{
public:
};

template<class DataTypes>
class ClosestPointForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ClosestPointForceField,DataTypes),SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherit;
    typedef defaulttype::ImageF DepthTypes;

    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecReal VecReal;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;
    typedef sofa::defaulttype::Vector4 Vector4;
    typedef sofa::defaulttype::Vector2 Vec2;
	
    typedef core::behavior::MechanicalState<DataTypes> MechanicalState;
    enum { N=DataTypes::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;

    typedef typename interactionforcefield::LinearSpring<Real> Spring;
    typedef helper::fixed_array <unsigned int,3> tri;
	
    //typedef typename Coord::value_type real;
    typedef std::pair<int,Real> Col_Value;
    typedef vector< Col_Value > CompressedValue;
    typedef vector< CompressedValue > CompressedMatrix;
	
public:
    ClosestPointForceField(core::behavior::MechanicalState<DataTypes> *mm = NULL);
    virtual ~ClosestPointForceField();

    core::behavior::MechanicalState<DataTypes>* getObject() { return this->mstate; }
	
    static std::string templateName(const ClosestPointForceField<DataTypes>* = NULL) { return DataTypes::Name();    }
    virtual std::string getTemplateName() const    { return templateName(this);    }

    const sofa::helper::vector< Spring >& getSprings() const {return springs.getValue();}

    // -- ForceField interface
    void reinit();
    void init();
    void addForce(const core::MechanicalParams* mparams,DataVecDeriv& f , const DataVecCoord& x , const DataVecDeriv& v);
    void addDForce(const core::MechanicalParams* mparams ,DataVecDeriv&   df , const DataVecDeriv&   dx);
    double getPotentialEnergy(const core::MechanicalParams* ,const DataVecCoord&) const { return m_potentialEnergy; }
    //void addKToMatrix( const core::MechanicalParams* mparams,const sofa::core::behavior::MultiMatrixAccessor* matrix);
    virtual void addKToMatrix(sofa::defaulttype::BaseMatrix *m, SReal kFactor, unsigned int &offset);

    Real getStiffness() const{ return ks.getValue(); }
    Real getDamping() const{ return kd.getValue(); }
    void setStiffness(Real _ks){ ks.setValue(_ks); }
    void setDamping(Real _kd){ kd.setValue(_kd); }
    Real getArrowSize() const{return showArrowSize.getValue();}
    void setArrowSize(float s){showArrowSize.setValue(s);}
    int getDrawMode() const{return drawMode.getValue();}
    void setDrawMode(int m){drawMode.setValue(m);}

    void draw(const core::visual::VisualParams* vparams);

    // -- Modifiers

    void clearSprings(int reserve=0)
    {
        sofa::helper::vector<Spring>& springs = *this->springs.beginEdit();
        springs.clear();
        if (reserve) springs.reserve(reserve);
        this->springs.endEdit();
    }

    void removeSpring(unsigned int idSpring)
    {
        if (idSpring >= (this->springs.getValue()).size())
            return;

        sofa::helper::vector<Spring>& springs = *this->springs.beginEdit();
        springs.erase(springs.begin() +idSpring );
        this->springs.endEdit();
    }

    void addSpring(int m1, SReal ks, SReal kd )
    {
        springs.beginEdit()->push_back(Spring(m1,-1,ks,kd,0));
        springs.endEdit();
    }

    void addSpring(const Spring & spring)
    {
        springs.beginEdit()->push_back(spring);
        springs.endEdit();
    }

    protected :

    int npoints;
    vector<Mat>  dfdx;
    vector<Mat>  dfdx1;
    VecCoord closestPos;
    vector<unsigned int>  cnt;
    double m_potentialEnergy;	
    Real min,max;

    /// Accumulate the spring force and compute and store its stiffness
    virtual void addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
    virtual void addStoredSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
    virtual void addSpringForceWeight(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, int ivis, const Spring& spring);
    /// Apply the stiffness, i.e. accumulate df given dx
    virtual void addSpringDForce(VecDeriv& df,const  VecDeriv& dx, int i, const Spring& spring, double kFactor, double bFactor);
    virtual void addSpringForceKLTA(double& potentialEnergy, VecDeriv& f, const  VecCoord& p,const VecDeriv& v, Coord& KLTtarget, int i, const Spring& spring, double coef);


    Data<Real> ks;
    Data<Real> kd;
    Data<Real> blendingFactor;
    Data<bool> projectToPlane;
    Data<sofa::helper::vector<Spring> > springs;

    typename sofa::core::objectmodel::ClosestPoint<DataTypes> *closestpoint;

    Data<Vector4> cameraIntrinsicParameters;
    Eigen::Matrix3f rgbIntrinsicMatrix;
	
    // source mesh data
    Data< helper::vector< tri > > sourceTriangles;
    Data< VecCoord > sourceNormals;
    Data< VecCoord > sourceSurfacePositions;
    Data< VecCoord > sourceSurfaceNormals;
    Data< VecCoord > sourceContourPositions;
    Data< helper::vector< Vec2 > > sourceContourNormals;

    helper::vector< Real > sourcew;

    Data< helper::vector< bool > > sourceVisible;  // flag visiblevertices
    Data< helper::vector< bool > > sourceBorder;

    vector< bool > sourceSurface;
    vector< bool > targetBackground;  // flag ignored vertices   

    std::vector<int> indices;
    std::vector<int> indicesTarget;
    Data< helper::vector<int> > indicesVisible;
	
    VecCoord f_;//WDataRefVecDeriv f(_f);
    VecCoord x_;//RDataRefVecCoord x(_x);
    VecCoord v_;	

    // target point cloud data
    Data< VecCoord > sourcePositions;
    Data< VecCoord > targetPositions;
    Data< VecCoord > targetNormals;
    Data< helper::vector< bool > > targetBorder;
    Data< helper::vector < double > > targetWeights;
    Data< VecCoord > targetContourPositions;
    Data< helper::vector < double > > sourceWeights;
    vector < double > combinedWeights;
    Data< helper::vector< double > > curvatures;

    int ind;
    Data< VecCoord > sourceVisiblePositions;

    cv::Mat depth,depth_1, depth00;
    cv::Mat color, ir, ig, ib, gray;
    cv::Mat color_1;
    cv::Mat depthMap;
    cv::Mat silhouetteMap;
    cv::Mat distimage, dotimage;
    vpKltOpencv tracker,tracker1;
    vpImage<unsigned char> vpI ;
    vpDisplayX display;
    Data<int> windowKLT;
    Data< VecCoord > targetKLTPositions;

    void mapKLTPointsTriangles ( helper::vector< tri > &triangles);
    void KLTPointsTo3D();

    sofa::helper::vector<Vector3> mappingkltcoef;
    sofa::helper::vector<int> mappingkltind;

    Data<Real> outlierThreshold;
    Data<bool> rejectBorders;

    Data<float> showArrowSize;
    Data<int> drawMode; //Draw Mode: 0=Line - 1=Cylinder - 2=Arrow
    Data<bool> drawColorMap;
    Data<bool> theCloserTheStiffer;

    // Number of iterations
    Data<int> niterations;
    Data<int> startimage;
    Data<int> startimageklt;
    Data<int> nimages;
    int npasses;
    Data<bool> useContour;
    Data<bool> useContourWeight;
    Data<bool> useDistContourNormal;
    Data<bool> useVisible;
    Data<bool> useRealData;
    Data<bool> drawContour;

    Data<bool> useKLTPoints;

    Data<std::string> dataPath;
	
    std::vector<bool>* visible;
    int iter_im;

    std::ofstream filerror;
    std::ofstream fileindicesVisible, fileVisible, filePCDMatch, filePCD;

    void resetSprings();
    void addForceMesh(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v );

    double computeError(Vector3 sourcePoint, Vector3 targetPoint);	
		
};


/*#if defined(SOFA_EXTERN_TEMPLATE) && !defined(ClosestPointForceField_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API ClosestPointForceField<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API ClosestPointForceField<defaulttype::Vec3fTypes>;
#endif
#endif*/


} //

} //

} // namespace sofa

#endif
