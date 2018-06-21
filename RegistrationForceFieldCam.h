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

#ifndef SOFA_RGBDTRACKING_REGISTRATIONFORCEFIELDCAM_H
#define SOFA_RGBDTRACKING_REGISTRATIONFORCEFIELDCAM_H

#include <opencv/cv.h>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

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
#include <visp/vpKltOpencv.h>
#include <SofaGeneralEngine/NormalsFromPoints.h>
//#include <sofa/helper/kdTree.inl>
#include <RGBDTracking/config.h>
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

#define GL_GLEXT_PROTOTYPES 1
#define GL4_PROTOTYPES 1
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <GL/glext.h>
#include <GL/glu.h>


#include <visp/vpIoTools.h>
#include <visp/vpImageIo.h>
#include <visp/vpParseArgv.h>
#include <visp/vpMatrix.h>

#include <string>
#include <boost/thread.hpp>
#include "ClosestPoint.h"
#include "ccd.h"
#include "p_helper.h"
#include "RGBDDataProcessing.h"

#include "MeshProcessing.h"
#include "RenderTextureAR.h"
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

namespace component
{

namespace forcefield
{

using helper::vector;
using namespace sofa::defaulttype;

template<class DataTypes>
class RegistrationForceFieldCamInternalData
{
public:
};

template<class DataTypes>
class RegistrationForceFieldCam : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RegistrationForceFieldCam,DataTypes),SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

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
	
    typedef core::topology::BaseMeshTopology::Edge Edge;

    typedef core::behavior::MechanicalState<DataTypes> MechanicalState;
    enum { N=DataTypes::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;

    typedef typename interactionforcefield::LinearSpring<Real> Spring;
    typedef helper::fixed_array <unsigned int,3> tri;
	
    //typedef typename Coord::value_type real;
    typedef std::pair<int,Real> Col_Value;
    typedef vector< Col_Value > CompressedValue;
    typedef vector< CompressedValue > CompressedMatrix;
	
    Data< double > errorfunction;

public:
    RegistrationForceFieldCam(core::behavior::MechanicalState<DataTypes> *mm = NULL);
    virtual ~RegistrationForceFieldCam();

    core::behavior::MechanicalState<DataTypes>* getObject() { return this->mstate; }
	
    static std::string templateName(const RegistrationForceFieldCam<DataTypes>* = NULL) { return DataTypes::Name();    }
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

    double timef;
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
    virtual void addSpringForceWeight(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
    /// Apply the stiffness, i.e. accumulate df given dx
    virtual void addSpringDForce(VecDeriv& df,const  VecDeriv& dx, int i, const Spring& spring, double kFactor, double bFactor);

    Data<Real> ks;
    Data<Real> kd;
    Data<Real> blendingFactor;
    Data<bool> projectToPlane;
    Data<sofa::helper::vector<Spring> > springs;
	
    typename sofa::core::objectmodel::RGBDDataProcessing<DataTypes>::SPtr rgbddataprocessing;
    typename sofa::core::objectmodel::MeshProcessing<DataTypes>::SPtr meshprocessing;
    typename sofa::core::objectmodel::RenderTextureAR<DataTypes>::SPtr rendertexturear;
    typename sofa::core::objectmodel::ClosestPoint<DataTypes>::SPtr closestpoint;
    typename sofa::core::objectmodel::DataIO<DataTypes>::SPtr dataio;
    //typename ImageConverter<DataTypes,DepthTypes>::SPtr imconv;
		
    VecCoord tpos;
	
    // source mesh data
    Data< helper::vector< tri > > sourceTriangles;
    Data< VecCoord > sourceNormals;
    Data< VecCoord > sourceSurfacePositions;
    Data< VecCoord > sourceSurfaceNormals;
    Data< VecCoord > sourceContourPositions;

    vector< bool > sourceVisible;  // flag visiblevertices
    vector< bool > sourceSurface;
    vector< bool > targetBackground;  // flag ignored vertices

    std::vector<int> indices;
    std::vector<int> indicesTarget;
    std::vector<int> indicesVisible;
	
    VecCoord f_ ;       //WDataRefVecDeriv f(_f);
    VecCoord  x_ ;			//RDataRefVecCoord x(_x);
    VecCoord v_;	

    // target point cloud data
    Data< VecCoord > sourcePositions;
    Data< VecCoord > targetPositions;
    Data< VecCoord > targetNormals;
    vector< bool > targetBorder;
    vector < double > targetWeights;
    Data< VecCoord > targetContourPositions;
    vector < double > sourceWeights;
    vector < double > combinedWeights;
	
    int ind;
    Data< VecCoord > sourceVisiblePositions;

    Data<float> showArrowSize;
    Data<int> drawMode; //Draw Mode: 0=Line - 1=Cylinder - 2=Arrow
    Data<bool> drawColorMap;
    Data<bool> theCloserTheStiffer;


    // Number of iterations
    Data<int> niterations;
    Data<int> nimages;
    int npasses;
    Data<bool> useContour;
    Data<bool> useVisible;
    Data<bool> useRealData;
    Data<bool> useSensor;
    Data<bool> drawSource;
    Data<bool> drawTarget;
    Data<bool> drawContour;
	
    int ntargetcontours;

    std::vector<bool>* visible;
    int iter_im;

    std::vector<cv::Point2f> normalsContour;

    void resetSprings();
    void addForceMesh(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v );

    double computeError(Vector3 sourcePoint, Vector3 targetPoint);	
		
};


/*#if defined(SOFA_EXTERN_TEMPLATE) && !defined(RegistrationForceFieldCam_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API RegistrationForceFieldCam<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API RegistrationForceFieldCam<defaulttype::Vec3fTypes>;
#endif
#endif*/


} //

} //

} // namespace sofa

#endif
