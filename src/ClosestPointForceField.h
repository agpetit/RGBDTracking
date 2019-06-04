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
//#include "KalmanFilter.h"
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
//#include <GL/glext.h>
#include <GL/glu.h>

#include <string>
#include <boost/thread.hpp>
#include "ClosestPoint.h"
#include "RGBDDataProcessing.h"
#include "MeshProcessing.h"


using namespace std;
using namespace cv;


namespace sofa {

namespace rgbdtracking {

using helper::vector;
using namespace sofa::defaulttype;

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

    typedef typename  component::interactionforcefield::LinearSpring<Real> Spring;
    typedef helper::fixed_array <unsigned int,3> tri;
	
    //typedef typename Coord::value_type real;
    typedef std::pair<int,Real> Col_Value;
    typedef std::vector< Col_Value > CompressedValue;
    typedef std::vector< CompressedValue > CompressedMatrix;
	
public:

    core::objectmodel::SingleLink<
        ClosestPointForceField<DataTypes>,
        DataIO<DataTypes>,
        BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_dataio ;
    core::objectmodel::SingleLink<
        ClosestPointForceField<DataTypes>,
        RGBDDataProcessing<DataTypes>,
        BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_rgbddataprocess ;
    core::objectmodel::SingleLink<
        ClosestPointForceField<DataTypes>,
        MeshProcessing<DataTypes>,
        BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_meshprocessing ;

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
    virtual void addKToMatrix(sofa::defaulttype::BaseMatrix *m, SReal kFactor, unsigned int &offset);

    // data accessors
    double getPotentialEnergy(const core::MechanicalParams* ,const DataVecCoord&) const { return m_potentialEnergy; }
    Real getStiffness() const{ return ks.getValue(); }
    Real getDamping() const{ return kd.getValue(); }
    void setStiffness(Real _ks){ ks.setValue(_ks); }
    void setDamping(Real _kd){ kd.setValue(_kd); }

    void draw(const core::visual::VisualParams* vparams);

    // -- Modifiers
    void clearSprings(int reserve=0)
    {
        sofa::helper::vector<Spring>& springs = *this->springs.beginEdit();
        springs.clear();
        if (reserve) {
            springs.reserve(reserve);
        }
        this->springs.endEdit();
    }
    void removeSpring(unsigned int idSpring)
    {
        if (idSpring >= (this->springs.getValue()).size()) {
            return;
        }

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
    void resetSprings() {
        this->clearSprings(l_meshprocessing->sourceVisiblePositions.getValue().size());
        for(unsigned int i = 0 ; i < l_meshprocessing->sourceVisiblePositions.getValue().size() ; i++) {
            this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());
        }
    }

protected :

    // internal data used for closest point
    size_t npoints;
    VecCoord closestPos;
    typename rgbdtracking::ClosestPoint<DataTypes> *closestpoint;

    vector<Mat>  dfdx, dfdx1;
    double m_potentialEnergy;
    Real min,max;

    /// Accumulate the spring force and compute and store its stiffness
    virtual void addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
    virtual void addSpringForceWeight(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, int ivis, const Spring& spring);
    /// Apply the stiffness, i.e. accumulate df given dx
    virtual void addSpringForceKLTA(double& potentialEnergy, VecDeriv& f, const  VecCoord& p,const VecDeriv& v, Coord& KLTtarget, int i, const Spring& spring, double coef);

    // component interface data
    Data<Real> ks;
    Data<Real> kd;
    Data<bool> projectToPlane;
    Data<sofa::helper::vector<Spring> > springs;

    //closestpoint parameters
    Data<Real> blendingFactor;
    Data<Real> outlierThreshold;
    Data<bool> rejectBorders;
    Data<bool> useContour;
    Data<bool> useVisible;
    Data<bool> useDistContourNormal;
    Data<Vector4> cameraIntrinsicParameters;

    // Number of iterations
    Data<bool> theCloserTheStiffer;
    Data<bool> useContourWeight;
    Data<int> startimage;
    Data<int> niterations;

    // KLT
    Data<int> startimageklt;
    Data<int> windowKLT;
    Data<bool> useKLTPoints;
    Data< VecCoord > targetKLTPositions;

    Eigen::Matrix3f rgbIntrinsicMatrix;
    vpKltOpencv tracker,tracker1;

    sofa::helper::vector<Vector3> mappingkltcoef;
    sofa::helper::vector<int> mappingkltind;

    VecCoord f_;//WDataRefVecDeriv f(_f);
    VecCoord x_;//RDataRefVecCoord x(_x); // not used anywhere
    VecCoord v_; // unused

    void mapKLTPointsTriangles ( helper::vector< tri > &triangles);
    void KLTPointsTo3D();

    void addForceMesh(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v );

    double computeError(Vector3 sourcePoint, Vector3 targetPoint);
		
}; // end class ClosestPointForceFields


/*#if defined(SOFA_EXTERN_TEMPLATE) && !defined(ClosestPointForceField_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API ClosestPointForceField<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API ClosestPointForceField<defaulttype::Vec3fTypes>;
#endif
#endif*/

} // rgbdtracking

} // namespace sofa

#endif
