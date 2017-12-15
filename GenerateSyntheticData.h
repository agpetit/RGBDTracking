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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_RGBDTRACKING_GENERATESYNTHETICDATA_H
#define SOFA_RGBDTRACKING_GENERATESYNTHETICDATA_H

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
#include <sofa/helper/kdTree.inl>

#include <set>

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


#include <opencv/cv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/legacy/legacy.hpp>

#include <visp/vpIoTools.h>
#include <visp/vpImageIo.h>
#include <visp/vpParseArgv.h>

#include <iostream>
#include <string>
#include <map>
//#include <XnCppWrapper.h>
#include <opencv2/opencv.hpp>
#include <boost/thread.hpp>
#include <sys/times.h>

#include "luaconfig.h"

#include "segmentation.h"

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
using cimg_library::CImg;


template<class DataTypes, class _ImageTypes>
class GenerateSyntheticDataInternalData
{
public:
};

template<class DataTypes, class _ImageTypes>
class GenerateSyntheticData : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(GenerateSyntheticData,DataTypes,_ImageTypes),SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherit;
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
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<typename DataTypes::VecCoord> DataVecCoord;
    typedef Data<typename DataTypes::VecDeriv> DataVecDeriv;

    typedef core::behavior::MechanicalState<DataTypes> MechanicalState;
    enum { N=DataTypes::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;

    typedef typename interactionforcefield::LinearSpring<Real> Spring;
    typedef helper::fixed_array <unsigned int,3> tri;
    typedef helper::kdTree<Coord> KDT;
    typedef typename KDT::distanceSet distanceSet;

public:
    GenerateSyntheticData(core::behavior::MechanicalState<DataTypes> *mm = NULL);
    virtual ~GenerateSyntheticData();

    core::behavior::MechanicalState<DataTypes>* getObject() { return this->mstate; }
	
	static std::string templateName(const GenerateSyntheticData<DataTypes,DepthTypes >* = NULL) { return DataTypes::Name()+ std::string(",")+DepthTypes::Name();    }
    virtual std::string getTemplateName() const    { return templateName(this);    }

    const sofa::helper::vector< Spring >& getSprings() const {return springs.getValue();}

    // -- ForceField interface
    void reinit();
    void init();
    void addForce(const core::MechanicalParams* /*mparams*/,DataVecDeriv& f , const DataVecCoord& x , const DataVecDeriv& v);
    void addDForce(const core::MechanicalParams* mparams ,DataVecDeriv&   df , const DataVecDeriv&   dx);
    double getPotentialEnergy(const core::MechanicalParams* ,const DataVecCoord&) const { return m_potentialEnergy; }
    void addKToMatrix( const core::MechanicalParams* mparams,const sofa::core::behavior::MultiMatrixAccessor* matrix);

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
	void writeImages();
	void writeImagesSynth();
	void writeData();

    vector<Mat>  dfdx;
    VecCoord closestPos;
    vector<unsigned int>  cnt;
    double m_potentialEnergy;

    Real min,max;

    /// Accumulate the spring force and compute and store its stiffness
    virtual void addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
	virtual void addSpringForce1(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
	virtual void addSpringForce2(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
	virtual void addSpringForce3(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
    virtual void addSpringForceSoft(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, int j, const Spring& spring, Real weightP);

	/// Apply the stiffness, i.e. accumulate df given dx
    virtual void addSpringDForce(VecDeriv& df,const  VecDeriv& dx, int i, const Spring& spring, double kFactor, double bFactor);

    Data<Real> ks;
    Data<Real> kd;
    Data<unsigned int> cacheSize;
    Data<Real> blendingFactor;
    Data<Real> outlierThreshold;
    Data<Real> normalThreshold;
    Data<bool> projectToPlane;
    Data<bool> rejectBorders;
    Data<bool> rejectOutsideBbox;
    defaulttype::BoundingBox targetBbox;
    Data<sofa::helper::vector<Spring> > springs;
	
    VecCoord tpos;

    // source mesh data
    Data< helper::vector< tri > > sourceTriangles;
    Data< VecCoord > sourceNormals;
	Data< VecCoord > sourceSurfacePositions;
    Data< VecCoord > sourceSurfaceNormals;
	Data< VecCoord > sourceSurfaceNormalsM;
    vector< distanceSet >  closestSource; // CacheSize-closest target points from source
    vector< Real > cacheDist;	vector< Real > cacheDist2; VecCoord previousX; // storage for cache acceleration
    KDT sourceKdTree;
    vector< bool > sourceBorder;
    vector< bool > sourceIgnored;  // flag ignored vertices
	vector< bool > sourceVisible;  // flag ignored vertices
	vector< bool > sourceSurface;
    vector< bool > targetIgnored;  // flag ignored vertices
	vector< bool > targetBackground;  // flag ignored vertices
    void initSource(); // built k-d tree and identify border vertices
	void initSourceVisible(); // built k-d tree and identify border vertices
	void initSourceSurface(); // built k-d tree and identify border vertices
	void updateSourceSurface(); // built k-d tree and identify border vertices

	
	std::vector<int> indices;
	std::vector<int> indicesVisible;
	std::vector<int> sourceSurfaceMapping;

	
	Data< VecCoord > sourceContourPositions;


    // target point cloud data
    Data< VecCoord > targetPositions;
	Data< VecCoord > targetGtPositions;
    Data< VecCoord > targetNormals;
    Data< helper::vector< tri > > targetTriangles;
    vector< distanceSet >  closestTarget; // CacheSize-closest source points from target
    KDT targetKdTree;
    vector< bool > targetBorder;
	vector < double > targetWeights;
    void initTarget();  // built k-d tree and identify border vertices
	
	int ind;
	
	Data< VecCoord > sourceVisiblePositions;

    Data<float> showArrowSize;
    Data<int> drawMode; //Draw Mode: 0=Line - 1=Cylinder - 2=Arrow
    Data<bool> drawColorMap;
    Data<bool> theCloserTheStiffer;
	
	//SoftKinetic softk;
	cv::Mat depth;	
	cv::Mat color;
	cv::Mat color_1;
	cv::Mat depthMap;
	cv::Mat silhouetteMap;


	// Number of iterations
	int niterations;
	int nimages;

    segmentation seg;

    cv::Mat foreground;
	bool pcl;
	bool disp;
	Data<bool> useContour;
	Data<bool> useVisible;
	Data<bool> useRealData;
	Data<bool> useGroundTruth;
	Data<bool> generateSynthData;
	Data<bool> useSensor;
	
	int ntargetcontours;
 
    std::vector<cv::Mat*> listimg;
    std::vector<cv::Mat*> listimgseg;
    std::vector<cv::Mat*> listdepth;
    std::vector<cv::Mat*> listrtt;
	std::vector<std::vector<Vec3d>*> listpcd;
	std::vector<std::vector<bool>*> listvisible;

    cv::Mat* imgl;
    cv::Mat* imglsg;
    cv::Mat* depthl;
	cv::Mat* rtt;
	std::vector<Vec3d>* pcd;
	std::vector<bool>* visible;
	
	
	int iter_im;

	double timeOverall;
	double timeTotal;
	double timeInt;
	double timei;
	
	double errorGroundTruth;

	ofstream timeFile;
	
	void resetSprings();
    void detectBorder(vector<bool> &border,const helper::vector< tri > &triangles);
	void setViewPointData();
	void getSourceVisible();

	void generateData(const core::MechanicalParams* /*mparams*/,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v );
  //
  // initialize paramters
  //
float sigma_p2;
float sigma_inf;
float sigma_factor;
float d_02;
float *h_A;
int pitchA;
	
		   
void findRTfromS(const float* h_Xc, const float* h_Yc, const float* h_S, float* h_R, float* h_t);

void printRT(const float* R, const float* t);

void cloud2dataC(const pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud,
		float **X, int &Xsize );
void cloud2data(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
		float **X, int &Xsize );
};


/*#if defined(SOFA_EXTERN_TEMPLATE) && !defined(GenerateSyntheticData_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API GenerateSyntheticData<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API GenerateSyntheticData<defaulttype::Vec3fTypes>;
#endif
#*/


} //

} //

} // namespace sofa

#endif
