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

#ifndef SOFA_RGBDTRACKING_DATAGENERATION_H
#define SOFA_RGBDTRACKING_DATAGENERATION_H


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
#include <sofa/component/topology/TopologyData.h>
#include <visp/vpKltOpencv.h>
#include <sofa/helper/kdTree.inl>

#include <set>
#include <RGBDTracking/config.h>
#include <visp/vpDisplayX.h>
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
#include <visp/vpMatrix.h>

#include <iostream>
#include <string>
#include <map>
//#include <XnCppWrapper.h>
#include <opencv2/opencv.hpp>
#include <boost/thread.hpp>

#include "luaconfig.h"
#include "p_helper.h"
#include "MeshProcessing.h"
#include "RenderTextureAR.h"


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
using cimg_library::CImg;
using namespace sofa::component::topology;


template<class DataTypes, class _ImageTypes>
class DataGenerationInternalData
{
public:
};

template<class DataTypes, class _ImageTypes>
class DataGeneration : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(DataGeneration,DataTypes,_ImageTypes),SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

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
	//sofa::component::forcefield::KalmanFilter kalman;
	Kalmanfilter kalman;
	
	int npoints;

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
    typedef helper::kdTree<Coord> KDT;
    typedef typename KDT::distanceSet distanceSet;
	
	double timef;
	cv::Rect rectRtt;
	
	//typedef typename Coord::value_type real;
	typedef std::pair<int,Real> Col_Value;
    typedef vector< Col_Value > CompressedValue;
    typedef vector< CompressedValue > CompressedMatrix;

    CompressedMatrix _stiffnesses;
	vpHomogeneousMatrix cMo;

public:
    DataGeneration(core::behavior::MechanicalState<DataTypes> *mm = NULL);
    virtual ~DataGeneration();

    core::behavior::MechanicalState<DataTypes>* getObject() { return this->mstate; }
	
	static std::string templateName(const DataGeneration<DataTypes,DepthTypes >* = NULL) { return DataTypes::Name()+ std::string(",")+DepthTypes::Name();    }
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
	void writeData();
    void computeTargetNormals();
	vector<Mat>  dfdx;
    VecCoord closestPos;
    vector<unsigned int>  cnt;
    double m_potentialEnergy;
	
	VecCoord displ;

    Real min,max;

    /// Accumulate the spring force and compute and store its stiffness
    virtual void addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
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
	
	typename MeshProcessing<DataTypes>::SPtr meshprocessing;
	typename RenderTextureAR<DataTypes>::SPtr rendertexturear;
	
	Data<Vector3> barycenter;

    // source mesh data
    Data< helper::vector< tri > > sourceTriangles;
    Data< VecCoord > sourceNormals;
	Data< VecCoord > sourceSurfacePositions;
    Data< VecCoord > sourceSurfaceNormals;
	Data< VecCoord > sourceSurfaceNormalsM;
	Data< VecCoord > sourceContourPositions;
    vector< bool > sourceBorder;
    vector< bool > sourceIgnored;  // flag ignored vertices
	vector< bool > sourceVisible;  // flag ignored vertices
	vector< bool > sourceSurface;
    void initSource(); // built k-d tree and identify border vertices
	void initSourceVisible(); // built k-d tree and identify border vertices
	void initSourceSurface(); // built k-d tree and identify border vertices

	std::vector<int> indices;
	std::vector<int> indicesTarget;
	std::vector<int> indicesVisible;
	std::vector<int> sourceSurfaceMapping;
	
	VecCoord f_ ;       //WDataRefVecDeriv f(_f);
    VecCoord  x_ ;			//RDataRefVecCoord x(_x);
    VecCoord v_;	

    // target point cloud data
	Data< VecCoord > sourcePositions;
	void normalizeWeights();
	
	std::vector<double> vonmisesstressGt;
	std::vector<double> elasticstrainsGt;
	std::vector<double> plasticstrainsGt;
	std::vector<double> totalstrainsGt;
	std::vector<double> elasticstrainsnodeGt;
	std::vector<double> plasticstrainsnodeGt;
	std::vector<double> totalstrainsnodeGt;

	
	int ind;
	
	Data< VecCoord > sourceVisiblePositions;

    Data<float> showArrowSize;
    Data<int> drawMode; //Draw Mode: 0=Line - 1=Cylinder - 2=Arrow
    Data<bool> drawColorMap;
    Data<bool> theCloserTheStiffer;
	
	//SoftKinetic softk;
	cv::Mat depth,depth_1, depthrend, depth00, depth01;	
	cv::Mat color, ir, ig, ib, gray;
	cv::Mat color_1,color_2, color_3, color_4, color_5, color_init;
	cv::Mat depthMap;
	cv::Mat silhouetteMap;
	
	// Number of iterations
	Data<int> niterations;
	Data<int> nimages;
	int npasses;
	
	Data<int> borderThdPCD;
	Data<int> borderThdSource;
	Data<int> windowKLT;

	
	// Paths
	Data<std::string> inputPath;
	Data<std::string> outputPath;
	Data<std::string> dataPath;

    cv::Mat foreground;
	cv::Mat foregroundbin;
	bool pcl;
	bool disp;
	Data<bool> useContour;
	Data<bool> useVisible;
	Data<bool> useRealData;
	Data<bool> useGroundTruth;
	Data<bool> generateSynthData;
	Data<bool> useSensor;
	Data<int> sensorType;
	Data<bool> useMassSpring;
    Data<bool> showStrainsPerElement;

	
	Data<Vector4> cameraIntrinsicParameters;
	Eigen::Matrix3f rgbIntrinsicMatrix;
	
	Data<Real> alphaIntensity;
	Data<Real> visibilityThreshold;
	Data<helper::vector<Real> > plasticStrainsN; ///< one plastic strain per element
	Data<helper::vector<Real> > elasticStrainsN; ///< one plastic strain per element
	Data<helper::vector<Real> > totalStrainsN; ///< one plastic strain per element
	Data<helper::vector<Real> > vonMisesStress; ///< one plastic strain per element
	Data<helper::vector<Real> > elasticStrainsPerNode; ///< one plastic strain per element
	Data<helper::vector<Real> > plasticStrainsPerNode; ///< one plastic strain per element
	Data<helper::vector<Real> > totalStrainsPerNode; ///< one plastic strain per element


	int ntargetcontours;
 
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


    cv::Mat* imgl;
	cv::Mat* imgklt;
    cv::Mat* imglsg;
    cv::Mat* depthl;
	cv::Mat* rtt;
	std::vector<Vec3d>* pcd;
	std::vector<bool>* visible;
	std::vector<double>* vm;
	std::vector<double>* ps;
	std::vector<double>* es;
	std::vector<double>* ts;
	std::vector<double>* esnode;
	std::vector<double>* psnode;
	std::vector<double>* tsnode;


	int iter_im;

	double timeOverall;
	double timeSeg;
	double timeRigid;
	double timeAddforce;
	double timeTotal;
	double timeSourceContour;
	double timeResolution;
	double timeInt;
	double timei,timeii;
	
	int timer;
	double timeSecondPass;
	double timeFirstPass;
	
	double timeOverall1;
	
	double errorMatching;
	double errorGroundTruth, errorGroundTruth_vM, errorGroundTruth_eS,errorGroundTruth_pS,errorGroundTruth_tS, errorGroundTruth_eSN, errorGroundTruth_pSN, errorGroundTruth_tSN, meanGt_vM, meanGt_eS, meanGt_pS, meanGt_tS, meanGt_eSN, meanGt_pSN, meanGt_tSN;

	ofstream timeFile;
	ofstream errorFile;
	ofstream errorGtFile;
	ifstream correspMassSpringFile;

	Eigen::Matrix4f transformation_matrix;
	
	cv::Mat rtd;

	void resetSprings();
    void detectBorder(vector<bool> &border,const helper::vector< tri > &triangles);
	void setViewPoint();
	void setViewPointData();
    void renderToTextureD(cv::Mat &_rtt);
	void renderToTextureDepth(cv::Mat &_rttd, cv::Mat &_rttdepth);
	void generateData(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v );
	//void addForceSoft(const core::MechancalParams* /*mparams*/,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v );

		
};


/*#if defined(SOFA_EXTERN_TEMPLATE) && !defined(DataGeneration_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API DataGeneration<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API DataGeneration<defaulttype::Vec3fTypes>;
#endif
#endif*/


} //

} //

} // namespace sofa

#endif
