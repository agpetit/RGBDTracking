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
	Kalmanfilter kalman;
    typedef defaulttype::ImageF DepthTypes;
	
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
    //typedef helper::kdTree<Coord> KDT;
    //typedef typename KDT::distanceSet distanceSet;
	
	double timef;
	cv::Rect rectRtt;
	
	//typedef typename Coord::value_type real;
	typedef std::pair<int,Real> Col_Value;
    typedef vector< Col_Value > CompressedValue;
    typedef vector< CompressedValue > CompressedMatrix;

    CompressedMatrix _stiffnesses;
	vpHomogeneousMatrix cMo;
	Data< VecReal > translation;
	Data< VecReal > rotation;
	
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
    void computeTargetNormals();

    vector<Mat>  dfdx;
	vector<Mat>  dfdx1;
    VecCoord closestPos;
    vector<unsigned int>  cnt;
    double m_potentialEnergy;
	
	VecCoord displ;

    Real min,max;

    /// Accumulate the spring force and compute and store its stiffness
    virtual void addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
	virtual void addStoredSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
	virtual void addSpringForceColorContour(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring);
	virtual void addSpringForceWeight(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v, int i, const Spring& spring);
	virtual void addSpringForceKLT(double& potentialEnergy, VecDeriv& f, const  VecCoord& p,const VecDeriv& v, Coord& KLTtarget, int i, const Spring& spring, double coef);
	virtual void addSpringForceKLTA(double& potentialEnergy, VecDeriv& f, const  VecCoord& p,const VecDeriv& v, Coord& KLTtarget, int i, const Spring& spring, double coef);

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
	
        typename sofa::core::objectmodel::RGBDDataProcessing<DataTypes>::SPtr rgbddataprocessing;
        typename sofa::core::objectmodel::MeshProcessing<DataTypes>::SPtr meshprocessing;
        typename sofa::core::objectmodel::RenderTextureAR<DataTypes>::SPtr rendertexturear;
        typename sofa::core::objectmodel::ClosestPoint<DataTypes>::SPtr closestpoint;
        typename sofa::core::objectmodel::DataIO<DataTypes>::SPtr dataio;
	//typename ImageConverter<DataTypes,DepthTypes>::SPtr imconv;

    Data<int> viewportWidth;
    Data<int> viewportHeight;
		
    VecCoord tpos;
	
	Data<Vector4> barycenter;

    // source mesh data
    Data< helper::vector< tri > > sourceTriangles;
    Data< VecCoord > sourceNormals;
	Data< VecCoord > sourceSurfacePositions;
    Data< VecCoord > sourceSurfaceNormals;
	Data< VecCoord > sourceSurfaceNormalsM;
	Data< VecCoord > sourceContourPositions;
    //vector< distanceSet >  closestSource; // CacheSize-closest target points from source
	//vector< distanceSet >  closestSourceContour; // CacheSize-closest source points from target

    //vector< bool > sourceBorder;
    //vector< bool > sourceIgnored;  // flag ignored vertices
	vector< bool > sourceVisible;  // flag ignored vertices
	vector< bool > sourceSurface;
    //vector< bool > targetIgnored;  // flag ignored vertices
	vector< bool > targetBackground;  // flag ignored vertices

	std::vector<int> indices;
	std::vector<int> indicesTarget;
	std::vector<int> indicesVisible;
	std::vector<int> sourceSurfaceMapping;
	
	VecCoord f_ ;       //WDataRefVecDeriv f(_f);
    VecCoord  x_ ;			//RDataRefVecCoord x(_x);
    VecCoord v_;	

    // target point cloud data
	Data< VecCoord > sourcePositions;
    Data< VecCoord > targetPositions;
	Data< VecCoord > targetGtPositions;
    Data< VecCoord > targetNormals;
	Data< VecCoord > targetKLTPositions;
	Data< VecCoord > targetCCDPositions;
    Data< helper::vector< tri > > targetTriangles;
    //vector< distanceSet >  closestTarget; // CacheSize-closest source points from target
    vector< bool > targetBorder;
	vector < double > targetWeights;
	Data< VecCoord > targetContourPositions;
    vector < double > sourceWeights;
	vector < double > combinedWeights;
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
	
	vpImage<unsigned char> vpI ; 
    vpDisplayX display;


	// Number of iterations
	Data<int> niterations;
	Data<int> nimages;
	Data<Real> sigmaWeight;
	int npasses;
	
	int sock_sender, err;
	
	Data<int> samplePCD;
	Data<int> offsetX, offsetY;
	Data<int> borderThdPCD;
	Data<int> borderThdSource;
	Data<int> windowKLT;
	Data<bool> useDistContourNormal;
		//**************************************************************
	point_struct center_pos;
	//**************************************************************
	
	// Paths
	Data<std::string> inputPath;
	Data<std::string> outputPath;
	Data<std::string> dataPath;
	Data<std::string> ipad;

    cv::Mat foreground;
	cv::Mat foregroundbin;
	bool pcl;
	bool disp;
	Data<bool> useContour;
	Data<bool> useVisible;
	Data<bool> useRealData;
	Data<bool> useGroundTruth;
	Data<bool> useIntensity;
	Data<bool> useKLTPoints;
	Data<bool> useCCD;
	Data<bool> generateSynthData;
	Data<bool> useSensor;
	Data<int> sensorType;
	Data<bool> useMassSpring;
    Data<bool> showStrainsPerElement;
	Data<bool> drawSource;
	Data<bool> drawTarget;
    Data<bool> drawContour;
	
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
	
	cv::Mat rtd;
    vpKltOpencv tracker,tracker1;

	sofa::helper::vector<Vector3> mappingkltcoef;
	sofa::helper::vector<int> mappingkltind;
	
	CCD ccd;
	std::vector<pointCCD> pointsCCD;
	std::vector<pointCCD> pointsCCDmin;
	std::vector<cv::Point2f> normalsContour;
	void initCCD();
	void updateCCD();

    std::vector<int> source2target_;
    std::vector<int> target2source_;
	std::vector<int> source2target_distances_;
    std::vector<int> target2source_distances_;
    pcl::CorrespondencesPtr correspondences_;
	std::vector<int> distances_;
	double determineErrorICP();
	
	void resetSprings();
	void setViewPoint();
	void addForceMesh(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v );
	//void addForceSoft(const core::MechanicalParams* /*mparams*/,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v );

    double computeError(Vector3 sourcePoint, Vector3 targetPoint);
	void initCamera();
	
    void mapKLTPointsTriangles ( helper::vector< tri > &triangles);
    void addPointInTriangle ( const int triangleIndex, const Real* baryCoords );
    void KLTPointsTo3D();
	void CCDPointsTo3D();
	
		
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
