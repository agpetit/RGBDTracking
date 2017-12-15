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

#ifndef SOFA_RGBDTRACKING_CLOSESTPOINT_H
#define SOFA_RGBDTRACKING_CLOSESTPOINT_H

#include <RGBDTracking/config.h>
#include <image/ImageTypes.h>
#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>

#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/helper/kdTree.inl>

#include <vector>
#include <opencv/cv.h>
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
using namespace sofa::component::topology;


template<class DataTypes>
class ClosestPointInternalData
{
public:
};

template<class DataTypes>
class ClosestPoint : public sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ClosestPoint,DataTypes),sofa::core::objectmodel::BaseObject);

    typedef sofa::core::objectmodel::BaseObject Inherit;
	
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
	
	typename core::behavior::MechanicalState<DataTypes> *mstate;
    
	enum { N=DataTypes::spatial_dimensions };
    typedef defaulttype::Mat<N,N,Real> Mat;

    typedef helper::fixed_array <unsigned int,3> tri;
    typedef helper::kdTree<Coord> KDT;
    typedef typename KDT::distanceSet distanceSet;
	
	double timef;
	cv::Rect rectRtt;
	
	int timer;
	
	//typedef typename Coord::value_type real;
	typedef std::pair<int,Real> Col_Value;
    typedef vector< Col_Value > CompressedValue;
    typedef vector< CompressedValue > CompressedMatrix;
	
	Data< double > errorfunction;

public:
    ClosestPoint();
    virtual ~ClosestPoint();
	
	/*static std::string templateName(const ClosestPoint<DataTypes >* = NULL) { return DataTypes::Name()+ std::string(",");    }
    virtual std::string getTemplateName() const    { return templateName(this);    }*/
	
    // -- baseobject interface
    void init();

    void findCorrespondences (pcl::PointCloud<pcl::PointXYZRGB>::Ptr source, pcl::PointCloud<pcl::PointXYZRGB>::Ptr target, std::vector<int>& correspondences, std::vector<int>& distances);
    void filterCorrespondences();
	void updateClosestPoints();
	void updateClosestPointsSoft();
    void updateClosestPointsPCL();
	void updateClosestPointsContours();
	void updateClosestPointsVisibleContours();
	void updateClosestPointsContoursNormals();

    vector<Mat>  dfdx;
	vector<Mat>  dfdx1;
    VecCoord closestPos;
    vector<unsigned int>  cnt;
	
	VecCoord displ;
    Real min,max;

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
    // source mesh data
    Data< helper::vector< tri > > sourceTriangles;
    Data< VecCoord > sourceNormals;
	Data< VecCoord > sourceSurfacePositions;
    Data< VecCoord > sourceSurfaceNormals;
	Data< VecCoord > sourceSurfaceNormalsM;
	Data< VecCoord > sourceContourPositions;
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
    Data< helper::vector< tri > > targetTriangles;
    vector< distanceSet >  closestTarget; // CacheSize-closest source points from target
    KDT targetKdTree;
    vector< bool > targetBorder;
	vector < double > targetWeights;
	Data< VecCoord > targetContourPositions;
    KDT targetContourKdTree;
    vector < double > sourceWeights;
	vector < double > combinedWeights;
	void initTarget();  // built k-d tree and identify border vertices
	void initTargetContour();  // built k-d tree and identify border vertices
	void normalizeWeights();

	int ind;
	Data< VecCoord > sourceVisiblePositions;

	// Number of iterations
	Data<int> niterations;
	Data<int> nimages;
	Data<Real> sigmaWeight;
	int npasses;
		
	Data<int> borderThdPCD;
	Data<int> borderThdSource;
	Data<bool> useDistContourNormal;

	bool pcl;
	bool disp;
	Data<bool> useContour;
	Data<bool> useVisible;
	
	Data<Vector4> cameraIntrinsicParameters;
	Eigen::Matrix3f rgbIntrinsicMatrix;
	
	int ntargetcontours;
	int iter_im;
	
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr source;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr source0;
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr target;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetP;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetGt;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr source_registered;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr source_registered0;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetContour;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetPointCloud;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr sourceSurfacePointCloud;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr sourceSurfacePointCloud_registered;

	Eigen::Matrix4f transformation_matrix;
	std::vector<cv::Point2f> normalsContour;

    std::vector<int> source2target_;
    std::vector<int> target2source_;
	std::vector<int> source2target_distances_;
    std::vector<int> target2source_distances_;
    pcl::CorrespondencesPtr correspondences_;
	std::vector<int> distances_;
	
	/*VecCoord getTargetPositions(VecCoord& _targetpositions){targetPositions.setValue(_targetpositions);}
	VecCoord setSourceVisiblePositions(VecCoord& _sourcevisiblepositions){sourceVisiblePositions.setValue(_sourcevisiblepositions);}
	VecCoord setSourceContourPositions(VecCoord& _sourcecontourpositions){sourceContourPositions.setValue(_sourcecontourpositions);}
	VecCoord setSourceSurfacePositions(VecCoord& _sourcesurfacepositions){sourcesurfacePositions.setValue(_sourcesurfacepositions);}
	VecCoord setTargetContourPositions(VecCoord& _targetcontourpositions){targetContourPositions.setValue(_targetcontourpositions);}*/
	
	vector< distanceSet > getClosestSource(){return closestSource;}
	vector< distanceSet > getClosestTarget(){return closestTarget;}

	vector< bool > getSourceIgnored(){return sourceIgnored;}
	vector< bool > getTargetIgnored(){return targetIgnored;}
	
	vector< int > getIndices(){return indices;}
	vector< int > getIndicesTarget(){return indicesTarget;}
	
    void detectBorder(vector<bool> &border,const helper::vector< tri > &triangles);


};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(ClosestPoint_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API ClosestPoint<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API ClosestPoint<defaulttype::Vec3fTypes>;
#endif
#endif


} //

} //

} // namespace sofa

#endif
