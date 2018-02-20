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

#define SOFA_RGBDTRACKING_DATAGENERATION_CPP

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/Mapping.inl>
#include <sofa/simulation/Simulation.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/gui/BaseGUI.h>
#include <sofa/gui/BaseViewer.h>
#include <sofa/gui/GUIManager.h>
#include <iostream>
#include <map>

#ifdef USING_OMP_PRAGMAS
    #include <omp.h>
#endif


#include <SofaLoader/MeshObjLoader.h>
#include <limits>
#include <set>
#include <iterator>
#include <sofa/helper/gl/Color.h>

#ifdef Success
  #undef Success
#endif

#include <pcl/common/common_headers.h>
#include <pcl/point_cloud.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/point_types.h>
#include <pcl/impl/point_types.hpp>
#include <pcl/features/normal_3d.h>
#include <pcl/common/projection_matrix.h>
#include <pcl/io/pcd_io.h>
#include <pcl/common/io.h>
#include <pcl/registration/transforms.h>
#include <pcl/search/kdtree.h>
#include <pcl/registration/transformation_estimation_svd.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/correspondence_rejection_sample_consensus.h>
#include <pcl/common/transforms.h>

#include <algorithm> 
#include "DataGeneration.h"

using std::cerr;
using std::endl;

namespace sofa
{

namespace component
{

namespace forcefield
{

    using namespace sofa::defaulttype;

      SOFA_DECL_CLASS(DataGeneration)

      // Register in the Factory
      int DataGenerationClass = core::RegisterObject("Compute forces based on closest points from/to a target surface/point set")
    #ifndef SOFA_FLOAT
        .add< DataGeneration<Vec3dTypes,ImageUC> >()
        .add< DataGeneration<Vec3dTypes,ImageUS> >()
        .add< DataGeneration<Vec3dTypes,ImageF> >()
    #endif
    #ifndef SOFA_DOUBLE
        .add< DataGeneration<Vec3fTypes,ImageUC> >()
        .add< DataGeneration<Vec3fTypes,ImageUS> >()
        .add< DataGeneration<Vec3fTypes,ImageF> >()
    #endif
    ;

    #ifndef SOFA_FLOAT
      template class SOFA_RGBDTRACKING_API DataGeneration<Vec3dTypes,ImageUC>;
      template class SOFA_RGBDTRACKING_API DataGeneration<Vec3dTypes,ImageUS>;
      template class SOFA_RGBDTRACKING_API DataGeneration<Vec3dTypes,ImageF>;
    #endif
    #ifndef SOFA_DOUBLE
      template class SOFA_RGBDTRACKING_API DataGeneration<Vec3fTypes,ImageUC>;
      template class SOFA_RGBDTRACKING_API DataGeneration<Vec3fTypes,ImageUS>;
      template class SOFA_RGBDTRACKING_API DataGeneration<Vec3fTypes,ImageF>;
    #endif

using namespace helper;


template <class DataTypes, class DepthTypes>
DataGeneration<DataTypes, DepthTypes>::DataGeneration(core::behavior::MechanicalState<DataTypes> *mm )
    : Inherit(mm)
    , depthImage(initData(&depthImage,DepthTypes(),"depthImage","depth map"))
    //, depthTransform(initData(&depthTransform, TransformType(), "depthTransform" , ""))
    , image(initData(&image,ImageTypes(),"image","image"))
	//, transform(initData(&transform, TransformType(), "transform" , ""))
	, cameraIntrinsicParameters(initData(&cameraIntrinsicParameters,Vector4(),"cameraIntrinsicParameters","camera parameters"))
    , ks(initData(&ks,(Real)0.0,"stiffness","uniform stiffness for the all springs."))
    , kd(initData(&kd,(Real)0.0,"damping","uniform damping for the all springs."))
    , cacheSize(initData(&cacheSize,(unsigned int)10,"cacheSize","number of closest points used in the cache to speed up closest point computation."))
    , normalThreshold(initData(&normalThreshold,(Real)0,"normalThreshold","suppress outliers when normal.closestPointNormal < threshold."))
    , projectToPlane(initData(&projectToPlane,false,"projectToPlane","project closest points in the plane defined by the normal."))
    , rejectBorders(initData(&rejectBorders,false,"rejectBorders","ignore border vertices."))
    , springs(initData(&springs,"spring","index, stiffness, damping"))
	, sourceSurfacePositions(initData(&sourceSurfacePositions,"sourceSurface","Points of the surface of the source mesh."))
	, sourcePositions(initData(&sourcePositions,"sourcePositions","Points of the surface of the source mesh."))
	, sourceTriangles(initData(&sourceTriangles,"sourceTriangles","Triangles of the source mesh."))
    , sourceNormals(initData(&sourceNormals,"sourceNormals","Normals of the source mesh."))
	, sourceSurfaceNormals(initData(&sourceSurfaceNormals,"sourceSurfaceNormals","Normals of the surface of the source mesh."))
	, showStrainsPerElement(initData(&showStrainsPerElement, false, "showStrainsPerElement", "  "))
	, plasticStrainsN(initData(&plasticStrainsN, "plasticStrains", "plastic strain per element"))
	, elasticStrainsN(initData(&elasticStrainsN, "elasticStrains", "elastic strain per element"))
	, totalStrainsN(initData(&totalStrainsN, "totalStrains", "total strain per element"))
	, plasticStrainsPerNode(initData(&plasticStrainsPerNode, "plasticStrainsPerNode", "plastic strain per node"))
	, elasticStrainsPerNode(initData(&elasticStrainsPerNode, "elasticStrainsPerNode", "elastic strain per node"))
	, totalStrainsPerNode(initData(&totalStrainsPerNode, "totalStrainsPerNode", "total strain per node"))
	, vonMisesStress(initData(&vonMisesStress, "vonMisesStress", "vonmisesstress per element"))
    , barycenter(initData(&barycenter,"barycenter","Barycenter of the mesh."))
	, showArrowSize(initData(&showArrowSize,0.01f,"showArrowSize","size of the axis."))
	, drawMode(initData(&drawMode,0,"drawMode","The way springs will be drawn:\n- 0: Line\n- 1:Cylinder\n- 2: Arrow."))
    , drawColorMap(initData(&drawColorMap,true,"drawColorMap","Hue mapping of distances to closest point"))
    , theCloserTheStiffer(initData(&theCloserTheStiffer,false,"theCloserTheStiffer","Modify stiffness according to distance"))
	, useContour(initData(&useContour,false,"useContour","Emphasize forces close to the target contours"))
	, useVisible(initData(&useVisible,true,"useVisible","Use the vertices of the viisible surface of the source mesh"))
	, visibilityThreshold(initData(&visibilityThreshold,(Real)0.001,"visibilityThreshold","Threshold to determine visible vertices"))
	, alphaIntensity(initData(&alphaIntensity,(Real)0.00004,"alphaIntensity","Weight of intensity features"))
	, useRealData(initData(&useRealData,true,"useRealData","Use real data"))
	, useGroundTruth(initData(&useGroundTruth,false,"useGroundTruth","Use the vertices of the visible surface of the source mesh"))
	//, useKalman(initData(&useKalman,false,"useKalman","Use the Kalman filter"))
	, sensorType(initData(&sensorType, 0,"sensorType","Type of the sensor"))
	, generateSynthData(initData(&generateSynthData, false,"generateSynthData","Generate synthetic data"))
	, niterations(initData(&niterations,3,"niterations","Number of iterations in the tracking process"))
	, nimages(initData(&nimages,1500,"nimages","Number of images to read"))
	, borderThdPCD(initData(&borderThdPCD,4,"borderThdPCD","border threshold on the target silhouette"))
	, borderThdSource(initData(&borderThdSource,7,"borderThdSource","border threshold on the source silhouette"))
	, inputPath(initData(&inputPath,"inputPath","Path for data readings",false))
	, outputPath(initData(&outputPath,"outputPath","Path for data writings",false))
{
	//softk.init();
	this->addAlias(&depthImage, "depthImage");
	depthImage.setGroup("depthImage");
	depthImage.setReadOnly(true); 
	this->addAlias(&image, "image");
	image.setGroup("image");
	image.setReadOnly(true); 
  	//niterations.setValue(2);
	nimages = 1500;
	pcl = false;
	disp = false;
	iter_im = 0;
	timeTotal = 0;
	color = cv::Mat::zeros(240,320, CV_8UC3);
    color_1 = color.clone();
    color_2 = color_1.clone();
    color_3 = color_2.clone();
 color_4 = color_3.clone(); 
 color_5 = color_4.clone(); 
	depth = cv::Mat::zeros(240,320,CV_32FC1);
	depth_1 = depth;
    listimg.resize(0);
    listimgseg.resize(0);
	listimgklt.resize(0);
    listdepth.resize(0);
	listrtt.resize(0);

	listrttstress.resize(0);
	listrttstressplast.resize(0);
	listpcd.resize(0);
	listvisible.resize(0);
	timef = 0;//(double)getTickCount();
	timei = 0;
	timeOverall = 0;
	timeResolution = 0;
	timer = 0;
	timeSecondPass = 0;
	timeFirstPass = 0;
	timeOverall1 = 0;
    // Tracker parameters
}

template <class DataTypes, class DepthTypes>
DataGeneration<DataTypes, DepthTypes>::~DataGeneration()
{
}


template <class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::reinit()
{

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
	this->clearSprings(x.size());
	
    for(unsigned int i=0;i<x.size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());	

}


template <class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::init()
{
	
	//pcddataprocessing.init();
	//**********************************************************************************

    this->Inherit::init();
    core::objectmodel::BaseContext* context = this->getContext();

    if(!(this->mstate)) this->mstate = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());

    // Get source normals
    if(!sourceNormals.getValue().size()) serr<<"normals of the source model not found"<<sendl;

    // add a spring for every input point
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();			//RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
	this->clearSprings(x.size());
	
	npoints = x.size();
	
    for(unsigned int i=0;i<x.size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());	
	
			//initCamera();
    bool opt_device = false;
    bool opt_display = true;
    bool use_cuda = true;
    bool opt_click_allowed = true;
    int start_image = 0;

    /*listimg.resize(0);
    listimgseg.resize(0);
    listdepth.resize(0);*/
	
	cv::Rect ROI(160, 120, 320, 240);

        sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());

        root->get(meshprocessing);
        root->get(rendertexturear);
	
	meshprocessing->visibilityThreshold.setValue(visibilityThreshold.getValue());
	meshprocessing->borderThdSource.setValue(borderThdSource.getValue());

	if (showStrainsPerElement.getValue())
		npasses = 2;
	else npasses = 1;
	

}

template <class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::resetSprings()
{
	
this->clearSprings(sourceVisiblePositions.getValue().size());	
for(unsigned int i=0;i<sourceVisiblePositions.getValue().size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());	
	
}


template<class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::initSource()
{
    // build k-d tree
    const VecCoord&  p = this->mstate->read(core::ConstVecCoordId::position())->getValue();
		
}

template<class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::initSourceSurface()
{
    // build k-d tree
    const VecCoord&  p = this->mstate->read(core::ConstVecCoordId::position())->getValue();
	
	ReadAccessor< Data< VecCoord > > pSurface(sourceSurfacePositions);
	ReadAccessor< Data< VecCoord > > nSurface(sourceSurfaceNormals);
	VecCoord nSurfaceM;
		
	sourceSurface.resize(p.size());
	Vector3 pos;
	Vector3 col;
		nSurfaceM.resize(p.size());

	for (unsigned int i=0; i<p.size(); i++)
	{
			sourceSurface[i] = false;
			for (unsigned int j=0; j<pSurface.size(); j++)
			{
				Real dist=(p[i]-pSurface[j]).norm();
				if ((float)dist < 0.001)
				{
					sourceSurface[i] = true;
					nSurfaceM[i] = nSurface[j];
				}
			}
	}
			
	sourceSurfaceNormalsM.setValue(nSurfaceM);
	
}


template<class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::initSourceVisible()
{
    // build k-d tree
	
	const VecCoord&  p = sourceVisiblePositions.getValue();
	
}

void writePCDToFile(string path, std::vector<Vec3d>& pcd, std::vector<bool>& visible)
{

//create the file stream
ofstream file(path.c_str(), ios::out | ios::binary );
double dvalue0,dvalue1,dvalue2;

for (int k = 0; k < pcd.size(); k++)
{
dvalue0 = pcd[k][0];
dvalue1 = pcd[k][1];
dvalue2 = pcd[k][2];

	file << dvalue0;	
	file << "\t";
	file << dvalue1;
	file << "\t";
	file << dvalue2;	
	file << "\t";
	if (visible[k])
	file << 1;
	else file << 0;
	file << "\n";
}

file.close();
		
}

void writeStressStrainToFile(string path1, string path2, std::vector<double>& vMStress, std::vector<double>& plasticStrains, std::vector<double>& elasticStrainsPerNode, std::vector<double>& plasticStrainsPerNode)
{

//create the file stream
ofstream file1(path1.c_str(), ios::out | ios::binary );
ofstream file2(path2.c_str(), ios::out | ios::binary );

double dvalue0,dvalue1,dvalue2,dvalue3;

for (int k = 0; k < vMStress.size(); k++)
{
dvalue0 = vMStress[k];
dvalue1 = plasticStrains[k];

	file1 << dvalue0;	
	file1 << "\t";
	file1 << dvalue1;
	file1 << "\n";
}

file1.close();

for (int k = 0; k < elasticStrainsPerNode.size(); k++)
{
	
	//std::cout << " elastic strains " << elasticStrainsPerNode[k] << std::endl;
	dvalue2 = elasticStrainsPerNode[k];
dvalue3 = plasticStrainsPerNode[k];
	file2 << dvalue2;	
	file2 << "\t";
	file2 << dvalue3;
	file2 << "\n";
}
file2.close();
	
}


int writeMatToFile(const cv::Mat &I, string path) {
 
    //load the matrix size
    int matWidth = I.size().width, matHeight = I.size().height;
 
    //read type from Mat
    int type = I.type();
 
    //declare values to be written
    float fvalue;
    double dvalue;
    Vec3f vfvalue;
    Vec3d vdvalue;
 
    //create the file stream
    ofstream file(path.c_str(), ios::out | ios::binary );
    if (!file)
        return -1;
 
    //write type and size of the matrix first
    file.write((const char*) &type, sizeof(type));
    file.write((const char*) &matWidth, sizeof(matWidth));
    file.write((const char*) &matHeight, sizeof(matHeight));
 
    //write data depending on the image's type
    switch (type)
    {
    default:
        cout << "Error: wrong Mat type: must be CV_32F, CV_64F, CV_32FC3 or CV_64FC3" << endl;
        break;
    // FLOAT ONE CHANNEL
    case CV_32F:
        //cout << "Writing CV_32F image" << endl;
        for (int i=0; i < matWidth*matHeight; ++i) {
            fvalue = I.at<float>(i);
            file.write((const char*) &fvalue, sizeof(fvalue));
        }
        break;
    // DOUBLE ONE CHANNEL
    case CV_64F:
        cout << "Writing CV_64F image" << endl;
        for (int i=0; i < matWidth*matHeight; ++i) {
            dvalue = I.at<double>(i);
            file.write((const char*) &dvalue, sizeof(dvalue));
        }
        break;
 
    // FLOAT THREE CHANNELS
    case CV_32FC3:
        cout << "Writing CV_32FC3 image" << endl;
        for (int i=0; i < matWidth*matHeight; ++i) {
            vfvalue = I.at<Vec3f>(i);
            file.write((const char*) &vfvalue, sizeof(vfvalue));
        }
        break;
 
    // DOUBLE THREE CHANNELS
    case CV_64FC3:
        cout << "Writing CV_64FC3 image" << endl;
        for (int i=0; i < matWidth*matHeight; ++i) {
            vdvalue = I.at<Vec3d>(i);
            file.write((const char*) &vdvalue, sizeof(vdvalue));
        }
        break;
 
    }
 
    //close file
    file.close();
 
    return 0;
}

int writeMatToFile0(const cv::Mat &I, string path) {

    //load the matrix size
    int matWidth = I.size().width, matHeight = I.size().height;

    //read type from Mat
    int type = I.type();

    //declare values to be written
    float fvalue;
    double dvalue;
    Vec3f vfvalue;
    Vec3d vdvalue;

    //create the file stream
    ofstream file(path.c_str(), ios::out | ios::binary );
    if (!file)
        return -1;

    //write type and size of the matrix first
    file.write((const char*) &type, sizeof(type));
    file.write((const char*) &matWidth, sizeof(matWidth));
    file.write((const char*) &matHeight, sizeof(matHeight));

    //write data depending on the image's type
    switch (type)
    {
    default:
        cout << "Error: wrong Mat type: must be CV_32F, CV_64F, CV_32FC3 or CV_64FC3" << endl;
        break;
        // FLOAT ONE CHANNEL
    case CV_32F:
        cout << "Writing CV_32F image" << endl;
        for (int i=0; i < matWidth*matHeight; ++i) {
            fvalue = I.at<float>(i);
            file.write((const char*) &fvalue, sizeof(fvalue));
        }
        break;
        // DOUBLE ONE CHANNEL
    case CV_64F:
        cout << "Writing CV_64F image" << endl;
        for (int i=0; i < matWidth*matHeight; ++i) {
            dvalue = I.at<double>(i);
            file.write((const char*) &dvalue, sizeof(dvalue));
        }
        break;

        // FLOAT THREE CHANNELS
    case CV_32FC3:
        cout << "Writing CV_32FC3 image" << endl;
        for (int i=0; i < matWidth*matHeight; ++i) {
            vfvalue = I.at<Vec3f>(i);
            file.write((const char*) &vfvalue, sizeof(vfvalue));
        }
        break;

        // DOUBLE THREE CHANNELS
    case CV_64FC3:
        cout << "Writing CV_64FC3 image" << endl;
        for (int i=0; i < matWidth*matHeight; ++i) {
            vdvalue = I.at<Vec3d>(i);
            file.write((const char*) &vdvalue, sizeof(vdvalue));
        }
        break;

    }

    //close file
    file.close();

    return 0;
}


template<class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::writeData()
{
	
	//std::string opath = "in/images4/img1%06d.png";
	std::string opath = outputPath.getValue() + "/img1%06d.png";
	std::string opathdepth = outputPath.getValue() + "/depth%06d.png";
	std::string opath1 = outputPath.getValue() + "/depthfile%06d.txt";
	std::string opath2 = outputPath.getValue() + "/depthim%06d.txt";
	std::string opath3 = outputPath.getValue() + "/stressstrain%06d.txt";
	std::string opath4 = outputPath.getValue() + "/stressstrainNode%06d.txt";
	std::string extensionfile = outputPath.getValue() + "/in/images4/extension.txt"; 
        std::vector<Vec3d> pcd1;
	std::vector<double> vmstress;
	std::vector<double> plsstrain;
	std::vector<double> plsstrainnode;
	std::vector<double> elsstrainnode;
	std::vector<bool> visible1;
	
	ofstream extfile(extensionfile.c_str(), ios::out | ios::binary );
	double extension;
	double pos = 0;
	double pos0;


	cv::Mat rtt,rtt1,rtt2,depthi, depthin;	
    for (int frame_count = 0 ;frame_count < listpcd.size(); frame_count++)
    {
		
        std::cout << " ok write " << frame_count << std::endl;
    	pcd1 = *listpcd[frame_count];
		vmstress = *listvm[frame_count];
		elsstrainnode = *listesnode[frame_count];
    	plsstrain = *listps[frame_count];
    	plsstrainnode = *listpsnode[frame_count];

		pos = pcd1[132][1];
		
		if (frame_count == 0) pos0 = pos;
		extension = pos - pos0;
		
		extfile << pos;
		extfile << "\t";
		extfile << extension;
		extfile << "\n";
		
		visible1 = *listvisible[frame_count];
		char buf1[FILENAME_MAX];
        sprintf(buf1, opath1.c_str(), frame_count);
        std::string filename1(buf1);
		
		char buf3[FILENAME_MAX];
        sprintf(buf3, opath3.c_str(), frame_count);
        std::string filename3(buf3);
		
		char buf4[FILENAME_MAX];
        sprintf(buf4, opath4.c_str(), frame_count);
        std::string filename4(buf4);

		writePCDToFile(filename1,pcd1,visible1);
		writeStressStrainToFile(filename3, filename4, vmstress, plsstrain, elsstrainnode, plsstrainnode);

		
		/*delete listimg[frame_count];
		delete listimgseg[frame_count];
		delete listdepth[frame_count];*/
    }
	extfile.close();


    //for (int frame_count = 1 ;frame_count < nimages.getValue()*niterations.getValue(); frame_count++)
    for (int frame_count = 1 ;frame_count < listrtt.size(); frame_count++)
    {	
                std::cout << " ok write rtt0 " << frame_count << std::endl;
                rtt = *listrtt[frame_count];
		cvtColor(rtt,rtt1 ,CV_RGB2BGR);
		cv::flip(rtt1,rtt2,0);
		//cv::flip(rtt1,rtt2,0);
		char buf4[FILENAME_MAX];
                sprintf(buf4, opath.c_str(), frame_count-1);
                std::string filename4(buf4);
                cv::imwrite(filename4,rtt2);
		
                depthi = *listdepth[frame_count];
		char buf5[FILENAME_MAX];
                sprintf(buf5, opathdepth.c_str(), frame_count-1);
                std::string filename5(buf5);
		
		char buf6[FILENAME_MAX];
                sprintf(buf6, opath2.c_str(), frame_count-1);
                std::string filename6(buf6);
                ofstream depthim(buf6, ios::out | ios::binary );

                for (int j = 0; j < depthi.cols; j++)
                for (int i = 0; i< depthi.rows; i++)
		{
		depthim << depthi.at<float>(i,j);
		depthim << "\n";
                }
                cv::flip(depthi,depthin,0);
                depthin.convertTo(depthi, CV_8UC1, 255);

                cv::Mat depth0 = depthi.clone();

                cv::imwrite(filename5,depth0);

                char buf1[FILENAME_MAX];
                sprintf(buf1, opath1.c_str(), frame_count-1);
                std::string filename1(buf1);

                writeMatToFile0(depthin,filename1);
	}
	
	
}

template<class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::setViewPointData()
{
	Eigen::Affine3f scene_sensor_pose (Eigen::Affine3f::Identity ());
	
	source.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
        const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
	
	pcl::PointXYZRGB newPoint;
	for (unsigned int i=0; i<x.size(); i++)
	{
	newPoint.z = x[i][2];
	newPoint.x = x[i][0];
	newPoint.y = x[i][1];
	//cout << "OK " <<  newPoint.z << " " << newPoint.x << " " << newPoint.y  << endl;
	newPoint.r = 0;
	newPoint.g = 0;
	newPoint.b = 0;
	source->points.push_back(newPoint);
	} 
	pcl::PointCloud<pcl::PointXYZRGB>& point_cloud = *source;
	
	scene_sensor_pose = Eigen::Affine3f (Eigen::Translation3f (point_cloud.sensor_origin_[0],
                                                             point_cloud.sensor_origin_[1],
                                                             point_cloud.sensor_origin_[2])) * 
															 Eigen::Affine3f (point_cloud.sensor_orientation_);
															 
  Eigen::Affine3f viewer_pose = scene_sensor_pose;
  Eigen::Vector3f pos_vector = viewer_pose * Eigen::Vector3f(0, 0, 0);  
  Eigen::Vector3f look_at_vector = viewer_pose.rotation () * Eigen::Vector3f(0, 0, 1) + pos_vector;
  Eigen::Vector3f up_vector = viewer_pose.rotation () * Eigen::Vector3f(0, -1, 0);
  /*viewer.setCameraPosition (pos_vector[0], pos_vector[1], pos_vector[2],
                            look_at_vector[0], look_at_vector[1], look_at_vector[2],
                            up_vector[0], up_vector[1], up_vector[2]);*/
							
	int hght = 480;
	int wdth = 640;
			
    sofa::gui::GUIManager::SetDimension(wdth,hght);
			
	sofa::gui::BaseGUI *gui = sofa::gui::GUIManager::getGUI();
    if (!gui)
    {
        std::cout << " no gui " << std::endl; 
    }
    sofa::gui::BaseViewer * viewer = gui->getViewer();
    if (!viewer)
    {
        std::cout << " no viewer " << std::endl; 

    }
    //viewer->getView(pos,orient);
	
			{	
			Vec3d position;
            Quat orientation;
			
			orientation[0] = point_cloud.sensor_orientation_.w ();
			orientation[1] = point_cloud.sensor_orientation_.x ();
			orientation[2] = point_cloud.sensor_orientation_.y ();
			orientation[3] = point_cloud.sensor_orientation_.z ();
			
			position[0] = point_cloud.sensor_origin_[0];
			position[1] = point_cloud.sensor_origin_[1];
			position[2] = point_cloud.sensor_origin_[2];

    //if (currentCamera)
			//currentCamera->setView(position, orientation);
			viewer->setView(position, orientation);

		    //glutPostRedisplay();
			        //glPushAttrib( GL_LIGHTING_BIT | GL_ENABLE_BIT | GL_LINE_BIT | GL_CURRENT_BIT);
        //glPopAttrib();
			}
	
}

template <class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::generateData(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
	
	int t = (int)this->getContext()->getTime();
	 //getImages();
	 //readData()
	 //timeOverall
	 double timeT;

	if (t == 0)
	{
		timeOverall = 0;
		timeTotal = 0;
		timeT = (double)getTickCount();
		//extractTargetContour();
	}
	else if (t > 0 && t%niterations.getValue() == 0)
	{
	double time0 = (double)getTickCount();
	double time1 = (double)getTickCount();
	}
	
    if(ks.getValue()==0) return;

    VecDeriv&        f = *_f.beginEdit();           //WDataRefVecDeriv f(_f);
    const VecCoord&  x = _x.getValue();			//RDataRefVecCoord x(_x);
    const VecDeriv&  v = _v.getValue();			//RDataRefVecDeriv v(_v);

	ReadAccessor< Data< VecCoord > > ssn(sourceSurfaceNormalsM);
	
	ReadAccessor< Data<helper::vector<Real> > > vM(vonMisesStress);
	ReadAccessor< Data<helper::vector<Real> > > eS(elasticStrainsN);
	ReadAccessor< Data<helper::vector<Real> > > pS(plasticStrainsN);
	ReadAccessor< Data<helper::vector<Real> > > tS(totalStrainsN);
	ReadAccessor< Data<helper::vector<Real> > > eSNode(elasticStrainsPerNode);
	ReadAccessor< Data<helper::vector<Real> > > pSNode(plasticStrainsPerNode);
	ReadAccessor< Data<helper::vector<Real> > > tSNode(totalStrainsPerNode);

    const vector<Spring>& s = this->springs.getValue();
    this->dfdx.resize(s.size());
    this->closestPos.resize(s.size());

	initSourceSurface();
		if (t < 2)
	setViewPointData();

                sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
                sofa::component::visualmodel::BaseCamera::SPtr currentCamera;// = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();
                root->get(currentCamera);
        if (t > 1)
        {
	double znear = currentCamera->getZNear();
	double zfar = currentCamera->getZFar();
        //meshprocessing->getSourceVisible(znear, zfar);
	sourceVisible = meshprocessing->sourceVisible;
	sourceVisiblePositions.setValue(meshprocessing->getSourceVisiblePositions());
	indicesVisible = meshprocessing->indicesVisible;
	depthMap = meshprocessing->depthMap.clone();	
			std::cout << " npoints " << std::endl;

	
			double maxz = 0;
	int kmaxz;
	
        double minx = 100;
	int kminx;
	
	double minx1 = 100;
	int kminx1;
	
	double maxx = 0;
	int kmaxx;
	
		double maxx1 = 0;
	int kmaxx1;
	
	double maxz2 = 0;
	int kmaxz2;
		double maxz3 = 0;
	int kmaxz3;
			double maxz4 = 0;
	int kmaxz4;
			double maxz5 = 0;
	int kmaxz5;
			double maxz6 = 0;
	int kmaxz6;
				double maxz7 = -1000;
	int kmaxz7;
			double maxz8 =  -1000;
	int kmaxz8;
			double maxz9 =  -1000;
	int kmaxz9;
			double maxz10 =  -1000;
	int kmaxz10;
	
	double minz = 100;
	int kminz;
		double minz2 = 100;
	int kminz2;
		double minz3 = 100;
	int kminz3;

	if (t > 0 && t%niterations.getValue() == 0)
	{
	double timertt = (double)getTickCount();
			std::cout << " npoints 12 " << std::endl;

	rtt = new cv::Mat;
	depthl = new cv::Mat;
	cv::Mat rtt_,rttdepth_;
        rendertexturear->renderToTextureDepth(rtt_,rttdepth_);
        *rtt = rtt_.clone();
        *depthl = rttdepth_.clone();

	listrtt.push_back(rtt);
	listdepth.push_back(depthl);
	timertt = ((double)getTickCount() - timertt)/getTickFrequency();
    cout << "Time RTT " << timertt << endl;
	
	pcd = new std::vector<Vec3d>;
	pcd->resize(0);
	
	Vec3d pcd0;
	
	visible = new std::vector<bool>;
	visible->resize(0);
	bool vis0;
	
	vm = new std::vector<double>;
	vm->resize(0);
	
	es = new std::vector<double>;
	es->resize(0);
	ps = new std::vector<double>;
	ps->resize(0);
	ts = new std::vector<double>;
	ts->resize(0);
	
	esnode = new std::vector<double>;
	esnode->resize(0);
	
	psnode = new std::vector<double>;
	psnode->resize(0);
	
	tsnode = new std::vector<double>;
	tsnode->resize(0);
	std::cout << " srouce visible" << sourceVisible.size() << std::endl;
	
	for (unsigned int i=0; i<x.size(); i++)
	{		
		if(sourceVisible[i]) 
		{
			vis0 = true;
			

		} 
		else vis0 = false;
		
	if (x[i][0] < minx)
	{
		minx = x[i][0];
		kminx = i;
	}
	else
	{
	if (x[i][0] < minx1)
	{
		minx1 = x[i][0];
		kminx1 = i;
	}
	}
	
    if (x[i][0] > maxx)
	{
		maxx = x[i][0];
		kmaxx = i;
	}
	else
	{
	if (x[i][0] > maxx1)
	{
		maxx1 = x[i][0];
		kmaxx1 = i;
	}
	}
	
    if (x[i][1] > maxz)
	{
		maxz = x[i][1];
		kmaxz = i;
	}
	else
	{
	
	if (x[i][1] > maxz2)
	{
		maxz2 = x[i][1];
		kmaxz2 = i;
	}
	else
	{
	
	if (x[i][1] > maxz3)
	{
		maxz3 = x[i][1];
		kmaxz3 = i;
	}
	else
	{
    if (x[i][1] > maxz4)
	{
		maxz4 = x[i][1];
		kmaxz4 = i;
	}
	else
	{
	if (x[i][1] > maxz5)
	{
		maxz5 = x[i][1];
		kmaxz5 = i;
	}
	else
	{
	if (x[i][1] > maxz6)
	{
		maxz6 = x[i][1];
		kmaxz6 = i;
	}
	else
	{
		
		//std::cout << " ok " << x[i][1]  << std::endl;
	if (x[i][1] > maxz7)
	{
		maxz7 = x[i][1];
		kmaxz7 = i;
	}
	else{
	if (x[i][1] > maxz8)
	{
		maxz8 = x[i][1];
		kmaxz8 = i;
	}
	else
	{
			if (x[i][1] > maxz9)
	{
		maxz9 = x[i][1];
		kmaxz9 = i;
	}
	}
	}
	}

	}
	}
	}
		
	}
		
	}
	
	if (x[i][1] < minz)
	{
		minz = x[i][1];
		kminz = i;
	}
	else
	{
	if (x[i][1] < minz2)
	{
		minz2 = x[i][1];
		kminz2 = i;
	}
	else
	{
	if (x[i][1] < minz3)
	{
		minz3 = x[i][1];
		kminz3 = i;
	}
	}
	}
		
		
		visible->push_back(vis0);
			
		{
	pcd0[0] = (double)x[i][0];
	pcd0[1] = (double)x[i][1];
	pcd0[2] = (double)x[i][2];
	pcd->push_back(pcd0);
		}
	}
	listpcd.push_back(pcd);
	listvisible.push_back(visible);
	
	std::cout << " " << eS.size() << " " << pS.size() << " " << tS.size() << std::endl;
	
	for (unsigned int i=0; i<vM.size(); i++)
		{
			vm->push_back((double)vM[i]);
			es->push_back((double)eS[i]);
			ps->push_back((double)pS[i]);
			ts->push_back((double)tS[i]);

		}
	listvm.push_back(vm);
	listps.push_back(ps);
	listes.push_back(es);
	listts.push_back(ts);


	
	for (unsigned int i=0; i<pSNode.size(); i++)
		{
			esnode->push_back((double)eSNode[i]);
			psnode->push_back((double)pSNode[i]);
			tsnode->push_back((double)tSNode[i]);

		}
		
	listesnode.push_back(esnode);
	listpsnode.push_back(psnode);
	listpsnode.push_back(tsnode);

		}
			
	//std::cout << " time " << t << " " << listpcd.size() << std::endl;
	
	//if (t == nimages.getValue()*niterations.getValue()-4*niterations.getValue() - 1)
	if (t == 400)
	{
		writeData();
		getchar();
	}
			std::cout << " npoints 11 " << std::endl;

	/*for (int kk = 0; kk < x.size(); kk++)
	{
		if ((x[kk][0]-x[kmaxz][0])*(x[kk][0]-x[kmaxz][0]) + (x[kk][1]-x[kmaxz][1])*(x[kk][1]-x[kmaxz][1]) + (x[kk][2]-x[kmaxz][2])*(x[kk][2]-x[kmaxz][2]) < 0.0004)
		{
		std::cout << " kk " << kk <<std::endl;	
		}
			
	}*/
	
	std::cout << " kminx " << kminx << " kmax x " << kmaxx << " kmax z " << kmaxz << " kmaxz2 " << kmaxz2 << " kmaxz3 " << kmaxz3 << " kminz " << kminz << " kminz2 " << kminz2 << " kminz3 " << kminz3 <<  std::endl;
	std::cout << " kmaxz4 " << kmaxz4 << " kmaxz5 " << kmaxz5 << " kmaxz6 " << kmaxz6 << " kmaxz7 " << kmaxz7 << " kmaxz8 " << kmaxz8 << " kmaxz9 " << kmaxz9 << std::endl;
	std::cout << " kminx1 " << kminx1 << " kmax x " << kmaxx1 <<  std::endl;

	int key=cvWaitKey(3);
	if ((char)key == 27) return;
	
	//computeTargetNormals();
	
	double time = (double)getTickCount();

	m_potentialEnergy = 0;
	
	//Bunny 
	/*int ind0 = 145;
	int ind1 = 132;
	int ind2 = 139;
	int ind3 = 126;
	int ind4 = 0;*/
	//res max
	/*int ind0 = 877;
	int ind1 = 148;
	int ind2 = 129;
	int ind3 = 810;
	int ind4 = 847;
	int ind5 = 10;
    int ind6 = 87;
	int ind7 = 147;
	int ind8 = 148;
	int ind9 = 153;
	int ind10 = 175;
	int ind11 = 487;
	int ind12 = 583;*/
		
		// resmin
	/*int ind0 = 201;
	int ind1 = 157;
	int ind2 = 186;
	int ind3 = 180;
	int ind4 = 175;
	int ind5 = 172;
    int ind6 = 183;
	int ind7 = 346;
	int ind8 = 232;
	int ind9 = 256;
	int ind10 = 247;
	int ind11 = 181;
	int ind12 = 262;
	*/
	
	//resmed
	int ind0 = 253;
	int ind1 = 157;
	int ind2 = 186;
	//int ind3 = 180;
	int ind3 = 305;
	int ind4 = 308;
	int ind5 = 306;
    int ind6 = 307;
	int ind7 = 323;
	int ind8 = 305;
	int ind9 = 294;
	int ind10 = 399;
	
	
	int ind11 = 252;
	int ind12 = 185;
	
	/* kminx 253 kmax x 186 kmax z 157 kmaxz2 180 kmaxz3 308 kminz 407 kminz2 408 kminz3 477
 kmaxz4 306 kmaxz5 307 kmaxz6 323 kmaxz7 305 kmaxz8 294 kmaxz9 399*/

	
	
		/*int ind0 = 201;
	int ind1 = 157;
	int ind2 = 186;
	int ind3 = 180;
	int ind4 = 181;
	int ind5 = 172;
    int ind6 = 183;
	int ind7 = 346;
	int ind8 = 262;
	int ind9 = 256;
	int ind10 = 240;*/
	/*int ind11 = 487;
	int ind12 = 583;*/
 /*kk 147
 kk 148
 kk 153
 kk 175
 kk 487
 kk 583
 kk 600
 kk 706
 kk 786
 kk 810
 kk 847
 kk 850
 kk 855
 kk 871*/
	
	/*int ind0 = 335;
	int ind1 = 15;
	int ind3 = 79;
	int ind2 = 399;
	//int ind3 = 399;
	int ind4 = 0;*/

	Vector3 trans0, trans1, trans2, trans3, trans4, trans5, trans6, trans7, trans8, trans9, trans10, trans11, trans12;
	
	
	if (t == 0)
		{
	/*trans0[0] = x[ind0][0] - 0.10;
	trans0[1] = x[ind0][1] - 0.0;
	trans0[2] = x[ind0][2] - 0.0;
	trans1[0] = x[ind1][0];
	trans1[1] = x[ind1][1] + 0.40;
	trans1[2] = x[ind1][2] - 0.0;
	trans2[0] = x[ind2][0] + 0.10;
	trans2[1] = x[ind2][1] ;
	trans2[2] = x[ind2][2] ;
	trans3[0] = x[ind3][0] - 0.08;
	trans3[1] = x[ind3][1] ;
	trans3[2] = x[ind3][2] ;
	trans4[0] = x[ind4][0] ;
	trans4[1] = x[ind4][1] ;
	trans4[2] = x[ind4][2] -0.1;*/
	
	
	//Bunny
		/*trans0[0] = x[ind0][0];
	trans0[1] = x[ind0][1] - 0.11;
	trans0[2] = x[ind0][2] - 0.0;
	trans1[0] = x[ind1][0];
	trans1[1] = x[ind1][1] + 0.20;
	trans1[2] = x[ind1][2] - 0.0;
	trans2[0] = x[ind2][0] + 0.08;
	trans2[1] = x[ind2][1] ;
	trans2[2] = x[ind2][2] ;
	trans3[0] = x[ind3][0] - 0.08;
	trans3[1] = x[ind3][1] ;
	trans3[2] = x[ind3][2] ;
	trans4[0] = x[ind4][0] ;
	trans4[1] = x[ind4][1] ;
	trans4[2] = x[ind4][2] -0.1;*/
	
	//pizza 1
	
	trans0[0] = x[ind0][0]+0.12;
	trans0[1] = x[ind0][1] - 0.1;
	trans0[2] = x[ind0][2] - 0.1;
	
	trans1[0] = x[ind1][0] + 0.0;
	trans1[1] = x[ind1][1] + 0.17;
	trans1[2] = x[ind1][2] + 0.0;
	
	trans2[0] = x[ind2][0] - 0.12;
	trans2[1] = x[ind2][1] -0.1;
	trans2[2] = x[ind2][2] -0.1;
	trans3[0] = x[ind3][0] - 0.0;
	trans3[1] = x[ind3][1] +0.17;
	trans3[2] = x[ind3][2] -0.0;
	
	trans4[0] = x[ind4][0] ;
	trans4[1] = x[ind4][1] +0.17;
	trans4[2] = x[ind4][2] +0.0;

	trans4[0] = x[ind4][0] ;
	trans4[1] = x[ind4][1] +0.17;
	trans4[2] = x[ind4][2] +0.0;
	
	trans4[0] = x[ind4][0] ;
	trans4[1] = x[ind4][1] +0.17;
	trans4[2] = x[ind4][2] +0.0;
	
	trans5[0] = x[ind5][0] ;
	trans5[1] = x[ind5][1] +0.17;
	trans5[2] = x[ind5][2] +0.0;

	trans6[0] = x[ind6][0] ;
	trans6[1] = x[ind6][1] +0.17;
	trans6[2] = x[ind6][2] +0.0;
	
	trans7[0] = x[ind7][0] ;
	trans7[1] = x[ind7][1] +0.17;
	trans7[2] = x[ind7][2] +0.0;
	
	trans8[0] = x[ind8][0] ;
	trans8[1] = x[ind8][1] +0.17;
	trans8[2] = x[ind8][2] +0.0;
	
		
	trans9[0] = x[ind9][0] ;
	trans9[1] = x[ind9][1] +0.17;
	trans9[2] = x[ind9][2] +0.0;
	
	trans10[0] = x[ind10][0] ;
	trans10[1] = x[ind10][1] +0.17;
	trans10[2] = x[ind10][2] +0.0;

	
	trans11[0] = x[ind11][0]+0.12;
	trans11[1] = x[ind11][1] - 0.1;
	trans11[2] = x[ind11][2] - 0.1;
	
	trans12[0] = x[ind12][0] - 0.12;
	trans12[1] = x[ind12][1] -0.1;
	trans12[2] = x[ind12][2] -0.1;
	
	closestPos[ind1]=trans1;
	closestPos[ind0]=trans0;
	closestPos[ind2]=trans2;
	closestPos[ind3]=trans3;
	closestPos[ind4]=trans4;
	closestPos[ind5]=trans5;
	closestPos[ind6]=trans6;
	closestPos[ind7]=trans7;
	closestPos[ind8]=trans8;
	closestPos[ind9]=trans9;
	closestPos[ind10]=trans10;
	closestPos[ind11]=trans11;
	closestPos[ind12]=trans12;
        }

        /*if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind0, s[ind0]);
	if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind1, s[ind1]);
	if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind2, s[ind2]);
	if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind3, s[ind3]);
	if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind4, s[ind4]);
	if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind5, s[ind5]);
	if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind6, s[ind6]);
	if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind7, s[ind7]);
	if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind8, s[ind8]);
        if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind9, s[ind9]);
        if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind10, s[ind10]);*/
	//if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind11, s[ind11]);
	//if (t > 50 && t < 300) this->addSpringForce(m_potentialEnergy,f,x,v, ind12, s[ind12]);

		/*else{
		if(!targetBorder[(int)closestSource[i].begin()->second])
        this->addSpringForce(m_potentialEnergy,f,x,v, i, s[i]);
		else this->addSpringForceContour(m_potentialEnergy,f,x,v, i, s[i]);}*/		
	
	
	//std::cout << " Error " << error << std::endl;
    _f.endEdit();
	
	time = ((double)getTickCount() - time)/getTickFrequency();
    cout << "Time addforce " << time << endl;
	timeOverall += time;
	
	timeii = (double)getTickCount();
	timeTotal =((double)getTickCount() - timeT)/getTickFrequency();
        cout << "Time total " << timeTotal << endl;
        }
		
}

template <class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::addForce(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
		generateData(mparams, _f, _x, _v);
}


template <class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
{
    int a = spring.m1;
    Coord u = this->closestPos[i]-p[a];
    Real d = u.norm();
    if( d>1.0e-4 )
    {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        potentialEnergy += elongation * elongation * spring.ks / 2;
        /*          serr<<"addSpringForce, p = "<<p<<sendl;
        serr<<"addSpringForce, new potential energy = "<<potentialEnergy<<sendl;*/
        Deriv relativeVelocity = -v[a];
        Real elongationVelocity = dot(u,relativeVelocity);
        Real forceIntensity;
        if(theCloserTheStiffer.getValue())
        {
            Real ks_max=spring.ks;
            Real ks_min=spring.ks/10;
            Real ks_mod = ks_min*(max-elongation)/(max-min)+ks_max*(elongation-min)/(max-min);
            forceIntensity = (Real)(ks_mod*elongation+spring.kd*elongationVelocity);
        }
        else {
			if (elongation < 0.02)
		forceIntensity = (Real)(spring.ks*elongation+spring.kd*elongationVelocity);
		else forceIntensity = (Real)(spring.ks*elongation+spring.kd*elongationVelocity);
		}
        Deriv force = u*forceIntensity;
        f[a]+=force;
        Mat& m = this->dfdx[i];
        Real tgt = forceIntensity * inverseLength;
        for( int j=0; j<N; ++j )
        {
            // anisotropic
            //for( int k=0; k<N; ++k ) m[j][k] = tgt * u[j] * u[k];

            // isotropic
            for( int k=0; k<N; ++k ) m[j][k] = ((Real)spring.ks-tgt) * u[j] * u[k];
            m[j][j] += tgt;
        }
			//dfdx1[i] = m;

    }
    else // null length, no force and no stiffness
    {
        Mat& m = this->dfdx[i];
        for( int j=0; j<N; ++j )
        {
            for( int k=0; k<N; ++k )
            {
                m[j][k] = 0;
            }
        }
    }
}

template<class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::addSpringDForce(VecDeriv& df,const  VecDeriv& dx, int i, const Spring& spring, double kFactor, double /*bFactor*/)
{
    const int a = spring.m1;
    const Coord d = -dx[a];
    Deriv dforce = this->dfdx[i]*d;
    dforce *= kFactor;
    df[a]+=dforce;
    //serr<<"addSpringDForce, a="<<a<<", b="<<b<<", dforce ="<<dforce<<sendl;
}

template <class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::addDForce(const core::MechanicalParams* mparams,DataVecDeriv& _df , const DataVecDeriv&  _dx )
{

    VecDeriv& df = *_df.beginEdit();		//WDataRefVecDeriv df(_df);
    const VecDeriv&  dx = _dx.getValue();	// RDataRefVecDeriv dx(_dx);

    double kFactor 		 =  mparams->kFactor();
    double bFactor       =  mparams->bFactor();

    if(ks.getValue()==0) return;

    const vector<Spring>& s = this->springs.getValue();

    //serr<<"addDForce, dx = "<<dx<<sendl;
    //serr<<"addDForce, df before = "<<f<<sendl;
    for (unsigned int i=0; i<s.size(); i++)
    {
        this->addSpringDForce(df,dx, i, s[i], kFactor, bFactor);
    }
    //serr<<"addDForce, df = "<<f<<sendl;
	
	
    _df.endEdit();

}

template<class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *m, SReal kFactor, unsigned int &offset)
{
	    if(ks.getValue()==0) return;

    double kFact = kFactor;
	
    const vector<Spring >& ss = this->springs.getValue();
    const unsigned int n = ss.size() < this->dfdx.size() ? ss.size() : this->dfdx.size();
    for (unsigned int e=0; e<n; e++)
    {
        const Spring& s = ss[e];
        unsigned p1 = offset+Deriv::total_size*s.m1;
        const Mat& mt = this->dfdx[e];
        for(int i=0; i<N; i++)
            for (int j=0; j<N; j++)
            {
                Real k = (Real)(mt[i][j]*kFact);
                m->add(p1+i,p1+j, -k);
            }
    }
}

            
template<class DataTypes, class DepthTypes>
void DataGeneration<DataTypes, DepthTypes>::draw(const core::visual::VisualParams* vparams)
{	
	
	int t = (int)this->getContext()->getTime();
		
	double timef = 0;
	timef = (double)getTickCount();
	//if (t > 0 && t%niterations.getValue() == niterations.getValue()-1)
	timeOverall += (timef - timei)/getTickFrequency();
	timeResolution += (timef - timei)/getTickFrequency();
	{
    cout <<" t " << t << " Time overall draw " << timei << " " << timef << " " << (timef - timei)/getTickFrequency() << endl;
	}
	
	ReadAccessor< Data< VecCoord > > x(*this->getMState()->read(core::ConstVecCoordId::position()));
		
        const VecCoord&  p0 = this->mstate->read(core::ConstVecCoordId::position())->getValue();

	sourcePositions.setValue(p0);
	
	if (t > 0 && t%niterations.getValue() == niterations.getValue()-1)
	{
	std::cout << " t " << t << " " << t/niterations.getValue() <<  " time overall " << timeOverall << std::endl;
	timeFile << t/niterations.getValue();	
	timeFile << "\t";
	timeFile << timeOverall;
	timeFile << "\t";
	timeFile << timeSeg;
	timeFile << "\t";
	timeFile << timeRigid;
	timeFile << "\t";
	timeFile << timeSourceContour;
	timeFile << "\t";
	timeFile << timeAddforce;
	timeFile << "\t";
	timeFile << timeResolution;
	timeFile << "\n";
	

	}
	
	  glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT | GL_ENABLE_BIT);
	  glDisable( GL_LIGHTING);

    {

		glEnable(GL_BLEND);
        glPointSize( 10);
        glBegin( GL_POINTS);
//std::cout << " size target " << targetPositions.getValue().size() << " dist size " << dists.size() << std::endl;

	 Vector3 coefs;
	 int index, id;
	 float xp0, yp0;
					//}
		  int kcp = 0;
for (unsigned int i=0; i<x.size(); i++)
            {
					if (t > 5)
					{
					if(sourceVisible[i])
						{
						
						if (useContour.getValue())
					if (sourceBorder[i])
						{
							//std::cout << " displ " << x[i][2] << std::endl;

                //sofa::helper::gl::Color::setHSVA(dists[i]*240./max,1.,.8,1.);
			   sofa::helper::gl::Color::setHSVA(140,1.,.8,1.5);
				//if (kcp== 0){
               glVertex3d(x[i][0],x[i][1],x[i][2]);
			   //kcp++;
					}
				}
				}
            }
			


		
		glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT | GL_ENABLE_BIT);
        glDisable( GL_LIGHTING);
		glBegin( GL_POINTS);
		
        glEnd();
        glPointSize( 1);

        glPopAttrib();
    }

}

}
}
} // namespace sofa


