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
#define SOFA_RGBDTRACKING_GENERATESYNTHETICDATA_CPP

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/Mapping.inl>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/gui/BaseGUI.h>
#include <sofa/gui/BaseViewer.h>
#include <sofa/gui/GUIManager.h>
#include <iostream>
#include <map>

#ifdef USING_OMP_PRAGMAS
    #include <omp.h>
#endif


#include <sofa/component/loader/MeshObjLoader.h>
#include <sofa/component/engine/NormalsFromPoints.h>
#include <limits>
#include <set>
#include <iterator>
#include <sofa/helper/gl/Color.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/component/visualmodel/InteractiveCamera.h>

#ifdef Success
  #undef Success
#endif

#include <pcl-1.7/pcl/common/common_headers.h>
#include <pcl-1.7/pcl/point_cloud.h>
#include <pcl-1.7/pcl/filters/extract_indices.h>
#include <pcl-1.7/pcl/point_types.h>
#include <pcl-1.7/pcl/impl/point_types.hpp>
#include <pcl-1.7/pcl/features/normal_3d.h>
#include <pcl-1.7/pcl/PointIndices.h>
#include <pcl-1.7/pcl/PolygonMesh.h>
#include <pcl-1.7/pcl/visualization/common/actor_map.h>
#include <pcl-1.7/pcl/visualization/common/common.h>
#include <pcl-1.7/pcl/visualization/point_cloud_geometry_handlers.h>
#include <pcl-1.7/pcl/visualization/point_cloud_color_handlers.h>
#include <pcl-1.7/pcl/visualization/point_picking_event.h>
#include <pcl-1.7/pcl/visualization/area_picking_event.h>
#include <pcl-1.7/pcl/visualization/interactor_style.h>
#include <pcl-1.7/pcl/visualization/pcl_visualizer.h>
#include <pcl-1.7/pcl/visualization/keyboard_event.h>
#include <pcl-1.7/pcl/correspondence.h>
#include <pcl-1.7/pcl/features/normal_3d_omp.h>
#include <pcl-1.7/pcl/features/shot_omp.h>
#include <pcl-1.7/pcl/features/board.h>
#include <pcl-1.7/pcl/console/parse.h>
#include <pcl-1.7/pcl/common/projection_matrix.h>
#include <pcl-1.7/pcl/common/pca.h>
#include <pcl-1.7/pcl/surface/reconstruction.h>
#include <pcl-1.7/pcl/io/pcd_io.h>
#include <pcl-1.7/pcl/common/io.h>
#include <pcl-1.7/pcl/registration/transforms.h>
#include <pcl-1.7/pcl/keypoints/sift_keypoint.h>
#include <pcl-1.7/pcl/keypoints/harris_3d.h>
#include <pcl-1.7/pcl/ModelCoefficients.h>
#include <pcl-1.7/pcl/sample_consensus/method_types.h>
#include <pcl-1.7/pcl/sample_consensus/model_types.h>
#include <pcl-1.7/pcl/segmentation/sac_segmentation.h>
#include <pcl-1.7/pcl/search/kdtree.h>
#include <pcl-1.7/pcl/segmentation/extract_clusters.h>
#include <pcl-1.7/pcl/features/fpfh_omp.h>
#include <pcl-1.7/pcl/features/feature.h>
#include <pcl-1.7/pcl/features/pfh.h>
#include <pcl-1.7/pcl/features/pfhrgb.h>
#include <pcl-1.7/pcl/features/3dsc.h>
#include <pcl-1.7/pcl/features/shot_omp.h>
#include <pcl-1.7/pcl/filters/filter.h>
#include <pcl-1.7/pcl/registration/transformation_estimation_svd.h>
#include <pcl-1.7/pcl/registration/icp.h>
#include <pcl-1.7/pcl/registration/correspondence_rejection_sample_consensus.h>
#include <pcl-1.7/pcl/search/kdtree.h>
#include <pcl-1.7/pcl/common/transforms.h>
#include <pcl-1.7/pcl/range_image/range_image.h>


#include "GenerateSyntheticData.h"

using std::cerr;
using std::endl;

namespace sofa
{

namespace component
{

namespace forcefield
{

    using namespace sofa::defaulttype;

      SOFA_DECL_CLASS(GenerateSyntheticData)

      // Register in the Factory
      int GenerateSyntheticDataClass = core::RegisterObject("Compute forces based on closest points from/to a target surface/point set")
    #ifndef SOFA_FLOAT
        .add< GenerateSyntheticData<Vec3dTypes,ImageUC> >()
        .add< GenerateSyntheticData<Vec3dTypes,ImageUS> >()
        .add< GenerateSyntheticData<Vec3dTypes,ImageF> >()
    #endif
    #ifndef SOFA_DOUBLE
        .add< GenerateSyntheticData<Vec3fTypes,ImageUC> >()
        .add< GenerateSyntheticData<Vec3fTypes,ImageUS> >()
        .add< GenerateSyntheticData<Vec3fTypes,ImageF> >()
    #endif
    ;

    #ifndef SOFA_FLOAT
      template class SOFA_RGBDTRACKING_API GenerateSyntheticData<Vec3dTypes,ImageUC>;
      template class SOFA_RGBDTRACKING_API GenerateSyntheticData<Vec3dTypes,ImageUS>;
      template class SOFA_RGBDTRACKING_API GenerateSyntheticData<Vec3dTypes,ImageF>;
    #endif
    #ifndef SOFA_DOUBLE
      template class SOFA_RGBDTRACKING_API GenerateSyntheticData<Vec3fTypes,ImageUC>;
      template class SOFA_RGBDTRACKING_API GenerateSyntheticData<Vec3fTypes,ImageUS>;
      template class SOFA_RGBDTRACKING_API GenerateSyntheticData<Vec3fTypes,ImageF>;
    #endif

using namespace helper;

template <class DataTypes, class DepthTypes>
GenerateSyntheticData<DataTypes, DepthTypes>::GenerateSyntheticData(core::behavior::MechanicalState<DataTypes> *mm )
    : Inherit(mm)
    , depthImage(initData(&depthImage,DepthTypes(),"depthImage","depth map"))
    //, depthTransform(initData(&depthTransform, TransformType(), "depthTransform" , ""))
    , image(initData(&image,ImageTypes(),"image","image"))
    //, transform(initData(&transform, TransformType(), "transform" , ""))
    , ks(initData(&ks,(Real)0.0,"stiffness","uniform stiffness for the all springs."))
    , kd(initData(&kd,(Real)0.0,"damping","uniform damping for the all springs."))
    , cacheSize(initData(&cacheSize,(unsigned int)10,"cacheSize","number of closest points used in the cache to speed up closest point computation."))
    , blendingFactor(initData(&blendingFactor,(Real)1,"blendingFactor","blending between projection (=0) and attraction (=1) forces."))
    , outlierThreshold(initData(&outlierThreshold,(Real)7,"outlierThreshold","suppress outliers when distance > (meandistance + threshold*stddev)."))
    , normalThreshold(initData(&normalThreshold,(Real)0,"normalThreshold","suppress outliers when normal.closestPointNormal < threshold."))
    , projectToPlane(initData(&projectToPlane,false,"projectToPlane","project closest points in the plane defined by the normal."))
    , rejectBorders(initData(&rejectBorders,false,"rejectBorders","ignore border vertices."))
    , rejectOutsideBbox(initData(&rejectOutsideBbox,false,"rejectOutsideBbox","ignore source points outside bounding box of target points."))
    , springs(initData(&springs,"spring","index, stiffness, damping"))
	, sourceSurfacePositions(initData(&sourceSurfacePositions,"sourceSurface","Points of the surface of the source mesh."))
    , sourceTriangles(initData(&sourceTriangles,"sourceTriangles","Triangles of the source mesh."))
    , sourceNormals(initData(&sourceNormals,"sourceNormals","Normals of the source mesh."))
	, sourceSurfaceNormals(initData(&sourceSurfaceNormals,"sourceSurfaceNormals","Normals of the surface of the source mesh."))
    , showArrowSize(initData(&showArrowSize,0.01f,"showArrowSize","size of the axis."))
    , drawMode(initData(&drawMode,0,"drawMode","The way springs will be drawn:\n- 0: Line\n- 1:Cylinder\n- 2: Arrow."))
    , drawColorMap(initData(&drawColorMap,true,"drawColorMap","Hue mapping of distances to closest point"))
	, useVisible(initData(&useVisible,true,"useVisible","Use the vertices of the visible surface of the source mesh"))
	, useGroundTruth(initData(&useGroundTruth,true,"useGroundTruth","Use the vertices of the visible surface of the source mesh"))	, generateSynthData(initData(generateSynthData&, false,"generateSynthData","Generate synthetic data"))
	, niterations(initData(niterations&,3,"nIterations","Number of iterations in the tracking process"))
	, nimages(initData(nimages&,400,"nImages","Number of images to read"))
{
	this->addAlias(&depthImage, "depthImage");
	depthImage.setGroup("depthImage");
	depthImage.setReadOnly(true); 
	this->addAlias(&image, "image");
	image.setGroup("image");
	image.setReadOnly(true); 
	pcl = false;
	disp = false;
	iter_im = 0;
	timeTotal = 0;
	color = cv::Mat::zeros(240,320, CV_8UC3); 
    listimg.resize(0);
    listimgseg.resize(0);
    listdepth.resize(0);
	listrtt.resize(0);
	listpcd.resize(0);
	listvisible.resize(0);
	
	paramSoft.sigma_p2 = 0.005f; // default parameters
    paramSoft.sigma_inf = 0.00001f;
    paramSoft.sigma_factor = 0.9f;
    paramSoft.d_02 = 0.01f;   

    //
    // initialize paramters
    //
	sigma_p2 = paramSoft.sigma_p2;
    sigma_inf = paramSoft.sigma_inf;
    sigma_factor = paramSoft.sigma_factor;
    d_02 = paramSoft.d_02;  
}

template <class DataTypes, class DepthTypes>
GenerateSyntheticData<DataTypes, DepthTypes>::~GenerateSyntheticData()
{
}

template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::initCamera()
{
    //softk.init();
}


template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::reinit()
{
    for (unsigned int i=0;i<springs.getValue().size();++i)
    {
        (*springs.beginEdit())[i].ks = (Real) ks.getValue();
        (*springs.beginEdit())[i].kd = (Real) kd.getValue();
    }
}

template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::init()
{

    this->Inherit::init();
    core::objectmodel::BaseContext* context = this->getContext();

    if(!(this->mstate)) this->mstate = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());

	    // Get source triangles
    if(!sourceTriangles.getValue().size()) {
        sofa::component::loader::MeshObjLoader *meshobjLoader;
        this->getContext()->get( meshobjLoader, core::objectmodel::BaseContext::Local);
        if (meshobjLoader) {sourceTriangles.virtualSetLink(meshobjLoader->triangles); sout<<"imported triangles from "<<meshobjLoader->getName()<<sendl;
		}
    }
    // Get source normals
    if(!sourceNormals.getValue().size()) serr<<"normals of the source model not found"<<sendl;

    // add a spring for every input point
    const VecCoord& x = *this->mstate->getX(); 			//RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
	this->clearSprings(x.size());
	
    for(unsigned int i=0;i<x.size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());	
	
	
			//initCamera();
	std::string configFile;	
	bool opt_device = false;
    bool opt_display = true;
    bool use_cuda = true;
    bool opt_click_allowed = true;
    int start_image = 0;
	
	if (!useRealData.getValue()) useGroundTruth.setValue(true);
	
	if (configFile.empty())
	configFile = vpIoTools::path("param/pizza.lua");

    /*listimg.resize(0);
    listimgseg.resize(0);
    listdepth.resize(0);*/
	
}

template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::resetSprings()
{
	
this->clearSprings(sourceVisiblePositions.getValue().size());	
for(unsigned int i=0;i<sourceVisiblePositions.getValue().size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());	
	
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::detectBorder(vector<bool> &border,const helper::vector< tri > &triangles)
{
    unsigned int nbp=border.size();
    unsigned int nbt=triangles.size();
    for(unsigned int i=0;i<nbp;i++) border[i]=false;

    if(!nbt) return;
    vector<vector< unsigned int> > ngbTriangles((int)nbp);
    for(unsigned int i=0;i<nbt;i++) for(unsigned int j=0;j<3;j++)	ngbTriangles[triangles[i][j]].push_back(i);
    for(unsigned int i=0;i<nbp;i++) if(ngbTriangles[i].size()==0) border[i]=true;
    for(unsigned int i=0;i<nbt;i++)
        for(unsigned int j=0;j<3;j++)
        {
            unsigned int id1=triangles[i][j],id2=triangles[i][(j==2)?0:j+1];
            if(!border[id1] || !border[id2]) {
                bool bd=true;
                for(unsigned int i1=0;i1<ngbTriangles[id1].size() && bd;i1++)
                    for(unsigned int i2=0;i2<ngbTriangles[id2].size() && bd;i2++)
                        if(ngbTriangles[id1][i1]!=i)
                            if(ngbTriangles[id1][i1]==ngbTriangles[id2][i2])
                                bd=false;
                if(bd) border[id1]=border[id2]=true;
            }
        }
}


template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::initSourceSurface()
{
    // build k-d tree
    const VecCoord&  p = *this->mstate->getX();
	
	ReadAccessor< Data< VecCoord > > pSurface(sourceSurfacePositions);
	ReadAccessor< Data< VecCoord > > nSurface(sourceSurfaceNormals);
	
    //const VecCoord&  pSurface = sourceSurfacePositions.getValue();
	
	//for (unsigned int i=0; i < sourceNormals.getValue().size(); i++)
	//std::cout << " source surface position " << sourceSurfaceNormals.getValue()[1][0] << std::endl;
	//std::cout << " source normals " << sourceT.getValue()[1][0] << std::endl;
	
	VecCoord nSurfaceM;
		
	sourceSurface.resize(p.size());
	Vector3 pos;
	Vector3 col;
	
	//sourceSurfaceMapping.resize(pSurface.size());
	nSurfaceM.resize(p.size());

	for (unsigned int i=0; i<p.size(); i++)
	{
			sourceSurface[i] = false;
			for (unsigned int j=0; j<pSurface.size(); j++)
			{
				Real dist=(p[i]-pSurface[j]).norm();
				if ((float)dist < 0.001)
				{
					//std::cout << " i " << i << " j " << j << "dist " << dist << std::endl; 
					sourceSurface[i] = true;
					nSurfaceM[i] = nSurface[j];
				}
			}
			//std::cout << " i " << (int)sourceSurface[i] << std::endl;
	}
			
	sourceSurfaceNormalsM.setValue(nSurfaceM);
	
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::updateSourceSurface()
{
    // build k-d tree
    const VecCoord&  p = *this->mstate->getX();
	
    const VecCoord&  pSurface = sourceSurfacePositions.getValue();
		
	sourceSurface.resize(p.size());
	Vector3 pos;
	Vector3 col;

	for (unsigned int i=0; i<p.size(); i++)
	{
		sourceSurface[i] = false;
			for (unsigned int j=0; j<pSurface.size(); j++)
			{
				Real dist=(p[i]-pSurface[j]).norm();
				if (dist < 10e-3f)
				{
					sourceSurface[i] = true;
				}
			}
	}
	
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::initSourceVisible()
{
    // build k-d tree
    //const VecCoord&  p = *this->mstate->getX();
	
	const VecCoord&  p = sourceVisiblePositions.getValue();
	
    sourceKdTree.build(p);
		
    // detect border
    if(sourceBorder.size()!=p.size()) 
	{ sourceBorder.resize(p.size()); 
	//detectBorder(sourceBorder,sourceTriangles.getValue()); 
	}
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


template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::writeImagesSynth()
{

	std::string opath3 = "out/imagesSynth30_S3_CoRot/rtt%06d.png";
 
	cv::Mat rtt,rtt1;

    //for (int frame_count = 1 ;frame_count < nimages*niterations; frame_count++)
	for (int frame_count = 1 ;frame_count < nimages-4; frame_count++)
    {	
		std::cout << " ok write rtt " << frame_count << std::endl;
		rtt = *listrtt[frame_count];
		cvtColor(rtt,rtt1 ,CV_RGB2BGR);
		char buf4[FILENAME_MAX];
        sprintf(buf4, opath3.c_str(), frame_count-1);
        std::string filename4(buf4);
		cv::imwrite(filename4,rtt1);
		//delete listrtt[frame_count];
	}
	
	
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::writeData()
{

	std::string opath = "in/images4/img1%06d.png";
	std::string opath1 = "in/images4/depthfile%06d.txt";
	std::string extensionfile = "in/images4/extension.txt";
 
    std::vector<Vec3d> pcd1;
	std::vector<bool> visible1;
	
	ofstream extfile(extensionfile.c_str(), ios::out | ios::binary );
	double extension;
	double pos = 0;
	double pos0;


	cv::Mat rtt,rtt1,rtt2;	
    for (int frame_count = 0 ;frame_count < listpcd.size(); frame_count++)
    {
		
		std::cout << " ok write " << frame_count << std::endl;
    	pcd1 = *listpcd[frame_count];
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

		writePCDToFile(filename1,pcd1,visible1);
		
		/*delete listimg[frame_count];
		delete listimgseg[frame_count];
		delete listdepth[frame_count];*/
    }
	extfile.close();


    //for (int frame_count = 1 ;frame_count < nimages*niterations; frame_count++)
	for (int frame_count = 1 ;frame_count < listrtt.size(); frame_count++)
    {	
		std::cout << " ok write rtt " << frame_count << std::endl;
		rtt = *listrtt[frame_count];
		cvtColor(rtt,rtt1 ,CV_RGB2BGR);
		cv::flip(rtt1,rtt2,0);
		char buf4[FILENAME_MAX];
        sprintf(buf4, opath.c_str(), frame_count-1);
        std::string filename4(buf4);
		cv::imwrite(filename4,rtt2);
		//delete listrtt[frame_count];
	}
	
	
}


void renderToTexture(cv::Mat &_rtt)
{
	_rtt.create(480,640, CV_8UC3);
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    //img.init(viewport[2], viewport[3], 1, 1, io::Image::UNORM8, io::Image::RGB);
    glReadBuffer(GL_FRONT);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3], GL_RGB, GL_UNSIGNED_BYTE, _rtt.data);
    glReadBuffer(GL_BACK);
	//cv::imwrite("rtt.png",_rtt);
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::renderToTextureD(cv::Mat &_rtt)
{
	_rtt.create(480,640, CV_8UC3);
	cv::Mat _rttd,_dd,_rttd_;
	_rttd.create(240,320, CV_8UC3);
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    //img.init(viewport[2], viewport[3], 1, 1, io::Image::UNORM8, io::Image::RGB);
    glReadBuffer(GL_FRONT);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3], GL_RGB, GL_UNSIGNED_BYTE, _rtt.data);
    glReadBuffer(GL_BACK);
	//cv::imwrite("rtt.png",_rtt);
	
	int hght = 480;
	int wdth = 640;
	//_rttd = _rtt;
		
	cv::resize(_rtt, _rttd_,_rttd.size());
	cv::flip( _rttd_,_rttd,0);
		
	//cv::namedWindow("rttd");
    //cv::imshow("rttd",_rtt);
	
	cv::Mat _rtd,_rtd0;
	_rtd.create(hght, wdth, CV_32F);
	_rtd0.create(hght,wdth, CV_8UC1);
	_dd.create(240,320, CV_8UC1);
	
    GLfloat depths[hght * wdth ];
	//depths = new GLfloat[240 * 320 ];
	
    glReadBuffer(GL_FRONT);
	glEnable(GL_DEPTH_TEST);

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
	std::cout << " no viewer 1" << std::endl; 
    glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3], GL_DEPTH_COMPONENT, GL_FLOAT, &depths);
    glReadBuffer(GL_BACK);
	
	int offset = 5;
	 
	 std::cout << " ok read " << std::endl;
	
	for (int j = 0; j < wdth; j++)
		for (int i = 0; i< hght; i++)
		{
			if ((double)(float)depths[j+i*wdth]	< 1)
			{
				_rtd0.at<uchar>(hght-i-1,j) = 255;
			}
			else _rtd0.at<uchar>(hght-i-1,j) = 0;
		}
		cv::resize(_rtd0, _dd, Size(320, 240));

	 std::cout << " ok read 1 " << std::endl;		
		
	for (int j = 0; j < 320; j++)
		for (int i = 0; i< 240; i++)
		{
			if (_dd.at<uchar>(i,j) == 0)
			{
				_rttd.at<Vec3b>(i,j)[0] = color_1.at<Vec3b>(i,j)[2];
				_rttd.at<Vec3b>(i,j)[2] = color_1.at<Vec3b>(i,j)[0];
				_rttd.at<Vec3b>(i,j)[1] = color_1.at<Vec3b>(i,j)[1];
			}
		}
		
		cv::namedWindow("dd");
		cv::imshow("dd",_rttd);
		_rtt = _rttd;
		
}


template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::initSegmentation()
{	
	cv::Mat mask,maskimg,mask0,roimask,mask1; // segmentation result (4 possible values)
    cv::Mat bgModel,fgModel; // the models (internally used)

	//cv::imwrite("color.png",color);
	
	cv::Mat downsampledbox,downsampled;
	
	//cv::pyrDown(color, downsampledbox, cv::Size(color.cols/2, color.rows/2));
	downsampledbox = color;
	
	IplImage temp;
    IplImage tempgs;
    cv::Mat imgs,imgklt;
    temp = downsampledbox;
    IplImage* temp1 = cvCloneImage(&temp);
	
	cout << "init seg" << endl;	
	
    const char* name = "image";

    cv::namedWindow(name);
    box = cvRect(0,0,1,1);

	//cv::imshow("image",	color);
	//tempm.resize(image.step1());
	  
    // Set up the callback
    cv::setMouseCallback(name, my_mouse_callback, (void*) &temp);
	
		/*for (int i = 0; i<color.rows; i++)
		  for (int j = 0; j<color.cols; j++)
	  std::cout << (int)color.at<Vec3b>(i,j)[0] << std::endl;*/
	  
	  std::cout << " Time " << (int)this->getContext()->getTime() << std::endl;

    // Main loop
    while(1)
    {
      if (destroy)
      {
        cvDestroyWindow(name); break;
      }
      cvCopyImage(&temp, temp1);

      if (drawing_box)
          draw_box(temp1, box);
		  
      cv::moveWindow(name, 200, 100);
      cvShowImage(name, temp1);
      //tempm.resize(image.step1());
      int key=cvWaitKey(10);
      if ((char)key == 27) break;
	
    }

    cvReleaseImage(&temp1);
    cv::setMouseCallback(name, NULL, NULL);
	
	    //cvDestroyWindow(name);
    cv::Rect rectangle(69,47,198,171);
    cv::Rect recttemp = rectangle;

    rectangle = box;
    seg.setRectangle(rectangle);

    int width = downsampled.cols;
    int height = downsampled.rows;

    foreground = cv::Mat(downsampled.size(),CV_8UC3,cv::Scalar(255,255,255));

    // GrabCut segmentation
	
	 if (useSensor.getValue())
	 getImages();
	 else readImages();
	
	//cv::pyrDown(color, downsampled, cv::Size(color.cols/2, color.rows/2));
	downsampled = color;

    seg.segmentationFromRect(downsampled,foreground);

    // draw rectangle on original image
    //cv::rectangle(image, rectangle, cv::Scalar(255,255,255),1);
    
	//cv::namedWindow("Image");
    //cv::imshow("Image",downsampled);

    // display result
    cv::namedWindow("Segmented Image");
    cv::imshow("Segmented Image",foreground);

    if (waitKey(20) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
   {
        cout << "esc key is pressed by user" << endl;
   }
   
   	/*imglsg = new cv::Mat;
    *imglsg = foreground;
	
    listimgseg.push_back(imglsg);
    delete imglsg;*/

}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::segment()
{

	cv::Mat downsampled;
	//cv::pyrDown(color, downsampled, cv::Size(color.cols/2, color.rows/2));
	downsampled = color;
	
	//cv::imwrite("downsampled.png",downsampled);
	
	int width = downsampled.cols;
    int height = downsampled.rows;
	//cv::imwrite("foreground0.png",foreground);

		
	seg.updateMask(foreground);
	
	seg.updateSegmentation(downsampled,foreground);
	//seg.updateSegmentationCrop(downsampled,foreground);
	//cv::imwrite("foreground1.png",foreground);

	
	/*cv::namedWindow("Image");
    cv::imshow("Image",downsampled);*/

    // display result
    cv::namedWindow("Segmented Image");
    cv::imshow("Segmented Image",foreground);
	
	//cv::imwrite("color.png",color);


	
	imglsg = new cv::Mat;
    *imglsg = foreground;
	
    listimgseg.push_back(imglsg);
    //delete imglsg;
	
	imgl = new cv::Mat;
    //*imgl = color;
	*imgl = downsampled;		
	depthl = new cv::Mat;
    *depthl = depth;

	
    listimg.push_back(imgl);
    listdepth.push_back(depthl);
	
	//delete depthl;
    //delete imgl;

}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::segmentSynth()
{

	cv::Mat downsampled;
	//cv::pyrDown(color, downsampled, cv::Size(color.cols/2, color.rows/2));
	downsampled = color;
	
	foreground = cv::Mat(downsampled.size(),CV_8UC4,cv::Scalar(255,255,255,0));
	
	cv::imwrite("downsampled0.png",downsampled);
	
	int width = downsampled.cols;
    int height = downsampled.rows;
	//cv::imwrite("foreground0.png",foreground);
	
	std::cout << " ok seg synth " << std::endl; 

	for (int j = 0; j < width; j++)
	  for (int i = 0; i < height; i++)
		{
				foreground.at<cv::Vec4b>(i,j)[0] = color.at<cv::Vec3b>(i,j)[0];
				foreground.at<cv::Vec4b>(i,j)[1] = color.at<cv::Vec3b>(i,j)[1];
				foreground.at<cv::Vec4b>(i,j)[2] = color.at<cv::Vec3b>(i,j)[2];
				
			if (color.at<cv::Vec3b>(i,j)[0] == 0 && color.at<cv::Vec3b>(i,j)[1] == 1 && color.at<cv::Vec3b>(i,j)[2] == 2)
			{
				foreground.at<cv::Vec4b>(i,j)[3] = 0;
			}
			else foreground.at<cv::Vec4b>(i,j)[3] = 255;
		}

	cv::Mat distanceMap(foreground.size(),CV_8U,cv::Scalar(255));
	cv::Mat dot(foreground.size(),CV_8U,cv::Scalar(0));
	
		//std::cout << " ok seg synth 1 " << std::endl; 

	seg.filter(foreground,distanceMap,dot);	
	
	//cv::imwrite("downsampled1.png",distanceMap);
	
	//distImage = distanceMap;
	//dotImage = dot;
	
	/*cv::namedWindow("Image");
    cv::imshow("Image",downsampled);*/

    // display result
    cv::namedWindow("Segmented Image");
    cv::imshow("Segmented Image",foreground);
	
	cv::imwrite("color.png",color);

	imglsg = new cv::Mat;
    *imglsg = foreground;
	
    listimgseg.push_back(imglsg);
    //delete imglsg;
	
	imgl = new cv::Mat;
    //*imgl = color;
	*imgl = downsampled;		

	
    listimg.push_back(imgl);	
	//delete depthl;
    //delete imgl;

}


template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::extractTargetContour()
{
	cv::Mat distImage, contourImage;
	targetContour.reset(new pcl::PointCloud<pcl::PointXYZRGB>); 
	//pcl::PointCloud<pcl::PointXYZRGB> target;
	
	//cv::imwrite("dist_image.png",seg.distImage);
	//cv::imwrite("dot_image.png",seg.dotImage);

	targetContour = PCDContour(depth,seg.distImage,seg.dotImage);
	
	/*cv::namedWindow("dist_image");
	cv::imshow("dist_image",seg.distImage);
	
	cv::namedWindow("dot_image");
	cv::imshow("dot_image",seg.dotImage);*/
		
	VecCoord targetpos;
	targetpos.resize(targetContour->size());

std::cout << " target contour size " << targetContour->size() << std::endl;
            Vector3 pos;
            Vector3 col;

	for (unsigned int i=0; i<targetContour->size(); i++)
	{
            pos[0] = (double)targetContour->points[i].x;
            pos[1] = (double)targetContour->points[i].y;
            pos[2] = (double)targetContour->points[i].z;
            targetpos[i]=pos;
			//std::cout << " x " << pos[0] << " y " << pos[1] << " z " << pos[2] << std::endl; 
	} 

    const VecCoord&  p = targetpos;
	targetContourPositions.setValue(p);	
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::setViewPoint()
{
	Eigen::Affine3f scene_sensor_pose (Eigen::Affine3f::Identity ());
	pcl::PointCloud<pcl::PointXYZRGB>& point_cloud = *target;
	
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

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::setViewPointData()
{
	Eigen::Affine3f scene_sensor_pose (Eigen::Affine3f::Identity ());
	
	source.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
	const VecCoord& x = *this->mstate->getX();
	
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

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::getSourceVisible()
{

	int hght = 480;
	int wdth = 640;
	
	std::cout << " no viewer 1" << std::endl; 
		
	sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
	sofa::component::visualmodel::BaseCamera::SPtr currentCamera;// = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();
	root->get(currentCamera);
	
	cv::Mat _rtd,_rtd0;
	_rtd.create(hght, wdth, CV_32F);
	_rtd0.create(hght,wdth, CV_8UC1);
	
    GLfloat depths[hght * wdth ];
	GLfloat depthsN[hght * wdth ];
	//depths = new GLfloat[240 * 320 ];
	
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    //img.init(viewport[2], viewport[3], 1, 1, io::Image::UNORM8, io::Image::RGB);
    //glReadBuffer(GL_FRONT);
	glEnable(GL_DEPTH_TEST);

    //glPixelStorei(GL_PACK_ALIGNMENT, 1);
	std::cout << " no viewer 1" << std::endl; 
			
	//glDepthMask(GL_TRUE);
	//glDepthFunc(GL_ALWAYS); // Change this to whatever kind of depth testing you want
    //glDepthRange(0.0f, 1.0f);

    glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3], GL_DEPTH_COMPONENT, GL_FLOAT, &depths);
    //glReadBuffer(GL_BACK);
	
	double znear = currentCamera->getZNear();
	double zfar = currentCamera->getZFar();
	
	//double zNear = 1e10;
	//double zFar = -1e10;
	
	std::cout << " znear " << znear << std::endl;	
	std::cout << " zfar " << zfar << std::endl;
	
	for (int j = 0; j < wdth; j++)
		for (int i = 0; i< hght; i++)
		{
			if ((double)(float)depths[j+i*wdth]	< 1)
			{
			_rtd0.at<uchar>(hght-i-1,j) = 255;//(int)100000*(1-depths[j+i*wdth]);
			
			double clip_z = (depths[j+i*wdth] - 0.5) * 2.0;
            depthsN[j+i*wdth] = 2*znear*zfar/(clip_z*(zfar-znear)-(zfar+znear));
			//std::cout << " depth " << world_z << std::endl;

			}
			else 
				{
				_rtd0.at<uchar>(hght-i-1,j) = 0;
				depthsN[j+i*wdth] = 0;
				}
				//std::cout << " depth " << (int)_rtd0.at<uchar>(i,j) << std::endl;
		//if (depths[j+i*319]	> 0)
			//_rtd0.at<uchar>(j,i) = 255;
		}
		
		depthMap = _rtd0;
		
		//cv::imwrite("depth.png", _rtd0);
	
	std::cout << " no viewer 2" << std::endl;
	
const VecCoord& x = *this->mstate->getX();
 	
Eigen::Matrix3f rgbIntrinsicMatrix;
rgbIntrinsicMatrix(0,0) = 2*300.34;
rgbIntrinsicMatrix(1,1) = 2*300.34;
rgbIntrinsicMatrix(0,2) = 2*160;
rgbIntrinsicMatrix(1,2) = 2*120;

sourceVisible.resize(x.size());
int nvisible = 0;


for (int k = 0; k < x.size(); k++)
{
	int x_u = (int)(x[k][0]*rgbIntrinsicMatrix(0,0)/x[k][2] + rgbIntrinsicMatrix(0,2));
    int x_v = (int)(x[k][1]*rgbIntrinsicMatrix(1,1)/x[k][2] + rgbIntrinsicMatrix(1,2));
	
	//std::cout << " depths " << x_u << " " << x_v << std::endl;
	//std::cout << " depths " << (float)depthsN[x_u+(hght-x_v-1)*wdth] << " " <<(float)x[k][2] << std::endl;

	
if ((float)abs(depthsN[x_u+(hght-x_v-1)*wdth]+(float)x[k][2]) < 1*10e-3f || (float)depthsN[x_u+(hght-x_v-1)*wdth] == 0)
	{
	sourceVisible[k] = true;
	nvisible ++;
	}
	else {
	sourceVisible[k] = false;	
	}
	
}

std::cout << " npoints " << x.size() << " " << nvisible << std::endl;

VecCoord sourceVis;
sourceVis.resize(nvisible);
indicesVisible.resize(nvisible);

std::cout << " nvisible " << nvisible << " xsize " << sourceVisible.size() <<  std::endl;

Vector3 pos;
int k = 0;
			
for (unsigned int i=0; i< x.size(); i++)
	{
		if (sourceVisible[i])
		{
            pos = x[i];
            sourceVis[k]=pos;
			indicesVisible[k] = i; 
			k++;			
		}
	}
sourceVisiblePositions.setValue(sourceVis);

//resetSprings();

}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::extractSourceContour()
{
	
	getSourceVisible();
	
	
	int hght = 480;
	int wdth = 640;
		
	double cannyTh1 = 350;
	double cannyTh2 = 10;
	cv::Mat contour,dist,dist0;
	
	//cv::imwrite("depthmap.png", depthMap);
	
	cv::Canny( depthMap, contour, cannyTh1, cannyTh2, 3);
    contour = cv::Scalar::all(255) - contour;

	cv::distanceTransform(contour, dist, CV_DIST_L2, 3);
//dt0 *= 5000;
//pow(dt0, 0.5, dt0);
    dist.convertTo(dist0, CV_8U, 1, 0);
	
	int ncontour = 0;
	pcl::PointCloud<pcl::PointXYZRGB> sourceContour;
	
	Eigen::Matrix3f rgbIntrinsicMatrix;

	rgbIntrinsicMatrix(0,0) = 2*300.34;
	rgbIntrinsicMatrix(1,1) = 2*300.34;
	rgbIntrinsicMatrix(0,2) = 2*160;
	rgbIntrinsicMatrix(1,2) = 2*120;
	
	const VecCoord& x = *this->mstate->getX();
	
	unsigned int nbs=x.size();
	
	cv::Mat contourpoints = cv::Mat::zeros(480,640, CV_8UC3);
	contourpoints= cv::Mat(480,640,CV_8UC3,cv::Scalar(255,255,255));
		for (int j = 0; j < wdth; j++)
		for (int i = 0; i< hght; i++)
		{
		contourpoints.at<Vec3b>(i,j)[0]= contour.at<uchar>(i,j);
		}
	
	pcl::PointXYZRGB newPoint;
	
	sourceBorder.resize(nbs);
	
	//cv::imwrite("dist0.png", dist0);
	
	int nsourcecontour = 0;
	
	sourceWeights.resize(0);
	
	double totalweights = 0;
	
	for (unsigned int i=0; i<nbs; i++)
	{
		int x_u = (int)(x[i][0]*rgbIntrinsicMatrix(0,0)/x[i][2] + rgbIntrinsicMatrix(0,2));
		int x_v = (int)(x[i][1]*rgbIntrinsicMatrix(1,1)/x[i][2] + rgbIntrinsicMatrix(1,2));
		int thickness = 1;
        int lineType = 2;

        circle( contourpoints,
         cv::Point(x_u,x_v),
         wdth/128.0,
         Scalar( 0, 0, 255 ),
         thickness,
         lineType );
				if (dist0.at<uchar>(x_v,x_u) < 5)
				{
				newPoint.z = x[i][2];
				newPoint.x = x[i][0];
				newPoint.y = x[i][1];
				newPoint.r = 0;
				newPoint.g = 0;
				newPoint.b = 0;
				sourceContour.push_back(newPoint);
				sourceBorder[i] = true;
				nsourcecontour++;
				}
				else sourceBorder[i] = false;
				
				//sourceWeights.push_back((double)1./(0.12*(1.0+sqrt(dist0.at<uchar>(x_v,x_u)))));
				
				sourceWeights.push_back((double)exp(-dist0.at<uchar>(x_v,x_u)/10.0));

				//std::cout << " weight " << (double)1./0.12*(1.0+sqrt(bvalue)) << std::endl;
				
				//targetWeights.push_back((double)0.4/(1.0+bvalue*bvalue));
				
				/*if (avalue > 0 && bvalue < 6)
				{
				targetWeights.push_back((double)3);
				}
				else targetWeights.push_back((double)3);*/
				
				//totalweights += (double)1./(0.12*(1.0+sqrt(bvalue)));
				totalweights += sourceWeights[i];
	
	}
	
	for (int i=0; i < sourceWeights.size();i++)
	{
		sourceWeights[i]*=((double)sourceWeights.size()/totalweights);
		//std::cout << " weights " << (double)targetWeights[i] << std::endl;
	}
	
	//cv::imwrite("contourpoints.png", contourpoints);
	//cv::imwrite("dist.png", dist);
	std::cout << " n source contour " << nsourcecontour << std::endl;

	
	//getchar();


		/*for (int j = 0; j < 320; j++)
		for (int i = 0;i<240; i++)
		{
			if (dist.at<uchar> > 3 && _rtd > 0)
				
				
				ncontour++;
		}*/
	
	VecCoord sourcecontourpos;
    sourcecontourpos.resize(sourceContour.size());
            Vector3 pos;

	for (unsigned int i=0; i<sourcecontourpos.size(); i++)
	{
            pos[0] = sourceContour[i].x;
            pos[1] = sourceContour[i].y;
            pos[2] = sourceContour[i].z;
            sourcecontourpos[i]=pos;
	}
    const VecCoord&  p = sourcecontourpos;
	sourceContourPositions.setValue(p);
	
}	
	
/*template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::extractSourceContour()
{
	Eigen::Affine3f scene_sensor_pose (Eigen::Affine3f::Identity ());
	pcl::PointCloud<pcl::PointXYZRGB>& point_cloud = *target;

	scene_sensor_pose = Eigen::Affine3f (Eigen::Translation3f (point_cloud.sensor_origin_[0],
                                                             point_cloud.sensor_origin_[1],
                                                             point_cloud.sensor_origin_[2])) * 
															 Eigen::Affine3f (point_cloud.sensor_orientation_);
															 
  Eigen::Affine3f viewer_pose = scene_sensor_pose;
  Eigen::Vector3f pos_vector = viewer_pose * Eigen::Vector3f(0, 0, 0);  
  Eigen::Vector3f look_at_vector = viewer_pose.rotation () * Eigen::Vector3f(0, 0, 1) + pos_vector;
  Eigen::Vector3f up_vector = viewer_pose.rotation () * Eigen::Vector3f(0, -1, 0);
							
	int hght = 480;
	int wdth = 640;
							
    //sofa::simulation::Node* root = dynamic_cast<simulation::Node*>(this->getContext());
	//if(root)
		{
		//sofa::component::visualmodel::InteractiveCamera* currentCamera = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();
			//if(currentCamera)
			
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
	cv::Mat _rtt;		
	//renderToTexture(_rtt);
	
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    //img.init(viewport[2], viewport[3], 1, 1, io::Image::UNORM8, io::Imag
		}
	cv::Mat _rtd,_rtd0;
	_rtd.create(hght, wdth, CV_32F);
	_rtd0.create(hght,wdth, CV_8UC1);
	
    GLfloat depths[hght * wdth ];
	//depths = new GLfloat[240 * 320 ];
	
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    //img.init(viewport[2], viewport[3], 1, 1, io::Image::UNORM8, io::Image::RGB);
    glReadBuffer(GL_FRONT);
	glEnable(GL_DEPTH_TEST);

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
	        std::cout << " no viewer 1" << std::endl; 

    glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3], GL_DEPTH_COMPONENT, GL_FLOAT, &depths);
    glReadBuffer(GL_BACK);
	
	std::cout << " no viewer 2" << std::endl;
 
	
	for (int j = 0; j < wdth; j++)
		for (int i = 0; i< hght; i++)
		{
			if ((double)(float)depths[j+i*wdth]	< 1)
			{
				//std::cout << " depth " << (double)(float)depths[i+j*hght] << std::endl;
				_rtd0.at<uchar>(hght-i-1,j) = 255;
			}
			else _rtd0.at<uchar>(hght-i-1,j) = 0;
				//std::cout << " depth " << (int)_rtd0.at<uchar>(i,j) << std::endl;
		//if (depths[j+i*319]	> 0)
			//_rtd0.at<uchar>(j,i) = 255;
		}
		
	double cannyTh1 = 350;
	double cannyTh2 = 10;
	cv::Mat contour,dist,dist0;
	
	cv::Canny( _rtd0, contour, cannyTh1, cannyTh2, 3);
    contour = cv::Scalar::all(255) - contour;

	cv::distanceTransform(contour, dist, CV_DIST_L2, 3);
//dt0 *= 5000;
//pow(dt0, 0.5, dt0);
    dist.convertTo(dist0, CV_8U, 1, 0);
	
	int ncontour = 0;
	pcl::PointCloud<pcl::PointXYZRGB> sourceContour;
	
	Eigen::Matrix3f rgbIntrinsicMatrix;

	rgbIntrinsicMatrix(0,0) = 2*300.34;
	rgbIntrinsicMatrix(1,1) = 2*300.34;
	rgbIntrinsicMatrix(0,2) = 2*160;
	rgbIntrinsicMatrix(1,2) = 2*120;
	
	const VecCoord& x = *this->mstate->getX();
	
	unsigned int nbs=x.size();
	
	cv::Mat contourpoints = cv::Mat::zeros(480,640, CV_8UC3);
	contourpoints= cv::Mat(480,640,CV_8UC3,cv::Scalar(255,255,255));
		for (int j = 0; j < wdth; j++)
		for (int i = 0; i< hght; i++)
		{
	contourpoints.at<Vec3b>(i,j)[0]= contour.at<uchar>(i,j);
		}
	
	pcl::PointXYZRGB newPoint;
	
	sourceBorder.resize(nbs);
	
	for (unsigned int i=0; i<nbs; i++)
	{
		int x_u = (int)(x[i][0]*rgbIntrinsicMatrix(0,0)/x[i][2] + rgbIntrinsicMatrix(0,2));
		int x_v = (int)(x[i][1]*rgbIntrinsicMatrix(1,1)/x[i][2] + rgbIntrinsicMatrix(1,2));
		int thickness = 1;
        int lineType = 2;

        circle( contourpoints,
         cv::Point(x_u,x_v),
         wdth/128.0,
         Scalar( 0, 0, 255 ),
         thickness,
         lineType );
				if (dist0.at<uchar>(x_v,x_u) < 6)
				{
				newPoint.z = x[i][2];
				newPoint.x = x[i][0];
				newPoint.y = x[i][1];
				newPoint.r = 0;
				newPoint.g = 0;
				newPoint.b = 0;
				sourceContour.push_back(newPoint);
				sourceBorder[i] = true;
				}
				else sourceBorder[i] = false;
	}
	
	//cv::imwrite("contourpoints.png", contourpoints);
	//cv::imwrite("dist.png", dist);

	
	//getchar();	
	VecCoord sourcecontourpos;
    sourcecontourpos.resize(sourceContour.size());
            Vector3 pos;

	for (unsigned int i=0; i<sourcecontourpos.size(); i++)
	{
            pos[0] = sourceContour[i].x;
            pos[1] = sourceContour[i].y;
            pos[2] = sourceContour[i].z;
            sourcecontourpos[i]=pos;
	}
    const VecCoord&  p = sourcecontourpos;
	sourceContourPositions.setValue(p);
		
}*/	

template<class DataTypes, class DepthTypes>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr GenerateSyntheticData<DataTypes, DepthTypes>::PCDFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	//pcl::PointCloud<pcl::PointXYZRGB> pointcloud;
	
	cv::Mat frgd;
	
	//cv::pyrDown(rgbImage, frgd, cv::Size(foreground.cols/2, foreground.rows/2));
	frgd = rgbImage;	

	//typename pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud;
	Eigen::Matrix3f rgbIntrinsicMatrix;
	/*rgbIntrinsicMatrix(0,0) = 120.53;
	rgbIntrinsicMatrix(1,1) = 120.146;
	rgbIntrinsicMatrix(0,2) = 80;
	rgbIntrinsicMatrix(1,2) = 60;*/
	
	//rgbIntrinsicMatrix(0,0) = 570.34;
	//rgbIntrinsicMatrix(1,1) = 570.34;
	rgbIntrinsicMatrix(0,0) = 300.34;
	rgbIntrinsicMatrix(1,1) = 300.34;
	rgbIntrinsicMatrix(0,2) = 157.25;
	rgbIntrinsicMatrix(1,2) = 117.75;
		
    //cv::imwrite("depth.png",depthImage);
    //cv::imwrite("RGB.png",frgd);

	float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
	float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy
	pcl::PointXYZRGB newPoint;
	int sample = 4;
	int offsetx = 3;
	int offsety = 0;
	for (int i=0;i<(int)(depthImage.rows-offsety)/sample;i++)
	{
		for (int j=0;j<(int)(depthImage.cols-offsetx)/sample;j++)
		{
			float depthValue = (float)depthImage.at<float>(sample*(i+offsety),sample*(j+offsetx));
			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)frgd.at<Vec4b>(sample*i,sample*j)[3];
			if (avalue > 0 && depthValue>0)                // if depthValue is not NaN
			{
				// Find 3D position respect to rgb frame:
				newPoint.z = depthValue;
				newPoint.x = (sample*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
				newPoint.y = (sample*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
				newPoint.r = frgd.at<cv::Vec4b>(sample*i,sample*j)[2];
				newPoint.g = frgd.at<cv::Vec4b>(sample*i,sample*j)[1];
				newPoint.b = frgd.at<cv::Vec4b>(sample*i,sample*j)[0];
				outputPointcloud->points.push_back(newPoint);
			}
			/*else
			{
				newPoint.z = std::numeric_limits<float>::quiet_NaN();
				newPoint.x = std::numeric_limits<float>::quiet_NaN();
				newPoint.y = std::numeric_limits<float>::quiet_NaN();
				newPoint.r = std::numeric_limits<unsigned char>::quiet_NaN();
				newPoint.g = std::numeric_limits<unsigned char>::quiet_NaN();
				newPoint.b = std::numeric_limits<unsigned char>::quiet_NaN();
				//outputPointcloud.push_back(newPoint);
			}*/
		}
	}
	/*const std::string file = "test_pcdf%03d.pcd";
    char buf[FILENAME_MAX];
    sprintf(buf, file.c_str(), frame);
    std::string filename(buf);*/
	//pcl::io::savePCDFileASCII (filename, outputPointcloud);
	//std::cout << "Saved " << outputPointcloud.points.size () << " data points to test_pcd.pcd." << std::endl;
		
	return outputPointcloud;
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::PCDFromRGBD2(cv::Mat& depthImage, cv::Mat& rgbImage, pcl::PointCloud<pcl::PointXYZRGB>& outputPointcloud)
{
	//pcl::PointCloud<pcl::PointXYZRGB> pointcloud;
	
	cv::Mat frgd;
	cv::pyrDown(rgbImage, frgd, cv::Size(foreground.cols/2, foreground.rows/2));

	//typename pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud;
	Eigen::Matrix3f rgbIntrinsicMatrix;
	rgbIntrinsicMatrix(0,0) = 166.53;
	rgbIntrinsicMatrix(1,1) = 166.146;
	rgbIntrinsicMatrix(0,2) = 80;
	rgbIntrinsicMatrix(1,2) = 60;
    //cv::imwrite("depth.png",depthImage);
    //cv::imwrite("RGB.png",frgd);


        float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
	float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy
	pcl::PointXYZRGB newPoint;
	int sample = 4;
	for (int i=0;i<(int)depthImage.rows/sample;i++)
	{
		for (int j=0;j<(int)depthImage.cols/sample;j++)
		{
			float depthValue = depthImage.at<float>(sample*i,sample*j);
			int avalue = (int)frgd.at<Vec4b>(sample*i,sample*j)[3];
			
			//std::cout << " x " << avalue << std::endl;  
			if (avalue > 0 && depthValue>0)                // if depthValue is not NaN
			{
				// Find 3D position respect to rgb frame:
				newPoint.z = depthValue;
				newPoint.x = (sample*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
				newPoint.y = (sample*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
				newPoint.r = frgd.at<cv::Vec4b>(sample*i,sample*j)[2];
				newPoint.g = frgd.at<cv::Vec4b>(sample*i,sample*j)[1];
				newPoint.b = frgd.at<cv::Vec4b>(sample*i,sample*j)[0];
				outputPointcloud.points.push_back(newPoint);

			}
			/*else
			{
				newPoint.z = std::numeric_limits<float>::quiet_NaN();
				newPoint.x = std::numeric_limits<float>::quiet_NaN();
				newPoint.y = std::numeric_limits<float>::quiet_NaN();
				newPoint.r = std::numeric_limits<unsigned char>::quiet_NaN();
				newPoint.g = std::numeric_limits<unsigned char>::quiet_NaN();
				newPoint.b = std::numeric_limits<unsigned char>::quiet_NaN();
				//outputPointcloud.push_back(newPoint);
			}*/
		}
	}
	/*const std::string file = "test_pcdf%03d.pcd";
    char buf[FILENAME_MAX];
    sprintf(buf, file.c_str(), frame);
    std::string filename(buf);*/
	//pcl::io::savePCDFileASCII (filename, outputPointcloud);
	//std::cout << "Saved " << outputPointcloud.points.size () << " data points to test_pcd.pcd." << std::endl;
}

template<class DataTypes, class DepthTypes>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr GenerateSyntheticData<DataTypes, DepthTypes>::PCDFromRGBD3(cv::Mat& depthImage, cv::Mat& rgbImage)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	//pcl::PointCloud<pcl::PointXYZRGB> pointcloud;
	
	cv::Mat frgd;
	
	//cv::pyrDown(rgbImage, frgd, cv::Size(foreground.cols/2, foreground.rows/2));
	frgd = rgbImage;	

	//typename pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud;
	Eigen::Matrix3f rgbIntrinsicMatrix;
	/*rgbIntrinsicMatrix(0,0) = 120.53;
	rgbIntrinsicMatrix(1,1) = 120.146;
	rgbIntrinsicMatrix(0,2) = 80;
	rgbIntrinsicMatrix(1,2) = 60;*/
	
	//rgbIntrinsicMatrix(0,0) = 570.34;
	//rgbIntrinsicMatrix(1,1) = 570.34;
	rgbIntrinsicMatrix(0,0) = 300.34;
	rgbIntrinsicMatrix(1,1) = 300.34;
	rgbIntrinsicMatrix(0,2) = 157.25;
	rgbIntrinsicMatrix(1,2) = 117.75;
		
    //cv::imwrite("depth.png",depthImage);
    //cv::imwrite("RGB.png",frgd);

	float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
	float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy
	pcl::PointXYZRGB newPoint;
	int sample = 4;
	int offsetx = 3;
	for (int i=0;i<(int)depthImage.rows/sample;i++)
	{
		for (int j=0;j<(int)(depthImage.cols-offsetx)/sample;j++)
		{
			float depthValue = depthImage.at<float>(sample*(i),sample*(j+offsetx));
			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)frgd.at<Vec4b>(sample*i,sample*j)[3];
			
			//std::cout << " x " << depthValue << std::endl;  
			if (depthValue>0)  // if depthValue is not NaN
			 if (avalue > 0) 
			  {
				// Find 3D position respect to rgb frame:
				newPoint.z = depthValue;
				newPoint.x = (sample*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
				newPoint.y = (sample*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
				newPoint.r = 1;
				newPoint.g = 0;
				newPoint.b = 0;
				outputPointcloud->points.push_back(newPoint);
			  }
			  /*else
			  {
								// Find 3D position respect to rgb frame:
				newPoint.z = depthValue;
				newPoint.x = (sample*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
				newPoint.y = (sample*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
				newPoint.r = 0;
				newPoint.g = 0;
				newPoint.b = 0;
				outputPointcloud->points.push_back(newPoint);  
			  }*/
			/*else
			{
				newPoint.z = std::numeric_limits<float>::quiet_NaN();
				newPoint.x = std::numeric_limits<float>::quiet_NaN();
				newPoint.y = std::numeric_limits<float>::quiet_NaN();
				newPoint.r = std::numeric_limits<unsigned char>::quiet_NaN();
				newPoint.g = std::numeric_limits<unsigned char>::quiet_NaN();
				newPoint.b = std::numeric_limits<unsigned char>::quiet_NaN();
				//outputPointcloud.push_back(newPoint);
			}*/
		}
	}
	/*const std::string file = "test_pcdf%03d.pcd";
    char buf[FILENAME_MAX];
    sprintf(buf, file.c_str(), frame);
    std::string filename(buf);*/
	//pcl::io::savePCDFileASCII (filename, outputPointcloud);
	//std::cout << "Saved " << outputPointcloud.points.size () << " data points to test_pcd.pcd." << std::endl;
		
	return outputPointcloud;
}

template<class DataTypes, class DepthTypes>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr GenerateSyntheticData<DataTypes, DepthTypes>::PCDContour(cv::Mat& depthImage, cv::Mat& distImage, cv::Mat& dotImage)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	//pcl::PointCloud<pcl::PointXYZRGB> pointcloud;
	
	cv::Mat distimg,dotimg;
	
	//cv::pyrDown(rgbImage, frgd, cv::Size(foreground.cols/2, foreground.rows/2));
	distimg = distImage;
	dotimg = dotImage;

	//typename pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud;
	Eigen::Matrix3f rgbIntrinsicMatrix;
	
	//rgbIntrinsicMatrix(0,0) = 570.34;
	//rgbIntrinsicMatrix(1,1) = 570.34;
	rgbIntrinsicMatrix(0,0) = 300.34;
	rgbIntrinsicMatrix(1,1) = 300.34;
	rgbIntrinsicMatrix(0,2) = 157.25;
	rgbIntrinsicMatrix(1,2) = 117.75;
			
    //cv::imwrite("depth.png",depthImage);
    //cv::imwrite("RGB.png",frgd);

	float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
	float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy
	pcl::PointXYZRGB newPoint;
	int sample = 1;
	int offsetx = 3;
	for (int i=0;i<(int)depthImage.rows/sample;i++)
	{
		for (int j=0;j<(int)(depthImage.cols-offsetx)/sample;j++)
		{
			float depthValue = depthImage.at<float>(sample*(i),sample*(j+offsetx));
			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)distimg.at<uchar>(sample*i,sample*(j));
			int dvalue = (int)dotimg.at<uchar>(sample*i,sample*(j));

			
			//std::cout << " x " << depthValue << std::endl;  
			if (dvalue == 0 && avalue == 4 && depthValue>0)                // if depthValue is not NaN
			{
				// Find 3D position respect to rgb frame:
				newPoint.z = depthValue;
				newPoint.x = (sample*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
				newPoint.y = (sample*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
				newPoint.r = 0;//frgd.at<cv::Vec4b>(sample*i,sample*j)[2];
				newPoint.g = 0;//frgd.at<cv::Vec4b>(sample*i,sample*j)[1];
				newPoint.b = 0;//frgd.at<cv::Vec4b>(sample*i,sample*j)[0];
				outputPointcloud->points.push_back(newPoint);

			}
			/*else
			{
				newPoint.z = std::numeric_limits<float>::quiet_NaN();
				newPoint.x = std::numeric_limits<float>::quiet_NaN();
				newPoint.y = std::numeric_limits<float>::quiet_NaN();
				newPoint.r = std::numeric_limits<unsigned char>::quiet_NaN();
				newPoint.g = std::numeric_limits<unsigned char>::quiet_NaN();
				newPoint.b = std::numeric_limits<unsigned char>::quiet_NaN();
				//outputPointcloud.push_back(newPoint);
			}*/
		}
	}
	/*const std::string file = "test_pcdf%03d.pcd";
    char buf[FILENAME_MAX];
    sprintf(buf, file.c_str(), frame);
    std::string filename(buf);*/
	//pcl::io::savePCDFileASCII (filename, outputPointcloud);
	//std::cout << "Saved " << outputPointcloud.points.size () << " data points to test_pcd.pcd." << std::endl;
		
	return outputPointcloud;
}

template<class DataTypes, class DepthTypes>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr GenerateSyntheticData<DataTypes, DepthTypes>::PCDContourFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage, cv::Mat& distImage, cv::Mat& dotImage)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	//pcl::PointCloud<pcl::PointXYZRGB> pointcloud;
	
	cv::Mat frgd;
	cv::Mat distimg,dotimg;
	
	//cv::pyrDown(rgbImage, frgd, cv::Size(foreground.cols/2, foreground.rows/2));
	distimg = distImage;
	dotimg = dotImage;
	
	targetBorder.resize(0);
	targetWeights.resize(0);
	
	//cv::pyrDown(rgbImage, frgd, cv::Size(foreground.cols/2, foreground.rows/2));
	frgd = rgbImage;	
	Eigen::Matrix3f rgbIntrinsicMatrix;
	rgbIntrinsicMatrix(0,0) = 300.34;
	rgbIntrinsicMatrix(1,1) = 300.34;
	rgbIntrinsicMatrix(0,2) = 157.25;
	rgbIntrinsicMatrix(1,2) = 117.75;
		
    //cv::imwrite("dist.png",distImage);
    //cv::imwrite("RGB.png",frgd);

	float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
	float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy
	pcl::PointXYZRGB newPoint;
	int sample = 4;
	int offsetx = 3;
	ntargetcontours = 0;
	int jj = 0;
	double totalweights = 0;
	for (int i=0;i<(int)depthImage.rows/sample;i++)
	{
		for (int j=0;j<(int)(depthImage.cols-offsetx)/sample;j++)
		{
			float depthValue = (float)depthImage.at<float>(sample*(i),sample*(j+offsetx));
			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)frgd.at<Vec4b>(sample*i,sample*j)[3];
			int bvalue = (int)distimg.at<uchar>(sample*i,sample*(j));
			int dvalue = (int)dotimg.at<uchar>(sample*i,sample*(j));
			
			if (dvalue == 0 && depthValue>0)                // if depthValue is not NaN
			{
				// Find 3D position respect to rgb frame:
				newPoint.z = depthValue;
				newPoint.x = (sample*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
				newPoint.y = (sample*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
				newPoint.r = frgd.at<cv::Vec4b>(sample*i,sample*j)[2];
				newPoint.g = frgd.at<cv::Vec4b>(sample*i,sample*j)[1];
				newPoint.b = frgd.at<cv::Vec4b>(sample*i,sample*j)[0];
				outputPointcloud->points.push_back(newPoint);
				
				//targetWeights.push_back((double)1./(0.12*(1.0+sqrt(bvalue))));
				
				targetWeights.push_back((double)exp(-bvalue/10.0));

				//std::cout << " weight " << (double)1./0.12*(1.0+sqrt(bvalue)) << std::endl;
				
				//targetWeights.push_back((double)0.4/(1.0+bvalue*bvalue));
				
				/*if (avalue > 0 && bvalue < 6)
				{
				targetWeights.push_back((double)3);
				}
				else targetWeights.push_back((double)3);*/
				
				//totalweights += (double)1./(0.12*(1.0+sqrt(bvalue)));
				totalweights += targetWeights[jj];
				
				jj++;
				
				if (avalue > 0 && bvalue < 6) {targetBorder.push_back(true);ntargetcontours++;}
					else targetBorder.push_back(false);
			}
			/*else
			{
				newPoint.z = std::numeric_limits<float>::quiet_NaN();
				newPoint.x = std::numeric_limits<float>::quiet_NaN();
				newPoint.y = std::numeric_limits<float>::quiet_NaN();
				newPoint.r = std::numeric_limits<unsigned char>::quiet_NaN();
				newPoint.g = std::numeric_limits<unsigned char>::quiet_NaN();
				newPoint.b = std::numeric_limits<unsigned char>::quiet_NaN();
				//outputPointcloud.push_back(newPoint);
			}*/
		}
	}
	
	for (int i=0; i < targetWeights.size();i++)
	{
		targetWeights[i]*=((double)targetWeights.size()/totalweights);
		//std::cout << " weights " << (double)targetWeights[i] << std::endl;
	}
		
	/*const std::string file = "test_pcdf%03d.pcd";
    char buf[FILENAME_MAX];
    sprintf(buf, file.c_str(), frame);
    std::string filename(buf);*/
	//pcl::io::savePCDFileASCII (filename, outputPointcloud);
	//std::cout << "Saved " << outputPointcloud.points.size () << " data points to test_pcd.pcd." << std::endl;
		
	return outputPointcloud;
}


template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::ContourFromRGBSynth(cv::Mat& rgbImage, cv::Mat& distImage, cv::Mat& dotImage)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	//pcl::PointCloud<pcl::PointXYZRGB> pointcloud;
	
	cv::Mat frgd;
	cv::Mat distimg,dotimg;
	
	//cv::pyrDown(rgbImage, frgd, cv::Size(foreground.cols/2, foreground.rows/2));
	distimg = distImage;
	dotimg = dotImage;
	frgd = rgbImage;
	
	targetBorder.resize(0);
	targetWeights.resize(0);
	
	const VecCoord& targetp = targetPositions.getValue();
 	
	Eigen::Matrix3f rgbIntrinsicMatrix;
	rgbIntrinsicMatrix(0,0) = 300.34;
	rgbIntrinsicMatrix(1,1) = 300.34;
	rgbIntrinsicMatrix(0,2) = 160;
	rgbIntrinsicMatrix(1,2) = 120;

	int nvisible = 0;
	double totalweights = 0;
	
	//cv::imwrite("dist.png",distimg);
	//cv::imwrite("dot.png",dotimg);
	//cv::imwrite("frgd.png",frgd);

	for (int k = 0; k < targetp.size(); k++)
	{
	int x_u = (int)(targetp[k][0]*rgbIntrinsicMatrix(0,0)/targetp[k][2] + rgbIntrinsicMatrix(0,2));
    int x_v = (int)(targetp[k][1]*rgbIntrinsicMatrix(1,1)/targetp[k][2] + rgbIntrinsicMatrix(1,2));	
	
	std::cout << " depths " << x_u << " " << x_v << std::endl;
	//std::cout << " depths " << (float)depthsN[x_u+(hght-x_v-1)*wdth] << " " <<(float)x[k][2] << std::endl;
	int avalue = (int)frgd.at<Vec4b>(x_v,x_u)[3];
	int bvalue = (int)distimg.at<uchar>(x_v,x_u);
	int dvalue = (int)dotimg.at<uchar>(x_v,x_u);
	
	std::cout << " depths " << x_u << " " << x_v << std::endl;

	if (avalue > 0)                // if depthValue is not NaN
		{
								
				targetWeights.push_back((double)exp(-bvalue/7.0));

				//std::cout << " weight " << (double)1./0.12*(1.0+sqrt(bvalue)) << std::endl;
				
				//targetWeights.push_back((double)0.4/(1.0+bvalue*bvalue));
				
				/*if (avalue > 0 && bvalue < 6)
				{
				targetWeights.push_back((double)3);
				}
				else targetWeights.push_back((double)3);*/
				
				//totalweights += (double)1./(0.12*(1.0+sqrt(bvalue)));
				totalweights += targetWeights[k];
								
				if (bvalue < 6) {targetBorder.push_back(true);ntargetcontours++;}
					else targetBorder.push_back(false);
		}
	
}
	
for (int i=0; i < targetWeights.size();i++)
	{
		targetWeights[i]*=((double)targetWeights.size()/totalweights);
		std::cout << " weights " << (double)targetWeights[i] << std::endl;
	}
		
}


template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::extractTargetPCD()
{
	
	//pcl::PointCloud<pcl::PointXYZRGB>::Ptr target(new pcl::PointCloud<pcl::PointXYZRGB>);
	targetP.reset(new pcl::PointCloud<pcl::PointXYZRGB>); 
 
	//pcl::PointCloud<pcl::PointXYZRGB> target;
	
	targetP = PCDFromRGBD(depth,foreground);
	
	//this->targetBackground.resize(target->size()); targetIgnored.fill(false);  

	//PCDFromRGBD2(depth,foreground,target);

    //pcl::PointCloud<pcl::PointXYZRGB> target0;
	//pcl::io::loadPCDFile ("/home/apetit/soft/catkin_ws/test_pcd12.pcd", target0);
  
 	//pcl::PointCloud<pcl::PointXYZRGB>::Ptr target;
		
	VecCoord targetpos;
	
	if (targetP->size() > 10)
	{
	target.reset(new pcl::PointCloud<pcl::PointXYZRGB>);	
	target = targetP;
	targetpos.resize(target->size());

std::cout << " target size " << target->size() << std::endl;
            Vector3 pos;
            Vector3 col;

	for (unsigned int i=0; i<target->size(); i++)
	{
            pos[0] = (double)target->points[i].x;
            pos[1] = (double)target->points[i].y;
            pos[2] = (double)target->points[i].z;
            targetpos[i]=pos;
			
			/*if(target->points[i].r == 0) 
				{targetBackground[i] = true;
				//std::cout << " ok fore " << std::endl;
				}
				else targetBackground[i] = false;*/
			//if(target->points[i].r == 1) std::cout << " ok back " << std::endl;

			//std::cout << " x " << pos[0] << " y " << pos[1] << " z " << pos[2] << std::endl; 
	} 

	
	/*for (unsigned int i=target->size(); i<targetpos.size(); i++)
	{
            pos[0] = target0[i].x;
            pos[1] = target0[i].y;
            pos[2] = target0[i].z;
            targetpos[i]=pos;
			//std::cout << " x " << pos[0] << " y " << pos[1] << " z " << pos[2] << std::endl;  
	}*/
//targetPositions.setValue(targetpos);
    // build k-d tree
    //const VecCoord&  p = targetPositions.getValue();
    const VecCoord&  p = targetpos;
	targetPositions.setValue(p);
}
    //targetKdTree.build(p);

	//targetKdTree = sourceKdTree;
	
    // detect border
    //if(targetBorder.size()!=p.size()) { targetBorder.resize(p.size()); detectBorder(targetBorder,targetTriangles.getValue()); }

}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::extractTargetPCDContour()
{
	
	//pcl::PointCloud<pcl::PointXYZRGB>::Ptr target(new pcl::PointCloud<pcl::PointXYZRGB>);
	targetP.reset(new pcl::PointCloud<pcl::PointXYZRGB>); 
	
    //cv::imwrite("dotimg.png",seg.dotImage);
 	
	targetP = PCDContourFromRGBD(depth,foreground, seg.distImage,seg.dotImage);
		
	VecCoord targetpos;
	
	if (targetP->size() > 10)
	{
	target.reset(new pcl::PointCloud<pcl::PointXYZRGB>);	
	target = targetP;
	targetpos.resize(target->size());


std::cout << " target size " << target->size() << std::endl;
            Vector3 pos;
            Vector3 col;

	for (unsigned int i=0; i<target->size(); i++)
	{
            pos[0] = (double)target->points[i].x;
            pos[1] = (double)target->points[i].y;
            pos[2] = (double)target->points[i].z;
            targetpos[i]=pos;
			
			/*if(target->points[i].r == 0) 
				{targetBackground[i] = true;
				//std::cout << " ok fore " << std::endl;
				}
				else targetBackground[i] = false;*/
			//if(target->points[i].r == 1) std::cout << " ok back " << std::endl;

			//std::cout << " x " << pos[0] << " y " << pos[1] << " z " << pos[2] << std::endl; 
	}

VecCoord targetContourpos;
targetContourpos.resize(ntargetcontours);
std::cout << " ntargetcontours " << ntargetcontours << std::endl;
int kk = 0;
	for (unsigned int i=0; i<target->size(); i++)
	{
		if (targetBorder[i]){
            pos[0] = (double)target->points[i].x;
            pos[1] = (double)target->points[i].y;
            pos[2] = (double)target->points[i].z;
            targetContourpos[kk]=pos;
			kk++;
		}
			//std::cout << " x " << pos[0] << " y " << pos[1] << " z " << pos[2] << std::endl;  
	}
//targetPositions.setValue(targetpos);
    // build k-d tree
    //const VecCoord&  p = targetPositions.getValue();
	const VecCoord&  p0 = targetpos;
	targetPositions.setValue(p0);
    const VecCoord&  p1 = targetContourpos;
	targetContourPositions.setValue(p1);
	}

}


template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::findCorrespondences (pcl::PointCloud<pcl::PointXYZRGB>::Ptr source, pcl::PointCloud<pcl::PointXYZRGB>::Ptr target, std::vector<int>& correspondences, std::vector<int>& distances)
{
  cout << "correspondence assignment..." << std::flush;
  correspondences.resize (source->size());

  // Use a KdTree to search for the nearest matches in feature space
  pcl::KdTreeFLANN<pcl::PointXYZRGB> descriptor_kdtree;
  descriptor_kdtree.setInputCloud (target);

  // Find the index of the best match for each keypoint, and store it in "correspondences_out"
  const int k = 1;
  std::vector<int> k_indices (k);
  std::vector<float> k_squared_distances (k);
  for (int i = 0; i < static_cast<int> (source->size ()); ++i)
  {
    descriptor_kdtree.nearestKSearch (*source, i, k, k_indices, k_squared_distances);
    correspondences[i] = k_indices[0];
	distances[i] = k_squared_distances[0];
  }
  cout << "OK" << endl;
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::filterCorrespondences ()
{
	
	const VecCoord& x = *this->mstate->getX();
    const VecCoord&  tp = targetPositions.getValue();
	
	std::cout << " tp size 00 " << tp.size() << std::endl;

    unsigned int nbs=x.size(),nbt=tp.size();
	
  cout << "correspondence rejection..." << std::flush;
  std::vector<std::pair<unsigned, unsigned> > correspondences;
  double dist, distmax;
  int cIdxmax;
  /*for (unsigned cIdx = 0; cIdx < source2target_.size (); ++cIdx)
  {
	  distmax = source2target_distances_[cIdx];
	  cIdxmax = source2target_[cIdx];
	    for (unsigned c2Idx = 0; c2Idx < target2source_.size (); ++c2Idx)
		{
			if(target2source_[c2Idx] = cIdx && target2source_distances_[cIdx] > distmax)
			{
				distmax = target2source_distances_[cIdx];
				cIdxmax = c2Idx;
			}
	correspondences.push_back(std::make_pair(cIdx, cIdxmax));
		}
		
  }*/

  	this->sourceIgnored.resize(nbs); sourceIgnored.fill(false);
    this->targetIgnored.resize(nbt); targetIgnored.fill(false);  
  
      //if(blendingFactor.getValue()<1) {
 
    //unsigned int count=0;
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
for (unsigned cIdx = 0; cIdx < source2target_.size (); ++cIdx)
{
	correspondences.push_back(std::make_pair(cIdx, source2target_[cIdx]));
	distances_.push_back(source2target_distances_[cIdx]);
		//closestTarget[cIdx].begin()->first = source2target_distances_[cIdx];
		//closestTarget[cIdx].begin()->second = source2target_[cIdx];
	
}
	
for (unsigned c2Idx = 0; c2Idx < target2source_.size (); ++c2Idx)
{
	int inds = target2source_[c2Idx];
	//closestSource[c2Idx].begin()->first = target2source_distances_[c2Idx];
	//closestSource[c2Idx].begin()->second = target2source_[c2Idx];
	  for (unsigned cIdx = 0; cIdx < correspondences.size(); ++cIdx)
	if(inds == correspondences[cIdx].first && c2Idx != correspondences[cIdx].second)
	{
	correspondences.push_back(std::make_pair(target2source_[c2Idx],c2Idx));
	distances_.push_back(target2source_distances_[c2Idx]);
	targetIgnored[cIdx]=false;
	//correspondences[cIdx].second = c2Idx;
	}
	else targetIgnored[cIdx]=true;
	  
}
  
    //if (target2source_[source2target_[cIdx]] == static_cast<int> (cIdx))
      //correspondences.push_back(std::make_pair(cIdx, source2target_[cIdx]));
	pcl::CorrespondencesPtr corr(new pcl::Correspondences);
	correspondences_ = corr;

  correspondences_->resize (correspondences.size());
  
  for (unsigned cIdx = 0; cIdx < correspondences.size(); ++cIdx)
  {
    (*correspondences_)[cIdx].index_query = correspondences[cIdx].first;
    (*correspondences_)[cIdx].index_match = correspondences[cIdx].second;
  }
//}
  
      // prune outliers
    if(outlierThreshold.getValue()!=0) {
        Real mean=0,stdev=0,count=0;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) {count++; stdev+=closestSource[i].begin()->first; mean+=(Real)(closestSource[i].begin()->first); }
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) {count++; stdev+=closestTarget[i].begin()->first; mean+=(Real)(closestTarget[i].begin()->first); }
        mean=mean/count; stdev=(Real)sqrt(stdev/count-mean*mean);
        mean+=stdev*outlierThreshold.getValue();
        mean*=mean;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(closestSource[i].begin()->first>mean) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(closestTarget[i].begin()->first>mean) targetIgnored[i]=true;
    }
    if(rejectBorders.getValue()) {
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(targetBorder[closestSource[i].begin()->second]) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(sourceBorder[closestTarget[i].begin()->second]) targetIgnored[i]=true;
    }

  /*pcl::registration::CorrespondenceRejectorSampleConsensus<pcl::PointXYZI> rejector;
  rejector.setInputCloud(source_keypoints_);
  rejector.setTargetCloud(target_keypoints_);
  rejector.setInputCorrespondences(correspondences_);
  rejector.getCorrespondences(*correspondences_);*/
  cout << "OK" << endl;
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::updateClosestPoints()
{
	const VecCoord& x = *this->mstate->getX();
	
	//const VecCoord& x = sourceVisiblePositions;
    const VecCoord&  tp = targetPositions.getValue();
	//const VecCoord&  tcp = targetContourPositions.getValue();
		int t = (int)this->getContext()->getTime();


	std::cout << " tp size 00 " << tp.size() << " x size " << x.size() << std::endl;

    unsigned int nbs=x.size(), nbt=tp.size();//, nbtc = tcp.size();

    distanceSet emptyset;
    if(nbs!=closestSource.size()) {initSource();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}
	
	/*if(nbtc!=closestSourceContour.size()) {initSource();  closestSourceContour.resize(nbtc);	
	closestSourceContour.fill(emptyset); 
	cacheDist.resize(nbtc); 
	cacheDist.fill((Real)0.); 
	cacheDist2.resize(nbtc); 
	cacheDist2.fill((Real)0.); 
	previousX.assign(x.begin(),x.end());}*/

	if(nbt!=closestTarget.size()) {initTarget();  closestTarget.resize(nbt);	closestTarget.fill(emptyset);}
	
    //if(nbt!=closestTarget.size()) {extractTargetPCD() ; closestTarget.resize(nbt);	closestTarget.fill(emptyset);}


    //if(nbs==0 || nbt==0) return;

    // closest target points from source points
    if(blendingFactor.getValue()<1) {

    //unsigned int count=0;
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbs;i++)
        {
            /*Real dx=(previousX[i]-x[i]).norm();
            //  closest point caching [cf. Simon96 thesis]
            if(dx>=cacheDist[i] || closestSource[i].size()==0)
            {
                targetKdTree.getNClosest(closestSource[i],x[i],this->cacheSize.getValue() );
                typename distanceSet::iterator it0=closestSource[i].begin(), it1=it0; it1++;
                typename distanceSet::reverse_iterator itn=closestSource[i].rbegin();
                cacheDist[i] =((itn->first)-(it0->first))*(Real)0.5;
                cacheDist2[i]=((it1->first)-(it0->first))*(Real)0.5;
                previousX[i]=x[i];
            }
            else if(dx>=cacheDist2[i]) // in the cache -> update N-1 distances 
            {
                targetKdTree.updateCachedDistances(closestSource[i],x[i]);
                //count++;
            }*/
			targetKdTree.getNClosest(closestSource[i],x[i],1);
						
        }
    //std::cout<<(Real)count*(Real)100./(Real)nbs<<" % cached"<<std::endl;
    }
	
	//getchar();
	
/*

    //unsigned int count=0;
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbs;i++)
        {

			if(sourceBorder[i]){
				
            //else if(dx>=cacheDist2[i]) // in the cache -> update N-1 distances 
            {
                //targetContourKdTree.updateCachedDistances(closestSourceContour[i],x[i]);
                //count++;
				//targetContourKdTree.getNClosest(closestSourceContour[i],x[i],1);
				
            }
			
			}
        }*/

		//std::cout << " tp size 2 " << tp.size() << std::endl;
		//std::cout << " tp size 3 " << targetBackground.size() << std::endl;
		
    // closest source points from target points
    if(blendingFactor.getValue()>0)
    {
        initSource();
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbt;i++)
		{
			//if(!targetBackground[i])
			{
            sourceKdTree.getNClosest(closestTarget[i],tp[i],1);
			}
			
		}
    }

    this->sourceIgnored.resize(nbs); sourceIgnored.fill(false);
    this->targetIgnored.resize(nbt); targetIgnored.fill(false);

    // prune outliers
    if(outlierThreshold.getValue()!=0) {
        Real mean=0,stdev=0,count=0;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) {count++; stdev+=closestSource[i].begin()->first; 
		mean+=(Real)(closestSource[i].begin()->first); 
		}
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) {count++; stdev+=closestTarget[i].begin()->first; mean+=(Real)(closestTarget[i].begin()->first); 
		//std::cout << " distances " << (double)(closestTarget[i].begin()->first) << std::endl;
		}
        mean=mean/count; stdev=(Real)sqrt(stdev/count-mean*mean);
        mean+=stdev*outlierThreshold.getValue();
        mean*=mean;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(closestSource[i].begin()->first>mean ) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(closestTarget[i].begin()->first>mean ) targetIgnored[i]=true;
		
		for(unsigned int i=0;i<nbs;i++) 
			if(closestSource[i].size()) 
				for(unsigned int j=0;j<nbt;j++)  
				if(j == closestSource[i].begin()->second && i == closestTarget[j].begin()->second) 
				{
					//sourceIgnored[i]=true;
					//targetIgnored[j] = true;
				}
			
		//if (t>2)	
		for(unsigned int i=0;i<nbs;i++) 
			if(closestSource[i].size())// && sourceBorder[i]) 
				for(unsigned int j=0;j<nbt;j++) 
				{
					if(i == closestTarget[j].begin()->second){ 
				if(closestSource[i].begin()->first < closestTarget[j].begin()->first)//&& i == closestTarget[j].begin()->second) 
				{
					//sourceIgnored[i]=true;
					//targetIgnored[j] = true;
				}
				//else targetIgnored[j] = true;
					
				}
				}

    }
    if(rejectBorders.getValue()) {
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(targetBorder[closestSource[i].begin()->second]) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(sourceBorder[closestTarget[i].begin()->second]) targetIgnored[i]=true;
    }
    /*if(normalThreshold.getValue()>(Real)-1. && sourceNormals.getValue().size()!=0 && targetNormals.getValue().size()!=0) {
        ReadAccessor< Data< VecCoord > > sn(sourceNormals);
        ReadAccessor< Data< VecCoord > > tn(targetNormals);
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(dot(sn[i],tn[closestSource[i].begin()->second])<normalThreshold.getValue()) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(dot(tn[i],sn[closestTarget[i].begin()->second])<normalThreshold.getValue()) targetIgnored[i]=true;
    }*/
}


template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::updateClosestPointsVisible()
{
	//const VecCoord& x = *this->mstate->getX();
	
	const VecCoord& x = sourceVisiblePositions.getValue();
    const VecCoord&  tp = targetPositions.getValue();
	//const VecCoord&  tcp = targetContourPositions.getValue();

	std::cout << " tp size 00 " << tp.size() << " x size () " << x.size() << std::endl;

    unsigned int nbs=x.size(), nbt=tp.size();//, nbtc = tcp.size();

    distanceSet emptyset;
    if(nbs!=closestSource.size()) {initSourceVisible();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}
	
	/*if(nbtc!=closestSourceContour.size()) {initSource();  closestSourceContour.resize(nbtc);	
	closestSourceContour.fill(emptyset); 
	cacheDist.resize(nbtc); 
	cacheDist.fill((Real)0.); 
	cacheDist2.resize(nbtc); 
	cacheDist2.fill((Real)0.); 
	previousX.assign(x.begin(),x.end());}*/

	if(nbt!=closestTarget.size()) {initTarget();  closestTarget.resize(nbt);	closestTarget.fill(emptyset);}
	
    //if(nbt!=closestTarget.size()) {extractTargetPCD() ; closestTarget.resize(nbt);	closestTarget.fill(emptyset);}


    //if(nbs==0 || nbt==0) return;

    // closest target points from source points
    if(blendingFactor.getValue()<1) {

    //unsigned int count=0;
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbs;i++)
        {
            /*Real dx=(previousX[i]-x[i]).norm();
            //  closest point caching [cf. Simon96 thesis]
            if(dx>=cacheDist[i] || closestSource[i].size()==0)
            {
                targetKdTree.getNClosest(closestSource[i],x[i],this->cacheSize.getValue() );
                typename distanceSet::iterator it0=closestSource[i].begin(), it1=it0; it1++;
                typename distanceSet::reverse_iterator itn=closestSource[i].rbegin();
                cacheDist[i] =((itn->first)-(it0->first))*(Real)0.5;
                cacheDist2[i]=((it1->first)-(it0->first))*(Real)0.5;
                previousX[i]=x[i];
            }
            else if(dx>=cacheDist2[i]) // in the cache -> update N-1 distances 
            {
                targetKdTree.updateCachedDistances(closestSource[i],x[i]);
                //count++;
            }*/
			
			//if(sourceVisible[i])
			targetKdTree.getNClosest(closestSource[i],x[i],1);
			
			
        }
    //std::cout<<(Real)count*(Real)100./(Real)nbs<<" % cached"<<std::endl;
    }
	
		std::cout << " ok target to source " << std::endl;
	
/*

    //unsigned int count=0;
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbs;i++)
        {

			if(sourceBorder[i]){
				
            //else if(dx>=cacheDist2[i]) // in the cache -> update N-1 distances 
            {
                //targetContourKdTree.updateCachedDistances(closestSourceContour[i],x[i]);
                //count++;
				//targetContourKdTree.getNClosest(closestSourceContour[i],x[i],1);
				
            }
			
			}
        }*/

		//std::cout << " tp size 2 " << tp.size() << std::endl;
		//std::cout << " tp size 3 " << targetBackground.size() << std::endl;
		
    // closest source points from target points
    if(blendingFactor.getValue()>0)
    {
        initSourceVisible();
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbt;i++)
		{
			//if(!targetBackground[i])
			{
            sourceKdTree.getNClosest(closestTarget[i],tp[i],1);
			}
			
		}
    }
	
    this->sourceIgnored.resize(nbs); sourceIgnored.fill(false);
    this->targetIgnored.resize(nbt); targetIgnored.fill(false);

    // prune outliers
    if(outlierThreshold.getValue()!=0) {
        Real mean=0,stdev=0,count=0;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) {count++; stdev+=closestSource[i].begin()->first; 
		mean+=(Real)(closestSource[i].begin()->first); 
		}
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) {count++; stdev+=closestTarget[i].begin()->first; mean+=(Real)(closestTarget[i].begin()->first); 
		//std::cout << " distances " << (double)(closestTarget[i].begin()->first) << std::endl;
		}
        mean=mean/count; stdev=(Real)sqrt(stdev/count-mean*mean);
        mean+=stdev*outlierThreshold.getValue();
        mean*=mean;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(closestSource[i].begin()->first>mean ) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(closestTarget[i].begin()->first>mean ) targetIgnored[i]=true;
		
		for(unsigned int i=0;i<nbs;i++) 
			if(closestSource[i].size()) 
				for(unsigned int j=0;j<nbt;j++)  
				if(j == closestSource[i].begin()->second && i == closestTarget[j].begin()->second) 
				{
					//sourceIgnored[i]=true;
					//targetIgnored[j] = true;
				}

    }
    if(rejectBorders.getValue()) {
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(targetBorder[closestSource[i].begin()->second]) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(sourceBorder[closestTarget[i].begin()->second]) targetIgnored[i]=true;
    }
    /*if(normalThreshold.getValue()>(Real)-1. && sourceNormals.getValue().size()!=0 && targetNormals.getValue().size()!=0) {
        ReadAccessor< Data< VecCoord > > sn(sourceNormals);
        ReadAccessor< Data< VecCoord > > tn(targetNormals);
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(dot(sn[i],tn[closestSource[i].begin()->second])<normalThreshold.getValue()) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(dot(tn[i],sn[closestTarget[i].begin()->second])<normalThreshold.getValue()) targetIgnored[i]=true;
    }*/
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::updateClosestPointsContours()
{
	const VecCoord& x = *this->mstate->getX();
    const VecCoord& tp = targetPositions.getValue();
	const VecCoord& xcp = sourceContourPositions.getValue();
	const VecCoord& tcp = targetContourPositions.getValue();

	std::cout << " tp size 00 " << tcp.size() << std::endl;
	
	int t = (int)this->getContext()->getTime();

    unsigned int nbs=x.size(), nbt=tp.size(), nbtc = tcp.size(), nbsc = xcp.size();

    distanceSet emptyset;
    if(nbs!=closestSource.size()) {initSource();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}
	
	/*if(nbtc!=closestSourceContour.size()) {initSource();  closestSourceContour.resize(nbtc);	
	closestSourceContour.fill(emptyset); 
	cacheDist.resize(nbtc); 
	cacheDist.fill((Real)0.); 
	cacheDist2.resize(nbtc); 
	cacheDist2.fill((Real)0.); 
	previousX.assign(x.begin(),x.end());}*/

	if(nbt!=closestTarget.size()) {initTarget();  initTargetContour(); closestTarget.resize(nbt);	closestTarget.fill(emptyset);}
	
	//if(nbsc!=closestSourceContour.size()) {initTargetContour();  closestSourceContour.resize(nbsc);	closestSourceContour.fill(emptyset);}
     
	 //initTargetContour();
     //std::cout << " ok targetcontour 0 " << std::endl;
	
    //if(nbt!=closestTarget.size()) {extractTargetPCD() ; closestTarget.resize(nbt);	closestTarget.fill(emptyset);}


    //if(nbs==0 || nbt==0) return;

    // closest target points from source points
    if(blendingFactor.getValue()<1) {

    //unsigned int count=0;
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbs;i++)
        {
            /*Real dx=(previousX[i]-x[i]).norm();
            //  closest point caching [cf. Simon96 thesis]
            if(dx>=cacheDist[i] || closestSource[i].size()==0)
            {
                targetKdTree.getNClosest(closestSource[i],x[i],this->cacheSize.getValue() );
                typename distanceSet::iterator it0=closestSource[i].begin(), it1=it0; it1++;
                typename distanceSet::reverse_iterator itn=closestSource[i].rbegin();
                cacheDist[i] =((itn->first)-(it0->first))*(Real)0.5;
                cacheDist2[i]=((it1->first)-(it0->first))*(Real)0.5;
                previousX[i]=x[i];
            }
            else if(dx>=cacheDist2[i]) // in the cache -> update N-1 distances 
            {
                targetKdTree.updateCachedDistances(closestSource[i],x[i]);
                //count++;
            }*/

			if(!sourceBorder[i])
			targetKdTree.getNClosest(closestSource[i],x[i],1 );
        }
    //std::cout<<(Real)count*(Real)100./(Real)nbs<<" % cached"<<std::endl;
    }
	
std::cout << " ok targetcontour 1 " << std::endl;

    //unsigned int count=0;
	
	indices.resize(0);
	
	
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        //for(int i=0;i<(int)nbs;i++)
		for(int i=0;i<(int)nbs;i++)
        {

			if(sourceBorder[i])// && t%niterations == 0)
				{
				
            //else if(dx>=cacheDist2[i]) // in the cache -> update N-1 distances 
            {
                //targetContourKdTree.updateCachedDistances(closestSourceContour[i],x[i]);
                //count++;
				
				//targetContourKdTree.getNClosest(closestSourceContour[i],xcp[i],1);
				
				targetContourKdTree.getNClosest(closestSource[i],x[i],1);
				double distmin = 1;
				double dist;
				int kmin;
				for (int k = 0; k < tcp.size(); k++)
				{
					dist = (x[i][0] - tcp[k][0])*(x[i][0] - tcp[k][0]) + (x[i][1] - tcp[k][1])*(x[i][1] - tcp[k][1]) + (x[i][2] - tcp[k][2])*(x[i][2] - tcp[k][2]);
					if (dist < distmin)
					{
						distmin = dist;
						kmin = k;
					}
				}
				indices.push_back(kmin);
				
				unsigned int id=closestSource[i].begin()->second;
				int id1 = indices[i];
				//std::cout << " tcp 0 id " << id << " " << kmin << std::endl;
				
            }
			
			}
        }
		
		std::cout << " ok targetcontour 2 " << std::endl;


		//std::cout << " tp size 2 " << tp.size() << std::endl;
		//std::cout << " tp size 3 " << targetBackground.size() << std::endl;
		
    // closest source points from target points
    if(blendingFactor.getValue()>0)
    {
        initSource();
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbt;i++)
		{
			//if(!targetBackground[i])
			{
            sourceKdTree.getNClosest(closestTarget[i],tp[i],1);
			}
			
		}
    }

    this->sourceIgnored.resize(nbs); sourceIgnored.fill(false);
    this->targetIgnored.resize(nbt); targetIgnored.fill(false);

    // prune outliers
    if(outlierThreshold.getValue()!=0) {
        Real mean=0,stdev=0,count=0;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size() ) {count++; stdev+=closestSource[i].begin()->first; 
		mean+=(Real)(closestSource[i].begin()->first); 
		}
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) {count++; stdev+=closestTarget[i].begin()->first; mean+=(Real)(closestTarget[i].begin()->first); 
		//std::cout << " distances " << (double)(closestTarget[i].begin()->first) << std::endl;
		}
        mean=mean/count; stdev=(Real)sqrt(stdev/count-mean*mean);
        mean+=stdev*outlierThreshold.getValue();
        mean*=mean;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size() ) if(closestSource[i].begin()->first>mean ) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(closestTarget[i].begin()->first>mean ) targetIgnored[i]=true;
		
		for(unsigned int i=0;i<nbs;i++) 
			if(closestSource[i].size() ) 
				for(unsigned int j=0;j<nbt;j++)  
				if(j == closestSource[i].begin()->second && i == closestTarget[j].begin()->second) 
				{
					//sourceIgnored[i]=true;
					//targetIgnored[j] = true;
				}

    }
    /*if(rejectBorders.getValue()) {
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(targetBorder[closestSource[i].begin()->second]) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(sourceBorder[closestTarget[i].begin()->second]) targetIgnored[i]=true;
    }*/
	
    /*if(normalThreshold.getValue()>(Real)-1. && sourceNormals.getValue().size()!=0 && targetNormals.getValue().size()!=0) {
        ReadAccessor< Data< VecCoord > > sn(sourceNormals);
        ReadAccessor< Data< VecCoord > > tn(targetNormals);
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(dot(sn[i],tn[closestSource[i].begin()->second])<normalThreshold.getValue()) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(dot(tn[i],sn[closestTarget[i].begin()->second])<normalThreshold.getValue()) targetIgnored[i]=true;
    }*/
}


template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::updateClosestPointsContoursNormals()
{
	const VecCoord& x = *this->mstate->getX();
    const VecCoord& tp = targetPositions.getValue();
	const VecCoord& xcp = sourceContourPositions.getValue();
	const VecCoord& tcp = targetContourPositions.getValue();
	const VecCoord&  ssn = sourceSurfaceNormalsM.getValue();

	std::cout << " tp size 00 " << tp.size() << " tcp size " << tcp.size() << std::endl;

    unsigned int nbs=x.size(), nbt=tp.size(), nbtc = tcp.size(), nssn = ssn.size();

    distanceSet emptyset;
    if(nbs!=closestSource.size()) {initSource();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}

	if(nbt!=closestTarget.size()) {initTarget();  /*initTargetContour();*/ closestTarget.resize(nbt);	closestTarget.fill(emptyset);}
					//std::cout << " tcp size () " << tcp.size() << std::endl;

    //if(nbt!=closestTarget.size()) {extractTargetPCD() ; closestTarget.resize(nbt);	closestTarget.fill(emptyset);}


    //if(nbs==0 || nbt==0) return;
	
	indices.resize(0);

    // closest target points from source points
    if(blendingFactor.getValue()<1) {

    //unsigned int count=0;
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbs;i++)
        {
				if(sourceBorder[i])
				{

				//targetContourKdTree.getNClosest(closestSource[i],x[i],1);
				double distmin = 1000;
				double dist,dist0;
				double sign;
				int kmin = -1;
				//std::cout << " tcp size () " << tcp.size() << std::endl;
				for (int k = 0; k < tcp.size(); k++)
				{
					//dist = (ssn[i][0] - tcp[k][0])*(ssn[i][0] - tcp[k][0]) + (ssn[i][1] - tcp[k][1])*(ssn[i][1] - tcp[k][1]) +(ssn[i][2] - tcp[k][2])*(ssn[i][2] - tcp[k][2]);
					double norm = sqrt(ssn[i][0]*ssn[i][0] + ssn[i][1]*ssn[i][1] + ssn[i][2]*ssn[i][2]);
					double dist_ = sqrt((ssn[i][0]*(tcp[k][1]-x[i][1])-ssn[i][1]*(tcp[k][0]-x[i][0]))*(ssn[i][0]*(tcp[k][1]-x[i][1])-ssn[i][1]*(tcp[k][0]-x[i][0]))+(ssn[i][2]*(tcp[k][0]-x[i][0])-ssn[i][0]*(tcp[k][2]-x[i][2]))*(ssn[i][2]*(tcp[k][0]-x[i][0])-ssn[i][0]*(tcp[k][2]-x[i][2]))+(ssn[i][1]*(tcp[k][2]-x[i][2])-ssn[i][2]*(tcp[k][1]-x[i][1]))*(ssn[i][1]*(tcp[k][2]-x[i][2])-ssn[i][2]*(tcp[k][1]-x[i][1])));
					dist = dist_/norm;
					sign = (x[i][0] - tcp[k][0])*(ssn[k][0]) + (x[i][1] - tcp[k][1])*(ssn[k][1]) + (x[i][2] - tcp[k][2])*(ssn[k][2]);
					dist0 = (x[i][0] - tcp[k][0])*(x[i][0] - tcp[k][0]) + (x[i][1] - tcp[k][1])*(x[i][1] - tcp[k][1]) + (x[i][2] - tcp[k][2])*(x[i][2] - tcp[k][2]);
					//dist = (x[i][0] - tcp[k][0])*(x[i][0] - tcp[k][0]) + (x[i][1] - tcp[k][1])*(x[i][1] - tcp[k][1]) + (x[i][2] - tcp[k][2])*(x[i][2] - tcp[k][2]);
					//std::cout << " tcp size () " << sqrt(dist0) << std::endl;
					if (sqrt(dist0) < distmin)// && sqrt(dist0) < 0.20)
					{
						distmin = sqrt(dist0);
						kmin = k;
					}
				}
				indices.push_back(kmin);
				
				//std::cout << " i " << i << " kmin " << kmin << std::endl;
				
				unsigned int id=closestSource[i].begin()->second;
				int id1 = indices[i];
				//targetKdTree.getNClosest(closestSource[i],x[i],1);
			}
			
			targetKdTree.getNClosest(closestSource[i],x[i],1);
			
        }
    //std::cout<<(Real)count*(Real)100./(Real)nbs<<" % cached"<<std::endl;
    }
	
		//std::cout << " tp size 2 " << tp.size() << std::endl;
		//std::cout << " tp size 3 " << targetBackground.size() << std::endl;
		
    // closest source points from target points
    if(blendingFactor.getValue()>0)
    {
        initSource();
#ifdef USING_OMP_PRAGMAS
        #pragma omp parallel for
#endif
        for(int i=0;i<(int)nbt;i++)
		{
			{
            sourceKdTree.getNClosest(closestTarget[i],tp[i],1);
			}
			
		}
    }

    this->sourceIgnored.resize(nbs); sourceIgnored.fill(false);
    this->targetIgnored.resize(nbt); targetIgnored.fill(false);

    // prune outliers
    if(outlierThreshold.getValue()!=0) {
        Real mean=0,stdev=0,count=0;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) {count++; stdev+=closestSource[i].begin()->first; 
		mean+=(Real)(closestSource[i].begin()->first); 
		}
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) {count++; stdev+=closestTarget[i].begin()->first; mean+=(Real)(closestTarget[i].begin()->first); 
		//std::cout << " distances " << (double)(closestTarget[i].begin()->first) << std::endl;
		}
        mean=mean/count; stdev=(Real)sqrt(stdev/count-mean*mean);
        mean+=stdev*outlierThreshold.getValue();
        mean*=mean;
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(closestSource[i].begin()->first>mean ) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(closestTarget[i].begin()->first>mean ) targetIgnored[i]=true;
		
		for(unsigned int i=0;i<nbs;i++) 
			if(closestSource[i].size()) 
				for(unsigned int j=0;j<nbt;j++)  
				if(j == closestSource[i].begin()->second && i == closestTarget[j].begin()->second) 
				{
					//sourceIgnored[i]=true;
					//targetIgnored[j] = true;
				}

    }
    if(rejectBorders.getValue()) {
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(targetBorder[closestSource[i].begin()->second]) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(sourceBorder[closestTarget[i].begin()->second]) targetIgnored[i]=true;
    }
    /*if(normalThreshold.getValue()>(Real)-1. && sourceNormals.getValue().size()!=0 && targetNormals.getValue().size()!=0) {
        ReadAccessor< Data< VecCoord > > sn(sourceNormals);
        ReadAccessor< Data< VecCoord > > tn(targetNormals);
        for(unsigned int i=0;i<nbs;i++) if(closestSource[i].size()) if(dot(sn[i],tn[closestSource[i].begin()->second])<normalThreshold.getValue()) sourceIgnored[i]=true;
        for(unsigned int i=0;i<nbt;i++) if(closestTarget[i].size()) if(dot(tn[i],sn[closestTarget[i].begin()->second])<normalThreshold.getValue()) targetIgnored[i]=true;
    }*/
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::updateClosestPointsPCL()
{
	const VecCoord& x = *this->mstate->getX();
    const VecCoord&  tp = targetPositions.getValue();
	
		std::cout << " tp size 00 " << tp.size() << std::endl;

    unsigned int nbs=x.size(),nbt=tp.size();

    distanceSet emptyset;
    if(nbs!=closestSource.size()) {initSource();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}
    if(nbt!=closestTarget.size()) {initTarget();  closestTarget.resize(nbt);	closestTarget.fill(emptyset);}

	
	source.reset(new pcl::PointCloud<pcl::PointXYZRGB>); 
	
	pcl::PointXYZRGB newPoint;
	for (unsigned int i=0; i<nbs; i++)
	{
	newPoint.z = x[i][0];
	newPoint.x = x[i][1];
	newPoint.y = x[i][2];
	newPoint.r = 0;
	newPoint.g = 0;
	newPoint.b = 0;
	source->points.push_back(newPoint);
	} 
	
  findCorrespondences (source, target, source2target_, source2target_distances_);
  findCorrespondences (target, source, target2source_, target2source_distances_);
  filterCorrespondences();
  
  
}

extern "C" {
  int sgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
             const float *alpha, const float *a, const int *lda, const float *b, const int *ldb,
             const float *beta, float *c, const int *ldc);
  int saxpy_(const int *n, const float *sa, const float *sx, const int *incx, float *sy, const int *incy);
  int sgemv_(const char *trans, const int *m, const int *n,
             const float *alpha, const float *a, const int *lda,
             const float *x, const int *incx, const float *beta, float *y, const int *incy);
  float sasum_(const int *n, const float *sx, const int *incx);
  int sscal_(const int *n, const float *sa, float *sx, const int *incx);
}

static void
updateA0(int rowsA, int colsA, int pitchA,
	const float* h_Xx, const float* h_Xy, const float* h_Xz, 
	const float* h_Yx, const float* h_Yy, const float* h_Yz,
	float* h_A,
	float sigma_p2){

#pragma omp parallel for
  for(int c=0; c<colsA; c++){

    float Xx = h_Xx[c];
    float Xy = h_Xy[c];
    float Xz = h_Xz[c];

    for(int r=0; r<rowsA; r++){

      float Yx = h_Yx[r];
      float Yy = h_Yy[r];
      float Yz = h_Yz[r];


      // #define Euclid(a,b,c) ((a)*(a)+(b)*(b)+(c)*(c))
      //     float tmp =
      //       Euclid(Xx - (R(0)*Yx + R(1)*Yy + R(2)*Yz + t(0)),
      //              Xy - (R(3)*Yx + R(4)*Yy + R(5)*Yz + t(1)),
      //              Xz - (R(6)*Yx + R(7)*Yy + R(8)*Yz + t(2)) );
    
      //     tmp = expf(-tmp/sigma_p^2)

      float tmpX = Xx - Yx;
      float tmpY = Xy - Yy;
      float tmpZ = Xz - Yz;

      tmpX *= tmpX;
      tmpY *= tmpY;
      tmpZ *= tmpZ;

      tmpX += tmpY;
      tmpX += tmpZ;

      tmpX /= sigma_p2;
      tmpX = expf(-tmpX);
      h_A[c * pitchA + r] = tmpX;

    }
  }
}


static void
updateA(int rowsA, int colsA, int pitchA,
	const float* h_Xx, const float* h_Xy, const float* h_Xz, 
	const float* h_Yx, const float* h_Yy, const float* h_Yz,
	const float* h_R, const float* h_t,
	float* h_A,
	float sigma_p2){

#pragma omp parallel for
  for(int c=0; c<colsA; c++){

    float Xx = h_Xx[c];
    float Xy = h_Xy[c];
    float Xz = h_Xz[c];

    for(int r=0; r<rowsA; r++){

      float Yx = h_Yx[r];
      float Yy = h_Yy[r];
      float Yz = h_Yz[r];

#define R(i) h_R[i]
#define t(i) h_t[i]

      // #define Euclid(a,b,c) ((a)*(a)+(b)*(b)+(c)*(c))
      //     float tmp =
      //       Euclid(Xx - (R(0)*Yx + R(1)*Yy + R(2)*Yz + t(0)),
      //              Xy - (R(3)*Yx + R(4)*Yy + R(5)*Yz + t(1)),
      //              Xz - (R(6)*Yx + R(7)*Yy + R(8)*Yz + t(2)) );
    
      //     tmp = expf(-tmp/sigma_p^2)

      float tmpX = Xx - (R(0)*Yx + R(1)*Yy + R(2)*Yz + t(0));
      float tmpY = Xy - (R(3)*Yx + R(4)*Yy + R(5)*Yz + t(1));
      float tmpZ = Xz - (R(6)*Yx + R(7)*Yy + R(8)*Yz + t(2));

      tmpX *= tmpX;
      tmpY *= tmpY;
      tmpZ *= tmpZ;

      tmpX += tmpY;
      tmpX += tmpZ;

      tmpX /= sigma_p2;
      tmpX = expf(-tmpX);

      h_A[c * pitchA + r] = tmpX;

    }
  }
}


static void
normalizeRowsOfA(int rowsA, int colsA, int pitchA,
		 float *h_A,
		 const float *h_C){
  
#pragma omp parallel for
  for(int c=0; c<colsA; c++)
    for(int r=0; r<rowsA; r++)
	{
      if(h_C[r] > 10e-7f)
	// each element in A is normalized C, then squre-rooted
	h_A[c * pitchA + r] = sqrtf( h_A[c * pitchA + r] / h_C[r]);
      else
	h_A[c * pitchA + r] = 1.0f/colsA; // ad_hoc code to avoid 0 division
	
	//std::cout << " h_A " << h_A[c * pitchA + r] << std::endl;
	
	}

}

static void
elementwiseDivision(int Xsize,
		    float* h_Xx, float* h_Xy, float* h_Xz,
		    const float* h_lambda){

#pragma omp parallel for
  for(int x=0; x<Xsize; x++){
    float l_lambda = h_lambda[x];
    h_Xx[x] /= l_lambda;
    h_Xy[x] /= l_lambda;
    h_Xz[x] /= l_lambda;
  }

}

static void
elementwiseMultiplication(int Xsize,
			  float* h_Xx, float* h_Xy, float* h_Xz,
			  const float* h_lambda){

#pragma omp parallel for
  for(int x=0; x<Xsize; x++){
    float l_lambda = h_lambda[x];
    h_Xx[x] *= l_lambda;
    h_Xy[x] *= l_lambda;
    h_Xz[x] *= l_lambda;
  }
}



static void
centeringXandY(int rowsA,
	       const float* h_Xc, const float* h_Yc,
	       const float* h_Xx, const float* h_Xy, const float* h_Xz,
	       const float* h_Yx, const float* h_Yy, const float* h_Yz,
	       float* h_XxCenterd, float* h_XyCenterd, float* h_XzCenterd,
	       float* h_YxCenterd, float* h_YyCenterd, float* h_YzCenterd){

  // do for both X and Y at the same time
  
#pragma omp parallel for
  for(int r=0; r<rowsA; r++){
    h_XxCenterd[r] = h_Xx[r] - h_Xc[0];
    h_XyCenterd[r] = h_Xy[r] - h_Xc[1];
    h_XzCenterd[r] = h_Xz[r] - h_Xc[2];
    
    h_YxCenterd[r] = h_Yx[r] - h_Yc[0];
    h_YyCenterd[r] = h_Yy[r] - h_Yc[1];
    h_YzCenterd[r] = h_Yz[r] - h_Yc[2];
  }

}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::cloud2data(const pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,
		float **X, int &Xsize )
{
  
  Xsize = cloud->size();
  
  float* h_X = new float [Xsize * 3];
  float* h_Xx = &h_X[Xsize*0];
  float* h_Xy = &h_X[Xsize*1];
  float* h_Xz = &h_X[Xsize*2];
  for (int i = 0; i < Xsize; i++)
  {
    h_Xx[i] = cloud->points[i].x;
    h_Xy[i] = cloud->points[i].y;
    h_Xz[i] = cloud->points[i].z;
  }
  
  *X = h_X;
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::cloud2dataC(const pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud,
		float **X, int &Xsize )
{
  
  Xsize = cloud->size();
  
  float* h_X = new float [Xsize * 3];
  float* h_Xx = &h_X[Xsize*0];
  float* h_Xy = &h_X[Xsize*1];
  float* h_Xz = &h_X[Xsize*2];
  for (int i = 0; i < Xsize; i++)
  {
    h_Xx[i] = cloud->points[i].x;
    h_Xy[i] = cloud->points[i].y;
    h_Xz[i] = cloud->points[i].z;
  }
  
  *X = h_X;
}

template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::updateCPSoft(const pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_target, 
	   const pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_source)
{
  
  int Xsize, Ysize;
  float *h_X, *h_Y;
  cloud2dataC(cloud_target, &h_X, Xsize);
  cloud2dataC(cloud_source, &h_Y, Ysize);
  
  // const for blas functions
  float onef = 1.0f;
  float zerof = 0.0f;
  int onei = 1;
  int threei = 3;
  char transN = 'n';
  char transT = 't';

  //
  // memory allocation
  //

   const float* h_Xx = &h_X[Xsize*0];
   const float* h_Xy = &h_X[Xsize*1];
   const float* h_Xz = &h_X[Xsize*2];
 
   const float* h_Yx = &h_Y[Ysize*0];
   const float* h_Yy = &h_Y[Ysize*1];
   const float* h_Yz = &h_Y[Ysize*2];


  // NOTE on matrix A
  // number of rows:     Ysize, or rowsA
  // number of columns : Xsize, or colsA
  // 
  //                    [0th in X] [1st]  ... [(Xsize-1)] 
  // [0th point in Y] [ A(0,0)     A(0,1) ... A(0,Xsize-1)      ] 
  // [1st           ] [ A(1,0)     A(1,1) ...                   ]
  // ...              [ ...                                     ]
  // [(Ysize-1)     ] [ A(Ysize-1, 0)     ... A(Ysize-1,Xsize-1)]
  //
  // 
  // CAUTION on matrix A
  // A is allcoated as a column-maijor format for the use of cublas.
  // This means that you must acces an element at row r and column c as:
  // A(r,c) = A[c * pitchA + r]


  int rowsA = Ysize;
  int colsA = Xsize;

  // pitchA: leading dimension of A, which is ideally equal to rowsA,
  //          but actually larger than that.
  //pitchA = (rowsA / 4 + 1) * 4;
  pitchA = rowsA;
  
  std::cout << " rowsA " << rowsA << " colsA " << colsA << " pitchA " << pitchA << std::endl;

  h_A = new float [pitchA*colsA];

  // a vector with all elements of 1.0f
  float* h_one = new float [max(Xsize,Ysize)];
  for(int t = 0; t < max(Xsize,Ysize); t++) h_one[t] = 1.0f;

  float *h_C = new float [rowsA]; // sum of a row in A
  float *h_lambda = new float [rowsA]; // weight of a row in A

  //
  // timer, notimer
  //

//#define START_TIMER(timer)
//#define STOP_TIMER(timer)




  // EM-ICP main loop
  int Titer = 1;
{

    fprintf(stderr, "%d iter. sigma_p2 %f  \n", Titer++, sigma_p2);
    // fprintf(stderr, "time %.10f [s]\n", cutGetTimerValue(timerTotal) / 1000.0f);


      //
      // UpdateA
      //

      //START_TIMER(timerUpdateA);

	  updateA0(rowsA, colsA, pitchA,h_Xx, h_Xy, h_Xz, h_Yx, h_Yy, h_Yz,h_A, sigma_p2);

      //STOP_TIMER(timerUpdateA);



      //
      // Normalization of A
      //

      // cublasSgemv (char trans, int m, int n, float alpha, const float *A, int lda,
      //              const float *x, int incx, float beta, float *y, int incy)
      // int sgemv_(char *trans, int *m, int *n,
      // 	     float *alpha, float *a, int *lda,
      // 	     float *x, int *incx, float *beta, float *y, int *incy);
      //    y = alpha * op(A) * x + beta * y,

      // A * one vector = vector with elements of row-wise sum
      //     h_A      *    h_one    =>  h_C
      //(rowsA*colsA) *  (colsA*1)  =  (rowsA*1)
      sgemv_(&transN,          // char trans
	     &rowsA, &colsA, // int m (rows of A), n (cols of A) ; not op(A)
	     &onef,         // float alpha
	     h_A, &pitchA,  // const float *A, int lda
	     h_one, &onei,     // const float *x, int incx
	     &zerof,         // float beta
	     h_C, &onei);      // float *y, int incy


      // void cublasSaxpy (int n, float alpha, const float *x, int incx, float *y, int incy)
      // int saxpy_(int *n, float *sa, float *sx, int *incx, float *sy, int *incy);
      // alpha * x + y => y
      // exp(-d_0^2/sigma_p2) * h_one + h_C => h_C
      {
	float alpha = expf(-d_02/sigma_p2);
	saxpy_(&rowsA, &alpha, h_one, &onei, h_C, &onei);
      }


      normalizeRowsOfA
	(rowsA, colsA, pitchA, h_A, h_C);
      
      
      sigma_p2 *= sigma_factor;
  }
 

  delete [] h_A;

  delete [] h_one;

  delete [] h_C;
  delete [] h_lambda;


}


template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::updateClosestPointsSoft()
{
	
    cout << "EM-ICP paramters" << endl
	 << "sigma_p2 " << paramSoft.sigma_p2 << endl
	 << "sigma_inf " << paramSoft.sigma_inf << endl
	 << "sigma_factor " << paramSoft.sigma_factor << endl
	 << "d_02 " << paramSoft.d_02 << endl;

	
	const VecCoord& x = *this->mstate->getX();
    const VecCoord&  tp = targetPositions.getValue();
	
	unsigned int nbs=x.size(), nbt=tp.size();//, nbtc = tcp.size();

    distanceSet emptyset;
    if(nbs!=closestSource.size()) {initSource();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}
	
	/*if(nbtc!=closestSourceContour.size()) {initSource();  closestSourceContour.resize(nbtc);	
	closestSourceContour.fill(emptyset); 
	cacheDist.resize(nbtc); 
	cacheDist.fill((Real)0.); 
	cacheDist2.resize(nbtc); 
	cacheDist2.fill((Real)0.); 
	previousX.assign(x.begin(),x.end());}*/

	if(nbt!=closestTarget.size()) {initTarget();  closestTarget.resize(nbt);	closestTarget.fill(emptyset);}
	

	source.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
	source_registered.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
	
	pcl::PointXYZRGB newPoint;
	for (unsigned int i=0; i<nbs; i++)
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
	
	
	updateCPSoft(source, target);

    this->sourceIgnored.resize(nbs); sourceIgnored.fill(false);
    this->targetIgnored.resize(nbt); targetIgnored.fill(false);

    // prune outliers
}

template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::determineRigidTransformation ()
{
	
	const VecCoord& x = *this->mstate->getX();
    const VecCoord&  tp = targetPositions.getValue();
	
    unsigned int nbs=x.size(),nbt=tp.size();
	
    distanceSet emptyset;
    if(nbs!=closestSource.size()) {initSource();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}
    if(nbt!=closestTarget.size()) {initTarget();  closestTarget.resize(nbt);	closestTarget.fill(emptyset);}

	
	source.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
	source_registered.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
	
	pcl::PointXYZRGB newPoint;
	for (unsigned int i=0; i<nbs; i++)
	{
	newPoint.z = x[i][2];
	newPoint.x = x[i][0];
	newPoint.y = x[i][1];
	newPoint.r = 0;
	newPoint.g = 0;
	newPoint.b = 0;
	source->points.push_back(newPoint);
	} 
	
	std::cout << " tp size 00 " << target->size() << std::endl;

  cout << "final registration..." << std::flush;
  pcl::Registration<pcl::PointXYZRGB, pcl::PointXYZRGB>::Ptr registration (new pcl::IterativeClosestPoint<pcl::PointXYZRGB, pcl::PointXYZRGB>);
  registration->setInputCloud(target);
  //registration->setInputCloud(source_segmented_);
  registration->setInputTarget (source);
  registration->setMaxCorrespondenceDistance(0.2);
  registration->setRANSACOutlierRejectionThreshold (0.1);
  registration->setTransformationEpsilon (0.000001);
  registration->setMaximumIterations (1000);
  /*registration->setMaxCorrespondenceDistance(0.1);
  registration->setRANSACOutlierRejectionThreshold (0.05);
  registration->setTransformationEpsilon (0.0001);
  registration->setMaximumIterations (20);*/
  	std::cout << " ok registration " << nbs << " " << nbt << std::endl;

  registration->align(*source_registered);
  	std::cout << " ok registration " << nbs << " " << nbt << std::endl;

  Eigen::Matrix4f transformation_matrix1 = registration->getFinalTransformation();
  transformation_matrix = transformation_matrix1.inverse();
  
std::cout << " ok registration " << std::endl;

Eigen::Matrix3f mat;
  for (int i = 0; i < 3; i++)
	  for (int j = 0; j< 3;j ++)
       mat(i,j) = transformation_matrix(i,j);
	   	   
Eigen::Vector3f ea = mat.eulerAngles(0,1,2);
Eigen::Quaternionf q(mat);

pcl::transformPointCloud(*source, *source_registered, transformation_matrix);

helper::WriteAccessor< Data<VecCoord> > x1 = *this->mstate->write(core::VecCoordId::position());
for (unsigned int i=0; i<nbs; i++)
{
	newPoint = source_registered->points[i];
	x1[i][0] = newPoint.x;
	x1[i][1] = newPoint.y;
	x1[i][2] = newPoint.z;
}

  
  	//for (unsigned int i=0; i<nbs; i++)
	{
		//this->mstate->applyTranslation(transformation_matrix(0,3),transformation_matrix(1,3),transformation_matrix(2,3));
		
		//this->mstate->applyTranslation(-0.01,0.0,0.0);
		//this->mstate->applyRotation(ea(0),ea(1),ea(2));

		/*(*this->mstate->getX())[i][0] = 0;
		(*this->mstate->getX())[i][1] = 0;
		(*this->mstate->getX())[i][2] = 0;*/

	}
  cout << "OK" <<  transformation_matrix(0,3) << " " << transformation_matrix(1,3) << " " << transformation_matrix(2,3) << endl;
}

template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::determineRigidTransformationVisible ()
{
	const VecCoord& x0 = *this->mstate->getX();
	const VecCoord& x = sourceVisiblePositions.getValue();
    const VecCoord&  tp = targetPositions.getValue();
	
    unsigned int nbs=x.size(),nbt=tp.size(), nbs0 = x0.size();

    distanceSet emptyset;
    //if(nbs!=closestSource.size()) {initSourceVisible();  closestSource.resize(nbs);	closestSource.fill(emptyset); cacheDist.resize(nbs); cacheDist.fill((Real)0.); cacheDist2.resize(nbs); cacheDist2.fill((Real)0.); previousX.assign(x.begin(),x.end());}
    //if(nbt!=closestTarget.size()) {initTarget();  closestTarget.resize(nbt);	closestTarget.fill(emptyset);}

	
	source.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
	source_registered.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
	
	pcl::PointXYZRGB newPoint;
	for (unsigned int i=0; i<nbs; i++)
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

   source0.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
   source_registered0.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
 
	
	for (unsigned int i=0; i<nbs0; i++)
	{
	newPoint.z = x0[i][2];
	newPoint.x = x0[i][0];
	newPoint.y = x0[i][1];
	//cout << "OK " <<  newPoint.z << " " << newPoint.x << " " << newPoint.y  << endl;
	newPoint.r = 0;
	newPoint.g = 0;
	newPoint.b = 0;
	source0->points.push_back(newPoint);
	} 
	
	std::cout << " tp size 00 " << target->size() << std::endl;

  cout << "final registration..." << std::flush;
  pcl::Registration<pcl::PointXYZRGB, pcl::PointXYZRGB>::Ptr registration (new pcl::IterativeClosestPoint<pcl::PointXYZRGB, pcl::PointXYZRGB>);
  registration->setInputCloud(target);
  //registration->setInputCloud(source_segmented_);
  registration->setInputTarget (source);
  registration->setMaxCorrespondenceDistance(0.2);
  registration->setRANSACOutlierRejectionThreshold (0.1);
  registration->setTransformationEpsilon (0.000001);
  registration->setMaximumIterations (1000);
  /*registration->setMaxCorrespondenceDistance(0.1);
  registration->setRANSACOutlierRejectionThreshold (0.05);
  registration->setTransformationEpsilon (0.0001);
  registration->setMaximumIterations (20);*/
  	std::cout << " ok registration " << nbs << " " << nbt << std::endl;

  registration->align(*source_registered);
  	std::cout << " ok registration " << nbs << " " << nbt << std::endl;

  Eigen::Matrix4f transformation_matrix1 = registration->getFinalTransformation();
  transformation_matrix = transformation_matrix1.inverse();
  
std::cout << " ok registration " << std::endl;

Eigen::Matrix3f mat;
  for (int i = 0; i < 3; i++)
	  for (int j = 0; j< 3;j ++)
       mat(i,j) = transformation_matrix(i,j);
	   	   
Eigen::Vector3f ea = mat.eulerAngles(0,1,2);
Eigen::Quaternionf q(mat);

pcl::transformPointCloud(*source0, *source_registered0, transformation_matrix);

helper::WriteAccessor< Data<VecCoord> > x1 = *this->mstate->write(core::VecCoordId::position());
for (unsigned int i=0; i<nbs0; i++)
{
	newPoint = source_registered0->points[i];
	x1[i][0] = newPoint.x;
	x1[i][1] = newPoint.y;
	x1[i][2] = newPoint.z;

}

  
  	//for (unsigned int i=0; i<nbs; i++)
	{
		//this->mstate->applyTranslation(transformation_matrix(0,3),transformation_matrix(1,3),transformation_matrix(2,3));
		
		//this->mstate->applyTranslation(-0.01,0.0,0.0);
		//this->mstate->applyRotation(ea(0),ea(1),ea(2));

		/*(*this->mstate->getX())[i][0] = 0;
		(*this->mstate->getX())[i][1] = 0;
		(*this->mstate->getX())[i][2] = 0;*/

	}
 // cout << "OK" <<  transformation_matrix(0,3) << " " << transformation_matrix(1,3) << " " << transformation_matrix(2,3) << endl;
}


template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::addForceSoft(const core::MechanicalParams* /*mparams*/,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
	
	int t = (int)this->getContext()->getTime();
	if (useSensor.getValue())
	 getImages();
	 else readImages();
	 
	 //timeOverall
	 double timeT;

	if (t == 0)
	{
        timeFile.open ("times8.txt");
		timeFile.clear();
		initSegmentation();
		extractTargetPCD();
		timeOverall = 0;
		timeTotal = 0;
		timeT = (double)getTickCount();
		//extractTargetContour();
	}
	else if (t > 0 && t%niterations == 0)
	{
    std::cout << " ok init seg " << std::endl;
	double time0 = (double)getTickCount();
	//getchar();
		segment();
		extractTargetPCD();
		//extractTargetPCDContour();
		//extractTargetContour();
			time0 = ((double)getTickCount() - time0)/getTickFrequency();
		cout << "Time segment + extract PCD " << time0 << endl;
		timeOverall += time0;
		double time1 = (double)getTickCount();
			if (t < 10) 
			determineRigidTransformation();
			time1 = ((double)getTickCount() - time1)/getTickFrequency();
			cout << "Rigid ICP " << time1 << endl;
			timeOverall += time1;
	}
	
	//extractSourceContour();
	if (t < 2)
	setViewPoint();
	
	
	//if (t > 0 && t == niterations + 1) getchar();
	
	if (t > 0 && t%niterations == 0)
	{
	double timertt = (double)getTickCount();
	rtt = new cv::Mat;
	cv::Mat rtt_;
    renderToTexture(rtt_);
	*rtt = rtt_;
	listrtt.push_back(rtt);
	timertt = ((double)getTickCount() - timertt)/getTickFrequency();
    cout << "Time RTT " << timertt << endl;
	}
		
	//if (t == nimages*niterations-1)
	
	std::cout << " time " << t << " " << nimages*niterations-6*niterations - 1 << std::endl;
	
	if (t == nimages*niterations-4*niterations - 1)
	{
		//getchar();
		writeImages();
	}
      int key=cvWaitKey(3);
      if ((char)key == 27) return;


    if(ks.getValue()==0) return;

    VecDeriv&        f = *_f.beginEdit();           //WDataRefVecDeriv f(_f);
    const VecCoord&  x = _x.getValue();			//RDataRefVecCoord x(_x);
    const VecDeriv&  v = _v.getValue();			//RDataRefVecDeriv v(_v);
    ReadAccessor< Data< VecCoord > > tn(targetNormals);
    ReadAccessor< Data< VecCoord > > tp(targetPositions);
	ReadAccessor< Data< VecCoord > > tcp(targetContourPositions);
	
	tpos.resize(tp.size());
	
    const vector<Spring>& s = this->springs.getValue();
    this->dfdx.resize(s.size());
    this->closestPos.resize(s.size());
	
	double time = (double)getTickCount();

	if(!pcl)
	{
		updateClosestPoints();
		//updateClosestPointsSoft();
	}
	else updateClosestPointsPCL();
	
	
	//if (t == 100*niterations +5) getchar();
	
    m_potentialEnergy = 0;

    // get attraction/ projection factors
    Real attrF=0.99;//(Real) blendingFactor.getValue();
    if(attrF<(Real)0.) attrF=(Real)0.;
    if(attrF>(Real)1.) attrF=(Real)1.;
    Real projF=((Real)1.-attrF);

    if(tp.size()==0)
        for (unsigned int i=0; i<s.size(); i++)
            closestPos[i]=x[i];
    else {
		//std::cout << " ok 2 " << tp.size() << std::endl;

        // count number of attractors
        cnt.resize(s.size()); cnt.fill(0);  
		if(attrF>0) 
			for (unsigned int i=0; i<tp.size(); i++) 
			{
				tpos[i] = tp[i];
				
				if(!targetIgnored[i])// && !targetBackground[i])
				{ 
					cnt[closestTarget[i].begin()->second]++;
				}
			}

        if(theCloserTheStiffer.getValue())
        {
            // find the min and the max distance value from source point to target point
            min=0;
            max=0;
            for (unsigned int i=0; i<x.size(); i++)
            {
                if(min==0 || min>closestSource[i].begin()->first) min=closestSource[i].begin()->first;
                if(max==0 || max<closestSource[i].begin()->first) max=closestSource[i].begin()->first;
            }
        }

        // compute targetpos = projF*closestto + attrF* sum closestfrom / count
		
		std::cout << " source size " << s.size() << " target size " << tp.size() << std::endl;
		
        // projection to point or plane
        if(projF>0) {
			//#pragma omp parallel for
            for (int i=0; i<s.size(); i++)
			{
				
				if(!sourceIgnored[i]) {
                    unsigned int id=closestSource[i].begin()->second;
                    if(projectToPlane.getValue() && tn.size()!=0)	closestPos[i]=(x[i]+tn[id]*dot(tp[id]-x[i],tn[id]))*projF;
                    else closestPos[i]=tp[id]*projF;
                    if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                }
                else {
                    closestPos[i]=x[i]*projF;
                    if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                }			
			}
        }
        else for (unsigned int i=0; i<s.size(); i++) { if(!cnt[i]) closestPos[i]=x[i]; else closestPos[i].fill(0); }

        // attraction

        if(attrF>0)
            for (unsigned int i=0; i<tp.size(); i++)
			 if(!targetIgnored[i] ) //&& !targetBackground[i])	
				{
                    unsigned int id=closestTarget[i].begin()->second;
					/*if(t>20*niterations){
					if(targetWeights[i] == 0.5)
                    closestPos[id]+=tp[i]*attrF/(Real)cnt[id];
					else closestPos[id]+= tp[i]*(0.2*attrF)/(Real)cnt[id];
                    sourceIgnored[id]=false;
					}
					else*/
					{                    
					closestPos[id]+=tp[i]*attrF/(Real)cnt[id];
                    sourceIgnored[id]=false;
					}
                }
                /*if(!targetIgnored[i] ) //&& !targetBackground[i])	
				{
			for (int j=0; j<s.size(); j++)
				{
                //if(!sourceIgnored[i]) {
				Real weight = (Real)h_A[i * pitchA + j];
				//std::cout << weight << std::endl;
				//closestPos[j] += ((Real)1/tp.size()) * tp[i];
                    if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                }
                }*/
}

h_A = new float [tp.size()*s.size()];
float *h_C = new float [s.size()];

#pragma omp parallel for
  for(unsigned int i=0; i<tp.size(); i++){

    float Xx = tp[i][0];
    float Xy = tp[i][1];
    float Xz = tp[i][2];

    for(unsigned int j=0; j<s.size(); j++){

      float Yx = x[j][0];
      float Yy = x[j][1];
      float Yz = x[j][2];


      // #define Euclid(a,b,c) ((a)*(a)+(b)*(b)+(c)*(c))
      //     float tmp =
      //       Euclid(Xx - (R(0)*Yx + R(1)*Yy + R(2)*Yz + t(0)),
      //              Xy - (R(3)*Yx + R(4)*Yy + R(5)*Yz + t(1)),
      //              Xz - (R(6)*Yx + R(7)*Yy + R(8)*Yz + t(2)) );
    
      //     tmp = expf(-tmp/sigma_p^2)

      float tmpX = Xx - Yx;
      float tmpY = Xy - Yy;
      float tmpZ = Xz - Yz;

      tmpX *= tmpX;
      tmpY *= tmpY;
      tmpZ *= tmpZ;

      tmpX += tmpY;
      tmpX += tmpZ;

      tmpX /= sigma_p2;
      tmpX = expf(-tmpX);
      h_A[i * s.size() + j] = tmpX;

    }
}

#pragma omp parallel for
for(int r=0; r<s.size(); r++)
  for(int c=0; c<tp.size(); c++)
	{
		h_C[r] += h_A[c * s.size() + r];
	}


#pragma omp parallel for
for(int c=0; c<tp.size(); c++)
 for(int r=0; r<s.size(); r++)
	{
      if(h_C[r] > 10e-7f)
	// each element in A is normalized C, then squre-rooted
	h_A[c * s.size() + r] = sqrtf( h_A[c * s.size() + r] / h_C[r]);
      else
	h_A[c * s.size() + r] = 1.0f/(double)tp.size(); // ad_hoc code to avoid 0 division
	
	//std::cout << " h_A " << h_A[c * s.size() + r] << std::endl;
	
	}
	
std::cout << " ok weight " << std::endl;

for (unsigned int i=0; i<s.size(); i++)
    {
        //serr<<"addForce() between "<<springs[i].m1<<" and "<<closestPos[springs[i].m1]<<sendl;
		    //if (t < 20*niterations) 
			for (unsigned int j=0; j<tp.size(); j++)
			{
			double fact = 0.005;//(double)0.0001/tp.size();	
			Real weight = (Real)(fact)*h_A[i * s.size() + j];
			//std::cout << weight << " pitch " << pitchA << std::endl;
			/*
			std::cout << weight << " pitch " << pitchA << std::endl;
			}
			{*/
			//closestPos[j] += ((Real)1/tp.size()) * tp[i];			
			this->addSpringForceSoft(m_potentialEnergy,f,x,v, i, j, s[i],weight);
			}
			//else this->addSpringForce3(m_potentialEnergy,f,x,v, i, s[i]);
		/*else{
		if(!targetBorder[(int)closestSource[i].begin()->second])
        this->addSpringForce(m_potentialEnergy,f,x,v, i, s[i]);
		else this->addSpringForce2(m_potentialEnergy,f,x,v, i, s[i]);}*/
    }
	
	//getchar();
    _f.endEdit();
	
	time = ((double)getTickCount() - time)/getTickFrequency();
    cout << "Time addforce " << time << endl;
	timeOverall += time;
	
	timei = (double)getTickCount();
	timeTotal =((double)getTickCount() - timeT)/getTickFrequency();
	cout << "Time total " << timeTotal << endl;
		
}

template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::generateData(const core::MechanicalParams* /*mparams*/,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
	
	int t = (int)this->getContext()->getTime();
	 double timeT;

	if (t == 0)
	{
		timeOverall = 0;
		timeTotal = 0;
		timeT = (double)getTickCount();
	}
	else if (t > 0 && t%niterations == 0)
	{
    std::cout << " ok init seg " << std::endl;
	double time0 = (double)getTickCount();
	double time1 = (double)getTickCount();
	}
	
    if(ks.getValue()==0) return;

    VecDeriv&        f = *_f.beginEdit();           //WDataRefVecDeriv f(_f);
    const VecCoord&  x = _x.getValue();			//RDataRefVecCoord x(_x);
    const VecDeriv&  v = _v.getValue();			//RDataRefVecDeriv v(_v);
    ReadAccessor< Data< VecCoord > > tn(targetNormals);
    ReadAccessor< Data< VecCoord > > tp(targetPositions);
	ReadAccessor< Data< VecCoord > > tcp(targetContourPositions);
	ReadAccessor< Data< VecCoord > > ssn(sourceSurfaceNormalsM);

    const vector<Spring>& s = this->springs.getValue();
    this->dfdx.resize(s.size());
    this->closestPos.resize(s.size());
		
	if (t < 2)
	setViewPointData();
	
	getSourceVisible();

	if (t > 0 && t%niterations == 0)
	{
	double timertt = (double)getTickCount();
	rtt = new cv::Mat;
	cv::Mat rtt_;
    renderToTexture(rtt_);
	*rtt = rtt_;
	listrtt.push_back(rtt);
	timertt = ((double)getTickCount() - timertt)/getTickFrequency();
    cout << "Time RTT " << timertt << endl;
	
	pcd = new std::vector<Vec3d>;
	pcd->resize(0);
	Vec3d pcd0;
	
	visible = new std::vector<bool>;
	visible->resize(0);
	bool vis0;
	
	for (unsigned int i=0; i<x.size(); i++)
	{		
		if(sourceVisible[i]) 
		{
			vis0 = true;
		} 
		else vis0 = false;
		
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
		}
			
	//std::cout << " time " << t << " " << listpcd.size() << std::endl;
	
	//if (t == nimages*niterations-4*niterations - 1)
	if (t == 1298)
	{
		writeData();
		getchar();
	}
	
	int key=cvWaitKey(1);
	if ((char)key == 27) return;
	
	//computeTargetNormals();
	
	double time = (double)getTickCount();

	m_potentialEnergy = 0;
	
	/*int ind0 = 145;
	int ind1 = 132;
	int ind2 = 139;
	int ind3 = 126;
	int ind4 = 0;*/
	
	int ind0 = 335;
	int ind1 = 15;
	int ind3 = 79;
	int ind2 = 399;
	//int ind3 = 399;
	int ind4 = 0;

	Vector3 trans0, trans1, trans2, trans3, trans4;
	
    if (t == 0)
		{
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
	
	trans0[0] = x[ind0][0]-0.05;
	trans0[1] = x[ind0][1] - 0.1;
	trans0[2] = x[ind0][2] - 0.0;
	
	trans1[0] = x[ind1][0] + 0.0;
	trans1[1] = x[ind1][1] - 0.1;
	trans1[2] = x[ind1][2] + 0.12;
	
	trans2[0] = x[ind2][0] - 0.0;
	trans2[1] = x[ind2][1] -0.1;
	trans2[2] = x[ind2][2] -0.12;
	trans3[0] = x[ind3][0] - 0.05;
	trans3[1] = x[ind3][1] -0.1;
	trans3[2] = x[ind3][2] -0.0;
	
	trans4[0] = x[ind4][0] ;
	trans4[1] = x[ind4][1] -0.15;
	trans4[2] = x[ind4][2] +0.0;
	
	closestPos[ind1]=trans1;
	closestPos[ind0]=trans0;
	closestPos[ind2]=trans2;
	closestPos[ind3]=trans3;
	closestPos[ind4]=trans4;
	}
	
			


	if (t > 50) this->addSpringForce(m_potentialEnergy,f,x,v, ind0, s[ind0]);
	if (t > 50) this->addSpringForce(m_potentialEnergy,f,x,v, ind1, s[ind1]);
	if (t > 50) this->addSpringForce(m_potentialEnergy,f,x,v, ind2, s[ind2]);
	if (t > 50) this->addSpringForce(m_potentialEnergy,f,x,v, ind3, s[ind3]);
	/*if (t > 200) this->addSpringForce(m_potentialEnergy,f,x,v, ind4, s[ind4]);*/

		/*else{
		if(!targetBorder[(int)closestSource[i].begin()->second])
        this->addSpringForce(m_potentialEnergy,f,x,v, i, s[i]);
		else this->addSpringForce2(m_potentialEnergy,f,x,v, i, s[i]);}*/		
	
	
	std::cout << " Error " << error << std::endl;
    _f.endEdit();
	
	time = ((double)getTickCount() - time)/getTickFrequency();
    cout << "Time addforce " << time << endl;
	timeOverall += time;
	
	timei = (double)getTickCount();
	timeTotal =((double)getTickCount() - timeT)/getTickFrequency();
	cout << "Time total " << timeTotal << endl;
		
}

template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::addForce(const core::MechanicalParams* /*mparams*/,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
generateData(_f, _x, _v);
}




template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
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
void GenerateSyntheticData<DataTypes, DepthTypes>::addSpringDForce(VecDeriv& df,const  VecDeriv& dx, int i, const Spring& spring, double kFactor, double /*bFactor*/)
{
    const int a = spring.m1;
    const Coord d = -dx[a];
    Deriv dforce = this->dfdx[i]*d;
    dforce *= kFactor;
    df[a]+=dforce;
    //serr<<"addSpringDForce, a="<<a<<", b="<<b<<", dforce ="<<dforce<<sendl;
}

template <class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::addDForce(const core::MechanicalParams* mparams,DataVecDeriv& _df , const DataVecDeriv&  _dx )
{

    VecDeriv& df = *_df.beginEdit();		//WDataRefVecDeriv df(_df);
    const VecDeriv&  dx = _dx.getValue();	// RDataRefVecDeriv dx(_dx);

    double kFactor       =  mparams->kFactor();
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
void GenerateSyntheticData<DataTypes, DepthTypes>::addKToMatrix(const core::MechanicalParams* mparams,const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    if(ks.getValue()==0) return;

    double kFact = mparams->kFactor();

    sofa::core::behavior::MultiMatrixAccessor::MatrixRef mat = matrix->getMatrix(this->mstate);
    if (!mat) return;
    const vector<Spring >& ss = this->springs.getValue();
    const unsigned int n = ss.size() < this->dfdx.size() ? ss.size() : this->dfdx.size();
    for (unsigned int e=0; e<n; e++)
    {
        const Spring& s = ss[e];
        unsigned p1 = mat.offset+Deriv::total_size*s.m1;
        const Mat& m = this->dfdx[e];
        for(int i=0; i<N; i++)
            for (int j=0; j<N; j++)
            {
                Real k = (Real)(m[i][j]*kFact);
                mat.matrix->add(p1+i,p1+j, -k);
            }
    }
}


// Function turn a cv::Mat into a texture, and return the texture ID as a GLuint for use
GLuint matToTexture(cv::Mat &mat, GLenum minFilter, GLenum magFilter, GLenum wrapFilter)
{
	// Generate a number for our textureID's unique handle
	GLuint textureID;
	glGenTextures(1, &textureID);
 
	// Bind to our texture handle
	glBindTexture(GL_TEXTURE_2D, textureID);
 
	// Catch silly-mistake texture interpolation method for magnification
	if (magFilter == GL_LINEAR_MIPMAP_LINEAR  ||
	    magFilter == GL_LINEAR_MIPMAP_NEAREST ||
	    magFilter == GL_NEAREST_MIPMAP_LINEAR ||
	    magFilter == GL_NEAREST_MIPMAP_NEAREST)
	{
		cout << "You can't use MIPMAPs for magnification - setting filter to GL_LINEAR" << endl;
		magFilter = GL_LINEAR;
	}
 
	// Set texture interpolation methods for minification and magnification
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);
 
	// Set texture clamping method
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapFilter);
 
	// Set incoming texture format to:
	// GL_BGR       for CV_CAP_OPENNI_BGR_IMAGE,
	// GL_LUMINANCE for CV_CAP_OPENNI_DISPARITY_MAP,
	// Work out other mappings as required ( there's a list in comments in main() )
	GLenum inputColourFormat = GL_BGR;
	if (mat.channels() == 1)
	{
		inputColourFormat = GL_LUMINANCE;
	}
 
	// Create the texture
	glTexImage2D(GL_TEXTURE_2D,     // Type of texture
	             0,                 // Pyramid level (for mip-mapping) - 0 is the top level
	             GL_RGB,            // Internal colour format to convert to
	             mat.cols,          // Image width  i.e. 640 for Kinect in standard mode
	             mat.rows,          // Image height i.e. 480 for Kinect in standard mode
	             0,                 // Border width in pixels (can either be 1 or 0)
	             inputColourFormat, // Input image format (i.e. GL_RGB, GL_RGBA, GL_BGR etc.)
	             GL_UNSIGNED_BYTE,  // Image data type
	             mat.ptr());        // The actual image data itself
 
	// If we're using mipmaps then generate them. Note: This requires OpenGL 3.0 or higher
	if (minFilter == GL_LINEAR_MIPMAP_LINEAR  ||
	    minFilter == GL_LINEAR_MIPMAP_NEAREST ||
	    minFilter == GL_NEAREST_MIPMAP_LINEAR ||
	    minFilter == GL_NEAREST_MIPMAP_NEAREST)
	{
		glGenerateMipmap(GL_TEXTURE_2D);
	}
 
	return textureID;
}

            
template<class DataTypes, class DepthTypes>
void GenerateSyntheticData<DataTypes, DepthTypes>::draw(const core::visual::VisualParams* vparams)
{
	
	int t = (int)this->getContext()->getTime();
		
	double timef = 0;
	timef = (double)getTickCount();
	//if (t > 0 && t%niterations == niterations-1)
	timeOverall += (timef - timei)/getTickFrequency();
	{
    cout <<" t " << t << " Time overall draw " << timei << " " << timef << " " << (timef - timei)/getTickFrequency() << endl;
	}
	
	if (t > 0 && t%niterations == niterations-1)
	{
	std::cout << " t " << t << " " << t/niterations <<  " time overall " << timeOverall << std::endl;
	timeFile << t/niterations;	
	timeFile << "\t";
	timeFile << timeOverall;
	timeFile << "\n";
	
	}
	
	sofa::gui::BaseGUI *gui = sofa::gui::GUIManager::getGUI();
    if (!gui)
    {
        std::cout << " no gui " << std::endl; 
    }
	
	/*sofa::gui::BaseViewer * viewer = gui->getViewer();
		std::string opath = "out/images16/img1%06d.png";
	
		if (t%niterations == 0)
		{
		char buf1[FILENAME_MAX];
        sprintf(buf1, opath.c_str(), iter_im);
        std::string filename1(buf1);
		viewer->setBackgroundImage( filename1);
		}*/
	
    if(ks.getValue()==0) return;

    if (this->closestPos.size()!=springs.getValue().size()) return;
    if (!vparams->displayFlags().getShowForceFields() && !drawColorMap.getValue()) return;

    ReadAccessor< Data< VecCoord > > x(*this->getMState()->read(core::ConstVecCoordId::position()));
    //const VecCoord& x = *this->mstate->getX();
    const vector<Spring>& springs = this->springs.getValue();

    /*if (vparams->displayFlags().getShowForceFields())
    {
        std::vector< Vector3 > points;
        for (unsigned int i=0; i<springs.size(); i++)
            if(!sourceIgnored[i])
            {
                Vector3 point1 = DataTypes::getCPos(x[springs[i].m1]);
                Vector3 point2 = DataTypes::getCPos(this->closestPos[i]);
				//std::cout << " pt " << point2[0] << " " << point2[1] << " " << point2[2] << std::endl;
                points.push_back(point1);
                points.push_back(point2);
            }

        const Vec<4,float> c(0,1,0.5,1);
        if (showArrowSize.getValue()==0 || drawMode.getValue() == 0)	vparams->drawTool()->drawLines(points, 1, c);
        else if (drawMode.getValue() == 1)	for (unsigned int i=0;i<points.size()/2;++i) vparams->drawTool()->drawCylinder(points[2*i+1], points[2*i], showArrowSize.getValue(), c);
        else if (drawMode.getValue() == 2)	for (unsigned int i=0;i<points.size()/2;++i) vparams->drawTool()->drawArrow(points[2*i+1], points[2*i], showArrowSize.getValue(), c);
        else serr << "No proper drawing mode found!" << sendl;
    }*/
	
        glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT | GL_ENABLE_BIT);
        glDisable( GL_LIGHTING);
    //if(drawColorMap.getValue())
	if(!disp)
    {
        /*std::vector< Real > dists(x.size());  for (unsigned int i=0; i<dists.size(); i++) dists[i]=0.;
        for (unsigned int i=0; i<springs.size(); i++)
            if(!sourceIgnored[i])
            {
                Vector3 point1 = DataTypes::getCPos(x[springs[i].m1]);
                Vector3 point2 = DataTypes::getCPos(this->closestPos[i]);
                dists[springs[i].m1]=(point2-point1).norm();
            }
        Real max=0; for (unsigned int i=0; i<dists.size(); i++) if(max<dists[i]) max=dists[i];

        glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT | GL_ENABLE_BIT);
        glDisable( GL_LIGHTING);

        ReadAccessor< Data< helper::vector< tri > > > t(sourceTriangles);

    //std::cout << " start draw 2" << std::endl;

        //if(t.size()) // mesh visu
        {
            glBegin( GL_TRIANGLES);
            for ( unsigned int i = 0; i < t.size(); i++)
            {
                for ( unsigned int j = 0; j < 3; j++)
                {
                    const unsigned int& indexP = t[i][j];
                    //sofa::helper::gl::Color::setHSVA(dists[indexP]*240./max,1.,.8,1.);
					sofa::helper::gl::Color::setHSVA(240,1.,.8,1.);
                    //glVertex3d(x[indexP][0],x[indexP][1],x[indexP][2]);
                }
            }
            glEnd();
        }*/

        //else // point visu
        /*{
            glPointSize( 10);
            glBegin( GL_POINTS);
            for (unsigned int i=0; i<x.size(); i++)
            {
				if(sourceVisible[i]){
                //sofa::helper::gl::Color::setHSVA(dists[i]*240./max,1.,.8,1.);
                glVertex3d(x[i][0],x[i][1],x[i][2]);}
            }
            glEnd();
            glPointSize( 1);

        }*/
      	ReadAccessor< Data< VecCoord > > xtarget(targetPositions);
		ReadAccessor< Data< VecCoord > > xtargetcontour(targetContourPositions);
		ReadAccessor< Data< VecCoord > > xsourcecontour(sourceContourPositions);
		ReadAccessor< Data< VecCoord > > ssn(sourceSurfaceNormalsM);
        glPointSize( 3);
        glBegin( GL_POINTS);
//std::cout << " size target " << targetPositions.getValue().size() << " dist size " << dists.size() << std::endl;
       
 for (unsigned int i=0; i< targetPositions.getValue().size(); i++)
          {
            //sofa::helper::gl::Color::setHSVA(dists[i]*240./max,1.,.8,1.);
			//if(targetBackground[i])
			//if(!targetBorder[i])
				{
			//sofa::helper::gl::Color::setHSVA(40,1.,.8,1.);
			//else 
			sofa::helper::gl::Color::setHSVA(50,1.,.8,1.);
			//glVertex3d(xtarget[i][0],xtarget[i][1],xtarget[i][2]);
			}
          }
		  
		  std::cout << "size " << sourceContourPositions.getValue().size() << std::endl;
		for (unsigned int i=0; i< sourceContourPositions.getValue().size(); i++)
          {
			sofa::helper::gl::Color::setHSVA(140,1.,.8,1.);
            //glVertex3d(xsourcecontour[i][0],xsourcecontour[i][1],xsourcecontour[i][2]);
          } 
		  
		//if (disp)  
		glBegin( GL_POINTS);

		/*if (t > 20)
		for (unsigned int i=0; i< x.size(); i++)
          {
			  //if (sourceBorder[i])
			 //if ( targetWeights[(int)closestSource[i].begin()->second] == 1)
			 if ( sourceVisible[i])
			  {
			sofa::helper::gl::Color::setHSVA(40,2.,.8,2.);
            glVertex3d(x[i][0],x[i][1],x[i][2]);
			  }
          }*/
		  
		for (unsigned int i=0; i< targetContourPositions.getValue().size(); i++)
		{
			sofa::helper::gl::Color::setHSVA(50,1.,.8,1.);
            //glVertex3d(xtargetcontour[i][0],xtargetcontour[i][1],xtargetcontour[i][2]);
		}

		if (!disp)
		for (unsigned int i=0; i< x.size(); i++)
          {
		   // if (sourceBorder[i])
			  {
			/*unsigned int id = closestSource[i].begin()->second;

			        std::vector< Vector3 > points;
                Vector3 point1 = DataTypes::getCPos(x[i]);
                Vector3 point2 = DataTypes::getCPos(xtargetcontour[id]);
				//std::cout << " pt " << point2[0] << " " << point2[1] << " " << point2[2] << std::endl;
                points.push_back(point1);
                points.push_back(point2);
							const Vec<4,float> c(0,1,0.5,1);
            //vparams->drawTool()->drawLines(points, 1, c);
			
			            //sofa::helper::gl::Color::setHSVA(dists[i]*240./max,1.,.8,1.);
			sofa::helper::gl::Color::setHSVA(100,2.,.8,2.);*/
            //glVertex3d(xtargetcontour[id][0],xtargetcontour[id][1],xtargetcontour[id][2]);
			        if (useVisible.getValue())
					if (t>20 && sourceVisible[i]){
			sofa::helper::gl::Color::setHSVA(140,2.,.8,2.);
            //glVertex3d(x[i][0],x[i][1],x[i][2]);
			}

			  }
          }
		  
		   const Vec<4,float> c(0,1,0.5,1);
		   int kk = 0;
		   std::cout << "size indices " << indices.size() << std::endl;
		  
		  	        std::vector< Vector3 > points;
					/*if (t>51){
        for (unsigned int i=0; i<x.size(); i++)
            if(sourceBorder[i])
            {
				points.resize(0);
                //Vector3 point1 = DataTypes::getCPos(x[springs[i].m1]);
				Vector3 point1 = x[i];
				int id=indices[kk];
				//std::cout << " i draw " << i << " ind " << indices[kk] << std::endl;

				kk++; 
				//closestPos[i]=tcp[id]*projF;
                //Vector3 point2 = DataTypes::getCPos(this->closestPos[i]);
				Vector3 point2 = xtargetcontour[id];
				//std::cout << " pt " << point2[0] << " " << point2[1] << " " << point2[2] << std::endl;
                points.push_back(point1);
                //points.push_back(point1+ssn[i]);
				points.push_back(point2);
			sofa::helper::gl::Color::setHSVA(140,5.,.8,4.);
            //glVertex3d(point1[0],point1[1],point1[2]);
			//if (t>30)
			//vparams->drawTool()->drawLines(points, 3, c);
            }
			
					}*/
		
			//for (unsigned int i=0;i<points.size()/2;++i) vparams->drawTool()->drawArrow(points[2*i+1], points[2*i], 0.01, c);
        glEnd();
        glPointSize( 1);

        glPopAttrib();
    }

}


}
}
} // namespace sofa


