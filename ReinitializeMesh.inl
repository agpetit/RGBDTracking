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

#define SOFA_RGBDTRACKING_ReinitializeMesh_INL

#include <sofa/component/loader/MeshObjLoader.h>
#include <sofa/component/engine/NormalsFromPoints.h>
#include <limits>
#include <set>
#include <iterator>
#include <sofa/helper/gl/Color.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/component/visualmodel/InteractiveCamera.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/simulation/common/Simulation.h>


#ifdef Success
  #undef Success
#endif

#include <pcl/common/common_headers.h>
#include <pcl/point_cloud.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/point_types.h>
#include <pcl/impl/point_types.hpp>

#include <algorithm>    // std::max
#include "ReinitializeMesh.h"


using std::cerr;
using std::endl;

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;
using namespace helper;


template <class DataTypes>
ReinitializeMesh<DataTypes>::ReinitializeMesh(core::behavior::MechanicalState<DataTypes> *mm )
 : Inherit(mm){
	iter_im = 0;
	color = cv::Mat::zeros(240,320, CV_8UC3);
	depth = cv::Mat::zeros(240,320,CV_32FC1);
	depth_1 = depth;
	Vector4 camParam = cameraIntrinsicParameters.getValue();
	
    rgbIntrinsicMatrix(0,0) = camParam[0];
	rgbIntrinsicMatrix(1,1) = camParam[1];
	rgbIntrinsicMatrix(0,2) = camParam[2];
	rgbIntrinsicMatrix(1,2) = camParam[3];
	
	rgbIntrinsicMatrix(0,0) = 275.34;
	rgbIntrinsicMatrix(1,1) = 275.34;
	//rgbIntrinsicMatrix(0,2) = 157.25;
	//rgbIntrinsicMatrix(1,2) = 117.75;
	rgbIntrinsicMatrix(0,2) = 162;
	rgbIntrinsicMatrix(1,2) = 120;
	
	rectRtt.x = 0;
	rectRtt.y = 0;
	rectRtt.height = 480;
	rectRtt.width = 640;
	//
    // initialize paramters
}

template <class DataTypes>
ReinitializeMesh<DataTypes>::~ReinitializeMesh()
{
}

template<class DataTypes>
void ReinitializeMesh<DataTypes>::getSourceVisible(double znear, double zfar)
{
	
	int t = (int)this->getContext()->getTime();

	//if (t%2 == 0)
    {

	int hght = 480;
	int wdth = 640;
			
	//sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
	//sofa::component::visualmodel::BaseCamera::SPtr currentCamera;// = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();
	//root->get(currentCamera);
	
	//double znear = currentCamera->getZNear();
	//double zfar = currentCamera->getZFar();
	//std::cout << " no viewer 2" << std::endl; 

	
	cv::Mat _rtd,_rtd0;
	_rtd.create(hght, wdth, CV_32F);
	depthrend.create(hght, wdth, CV_32F);
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
	//glDepthMask(GL_TRUE);
	//glDepthFunc(GL_ALWAYS); // Change this to whatever kind of depth testing you want
    //glDepthRange(0.0f, 1.0f);

	double time3 = (double)getTickCount();

	/*viewport[0] = rectRtt.x;
	viewport[1] = rectRtt.y;
	viewport[2] = rectRtt.width;
	viewport[3] = rectRtt.height;*/
	
	GLfloat depths1[rectRtt.width * rectRtt.height];
		
	//std::cout << "viewport " << viewport[0] << " " << viewport[1] << " " << viewport[2] << std::endl;
	
    glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3], GL_DEPTH_COMPONENT, GL_FLOAT, &depths);
	
	//time3 = ((double)getTickCount() - time3)/getTickFrequency();
	//cout << "Rigid glreadpixel " << time3 << endl;
    //glReadBuffer(GL_BACK);
	
	//double zNear = 1e10;
	//double zFar = -1e10;
	
	//std::cout << " znear " << znear << std::endl;	
	//std::cout << " zfar " << zfar << std::endl;
	
	std::vector<cv::Point> ptfgd;
	ptfgd.resize(0);
	cv::Point pt;

	for (int j = 0; j < wdth; j++)
		for (int i = 0; i< hght; i++)
		{
			if ((double)(float)depths[j+i*wdth]	< 1){
			//if (j >= rectRtt.x && j < rectRtt.x + rectRtt.width && i >= rectRtt.y && i < rectRtt.y + rectRtt.height) {
			//if ((double)(float)depths1[j-rectRtt.x+(i-rectRtt.y)*(rectRtt.width)]	< 1){
			//if ((double)(float)depths1[j-rectRtt.x+(i-rectRtt.y)*(rectRtt.width)]	< 1){
			_rtd0.at<uchar>(hght-i-1,j) = 255;//(int)100000*(1-depths[j+i*wdth]);
			
			double clip_z = (depths[j+i*wdth] - 0.5) * 2.0;
			//double clip_z = (depths1[j-rectRtt.x+(i-rectRtt.y)*(rectRtt.width)] - 0.5) * 2.0;
            depthsN[j+i*wdth] = 2*znear*zfar/(clip_z*(zfar-znear)-(zfar+znear));
			//std::cout << " depth " << world_z << std::endl;
			pt.x = j;
			pt.y = i;
			ptfgd.push_back(pt);
			}
			else 
			{
				_rtd0.at<uchar>(hght-i-1,j) = 0;
				depthsN[j+i*wdth] = 0;
			}
			
			depthrend.at<float>(hght-i-1,j) = depths[j+i*wdth];

				//std::cout << " depth " << (int)_rtd0.at<uchar>(i,j) << std::endl;
		//if (depths[j+i*319]	> 0)
			//_rtd0.at<uchar>(j,i) = 255;
		//}
		/*else 
				{
				_rtd0.at<uchar>(hght-i-1,j) = 0;
				depthsN[j+i*wdth] = 0;
				}*/
		
		}
		
		rtd = _rtd0;
		
		//cv::imwrite("rtd0.png", _rtd0);
				
		int t = (int)this->getContext()->getTime();
		
		/*if (t > 0)
		cv::imwrite("depth.png", _rtd0);*/

		if(t > 0)
		{
		rectRtt = cv::boundingRect(ptfgd);
	    //std::cout << " ptfgd " << frame_count << " rect1 " << rectangle.x << " " << rectangle.y << " rect2 " << rectangle.width << " " << rectangle.height << std::endl;

	    rectRtt.x -= 10;
	    rectRtt.y -= 10;
	    rectRtt.height += 20;
	    rectRtt.width += 20;
		
		}
		
		depthMap = _rtd0;
		
		//cv::imwrite("depth.png", _rtd0);
	
const VecCoord& x = *this->mstate->getX();
 	
Eigen::Matrix3f rgbIntrinsicMatrixV;
rgbIntrinsicMatrixV(0,0) = 2*rgbIntrinsicMatrix(0,0);
rgbIntrinsicMatrixV(1,1) = 2*rgbIntrinsicMatrix(1,1);
/*rgbIntrinsicMatrix(0,2) = 2*157.25;
rgbIntrinsicMatrix(1,2) = 2*117.75;*/
rgbIntrinsicMatrixV(0,2) = 2*rgbIntrinsicMatrix(0,2);
rgbIntrinsicMatrixV(1,2) = 2*rgbIntrinsicMatrix(1,2);

sourceVisible.resize(x.size());
int nvisible = 0;

for (int k = 0; k < x.size(); k++)
{
   //std::cout << " xk0 " << x[k][0] << " " << x[k][1] << " xk2 " << x[k][2] << std::endl; 
	int x_u = (int)(x[k][0]*rgbIntrinsicMatrixV(0,0)/x[k][2] + rgbIntrinsicMatrixV(0,2));
    int x_v = (int)(x[k][1]*rgbIntrinsicMatrixV(1,1)/x[k][2] + rgbIntrinsicMatrixV(1,2));
	
	//std::cout << " depths " << x_u << " " << x_v << " " << x[k][2] << " " << depthsN[x_u+(hght-x_v-1)*wdth] << std::endl;
	//std::cout << " depths0 " << (float)depthsN[x_u+(hght-x_v-1)*wdth] << " " << hght * wdth << std::endl;
	
if ((float)abs(depthsN[x_u+(hght-x_v-1)*wdth]+(float)x[k][2]) < visibilityThreshold.getValue() || (float)depthsN[x_u+(hght-x_v-1)*wdth] == 0)
	{
	sourceVisible[k] = true;
	nvisible ++;
	}
	else {
	sourceVisible[k] = false;	
	}
	
}

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
	
	std::cout << " npoints " << x.size() << " " << nvisible << std::endl;

sourceVisiblePositions.setValue(sourceVis);
		
	}
}

template<class DataTypes>
void ReinitializeMesh<DataTypes>::updateSourceVisible()
{
// build k-d tree
const VecCoord&  x = *this->mstate->getX();
			
VecCoord sourceVis;
sourceVis.resize(sourceVisiblePositions.getValue().size());

Vector3 pos;
int k = 0;
			
for (unsigned int i=0; i< x.size(); i++)
	{
		if (sourceVisible[i])
		{
            pos = x[i];
            sourceVis[k]=pos;
			k++;			
		}
	}
sourceVisiblePositions.setValue(sourceVis);
	
}

template<class DataTypes>
void ReinitializeMesh<DataTypes>::extractSourceContour()
{
	
	/*if (!useVisible.getValue())
	getSourceVisible();*/
	
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
	

	Eigen::Matrix3f rgbIntrinsicMatrixV;
	rgbIntrinsicMatrixV(0,0) = 2*rgbIntrinsicMatrix(0,0);
	rgbIntrinsicMatrixV(1,1) = 2*rgbIntrinsicMatrix(1,1);
	/*rgbIntrinsicMatrix(0,2) = 2*157.25;
	rgbIntrinsicMatrix(1,2) = 2*117.75;*/
	rgbIntrinsicMatrixV(0,2) = 2*rgbIntrinsicMatrix(0,2);
	rgbIntrinsicMatrixV(1,2) = 2*rgbIntrinsicMatrix(1,2);
	
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
		int x_u = (int)(x[i][0]*rgbIntrinsicMatrixV(0,0)/x[i][2] + rgbIntrinsicMatrixV(0,2));
		int x_v = (int)(x[i][1]*rgbIntrinsicMatrixV(1,1)/x[i][2] + rgbIntrinsicMatrixV(1,2));
		int thickness = 1;
        int lineType = 2;

				if (dist0.at<uchar>(x_v,x_u) < borderThdSource.getValue()/*7*/)
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
				        circle( contourpoints,
         cv::Point(x_u,x_v),
         wdth/128.0,
         Scalar( 0, 0, 255 ),
         thickness,
         lineType );
				}
				else sourceBorder[i] = false;
				
				//sourceWeights.push_back((double)1./(0.12*(1.0+sqrt(dist0.at<uchar>(x_v,x_u)))));
				
				sourceWeights.push_back((double)exp(-dist0.at<uchar>(x_v,x_u)/sigmaWeight.getValue()));

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
	
	
	std::cout << " n source " << nbs << " n source contour " << nsourcecontour << std::endl;
	
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

template<class DataTypes>
void ReinitializeMesh<DataTypes>::extractSourceVisibleContour()
{
	
	int hght = 480;
	int wdth = 640;
		
	double cannyTh1 = 350;
	double cannyTh2 = 10;
	cv::Mat contour,dist,dist0;
	
	//cv::imwrite("depthmap.png", depthMap);
	
	cv::Canny( depthMap, contour, cannyTh1, cannyTh2, 3);
    contour = cv::Scalar::all(255) - contour;

	cv::distanceTransform(contour, dist, CV_DIST_L2, 3);
	
   // cv::imwrite("contour.png", contour);

//dt0 *= 5000;
//pow(dt0, 0.5, dt0);
    dist.convertTo(dist0, CV_8U, 1, 0);
	
	int ncontour = 0;
	pcl::PointCloud<pcl::PointXYZRGB> sourceContour;
	
	Eigen::Matrix3f rgbIntrinsicMatrixV;
	rgbIntrinsicMatrixV(0,0) = 2*rgbIntrinsicMatrix(0,0);
	rgbIntrinsicMatrixV(1,1) = 2*rgbIntrinsicMatrix(1,1);
	/*rgbIntrinsicMatrix(0,2) = 2*157.25;
	rgbIntrinsicMatrix(1,2) = 2*117.75;*/
	rgbIntrinsicMatrixV(0,2) = 2*rgbIntrinsicMatrix(0,2);
	rgbIntrinsicMatrixV(1,2) = 2*rgbIntrinsicMatrix(1,2);
	
	const VecCoord& x =  *this->mstate->getX();
	
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
		int x_u = (int)(x[i][0]*rgbIntrinsicMatrixV(0,0)/x[i][2] + rgbIntrinsicMatrixV(0,2));
		int x_v = (int)(x[i][1]*rgbIntrinsicMatrixV(1,1)/x[i][2] + rgbIntrinsicMatrixV(1,2));
		int thickness = 1;
        int lineType = 2;

        circle( contourpoints,
         cv::Point(x_u,x_v),
         wdth/128.0,
         Scalar( 0, 0, 255 ),
         thickness,
         lineType );
				if (dist0.at<uchar>(x_v,x_u) < borderThdSource.getValue() /*6*/ & sourceVisible[i])
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
				
				sourceWeights.push_back((double)exp(-dist0.at<uchar>(x_v,x_u)/sigmaWeight.getValue()));

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
	std::cout << " n source " << nbs << " n source contour " << nsourcecontour << " " << borderThdSource.getValue() << std::endl;
	
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

template<class DataTypes>
void ReinitializeMesh<DataTypes>::updateSourceVisibleContour()
{
    // build k-d tree
    const VecCoord&  x = *this->mstate->getX();
			std::cout <<" size contour " << sourceContourPositions.getValue().size()<< std::endl;

	
		VecCoord sourcecontourpos;
    sourcecontourpos.resize(sourceContourPositions.getValue().size());
            Vector3 pos;
			int k = 0;


	for (unsigned int i=0; i<x.size(); i++)
	{
		if (sourceBorder[i])
		{
            pos = x[i];
            sourcecontourpos[k]=pos;	
			k++;		
	    }
	}
    const VecCoord&  p = sourcecontourpos;
	sourceContourPositions.setValue(p);
	
}
	
/*template<class DataTypes, class DepthTypes>
void ClosestPointRegistrationForceFieldCam<DataTypes, DepthTypes>::extractSourceContour()
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




}
}
} // namespace sofa

