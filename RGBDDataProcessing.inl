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

#define SOFA_RGBDTRACKING_RGBDDATAPROCESSING_INL

#include <limits>
#include <iterator>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/io.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/principal_curvatures.h>
#include <pcl/gpu/features/features.hpp>
#include "DataSource.hpp"

#include <sofa/helper/gl/Color.h>
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseVisual/BaseCamera.h>
#include <SofaBaseVisual/InteractiveCamera.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/simulation/Simulation.h>
#include <pcl/keypoints/sift_keypoint.h>

#include <pcl/features/fpfh_omp.h>
#include <pcl/features/pfh.h>
#include <pcl/features/pfhrgb.h>
#include <pcl/features/3dsc.h>
#include <pcl/features/shot_omp.h>
#include <pcl/kdtree/kdtree_flann.h>

#include "ImageConverter.h"
#ifdef Success
  #undef Success
#endif

#include <algorithm>

#include "RGBDDataProcessing.h"

using std::cerr;
using std::endl;

namespace sofa
{

namespace core
{

namespace objectmodel
{

using namespace sofa::defaulttype;
using namespace helper;

template <class DataTypes>
RGBDDataProcessing<DataTypes>::RGBDDataProcessing( )
 : Inherit()
 	, cameraIntrinsicParameters(initData(&cameraIntrinsicParameters,Vector4(),"cameraIntrinsicParameters","camera parameters"))
	, useContour(initData(&useContour,false,"useContour","Emphasize forces close to the target contours"))
	, useRealData(initData(&useRealData,true,"useRealData","Use real data"))
	, useSensor(initData(&useSensor,true,"useSensor","Use real data"))
	, sensorType(initData(&sensorType, 0,"sensorType","Type of the sensor"))
        , niterations(initData(&niterations,1,"niterations","Number of iterations in the tracking process"))
	, nimages(initData(&nimages,1500,"nimages","Number of images to read"))
	, samplePCD(initData(&samplePCD,4,"samplePCD","Sample step for the point cloud"))
	, offsetX(initData(&offsetX,3,"offsetX","offset along x for the point cloud"))
	, offsetY(initData(&offsetY,0,"offsetY","offset along y for the point cloud"))
        , sigmaWeight(initData(&sigmaWeight,(Real)4,"sigmaWeight","sigma weights"))
	, borderThdPCD(initData(&borderThdPCD,4,"borderThdPCD","border threshold on the target silhouette"))
	, inputPath(initData(&inputPath,"inputPath","Path for data readings",false))
	, outputPath(initData(&outputPath,"outputPath","Path for data writings",false))
	, useDistContourNormal(initData(&useDistContourNormal,false,"outputPath","Path for data writings"))
	, windowKLT(initData(&windowKLT,1500,"nimages","Number of images to read"))
        , segNghb(initData(&segNghb,8,"segnghb","Neighbourhood for segmentation"))
        , segImpl(initData(&segImpl,1,"segimpl","Implementation mode for segmentation (CUDA, OpenCV)"))
        , segMsk(initData(&segMsk,1,"segmsk","Mask type for segmentation"))
        , scaleImages(initData(&scaleImages,1,"downscaleimages","Down scaling factor on the RGB and depth images"))
        , displayImages(initData(&displayImages,true,"displayimages","Option to display RGB and Depth images"))
        , displayDownScale(initData(&displayDownScale,1,"downscaledisplay","Down scaling factor for the RGB and Depth images to be displayed"))
        , saveImages(initData(&saveImages,false,"saveimages","Option to save RGB and Depth images on disk"))
        , displaySegmentation(initData(&displaySegmentation,true,"displaySegmentation","Option to display the segmented image"))
        , drawPointCloud(initData(&drawPointCloud,false,"drawPointCloud"," "))
        , displayBackgroundImage(initData(&displayBackgroundImage,false,"displayBackgroundImage"," "))
        , useCurvature(initData(&useCurvature,false,"useCurvature"," "))
        , scaleSegmentation(initData(&scaleSegmentation,1,"downscalesegmentation","Down scaling factor on the RGB image for segmentation"))
        , imagewidth(initData(&imagewidth,640,"imagewidth","Width of the RGB-D images"))
        , imageheight(initData(&imageheight,480,"imageheight","Height of the RGB-D images"))
        , targetPositions(initData(&targetPositions,"targetPositions","Points of the target point cloud."))
        , targetContourPositions(initData(&targetContourPositions,"targetContourPositions","Contour points of the target point cloud."))
        , targetWeights(initData(&targetWeights,"targetWeights","Weights for the points of the target point cloud."))
        , targetBorder(initData(&targetBorder,"targetBorder","Boolean for the points of the target point cloud that belong to the contour."))
        , cameraPosition(initData(&cameraPosition,"cameraPosition","Position of the camera w.r.t the point cloud"))
        , cameraOrientation(initData(&cameraOrientation,"cameraOrientation","Orientation of the camera w.r.t the point cloud"))
        , cameraChanged(initData(&cameraChanged,false,"cameraChanged","If the camera has changed or not"))
        , curvatures(initData(&curvatures,"curvatures","curvatures."))
        , useSIFT3D(initData(&useSIFT3D,false,"useSIFT3D"," "))
        , stopatinit(initData(&stopatinit,false,"stopatinit","stopatinit."))
{
	this->f_listening.setValue(true); 
	iter_im = 0;
        timeSegmentation = 0;
        timePCD = 0;

}

template <class DataTypes>
RGBDDataProcessing<DataTypes>::~RGBDDataProcessing()
{
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::init()
{		
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);		
        seg.init(segNghb.getValue(), segImpl.getValue(), segMsk.getValue());

	if(displayImages.getValue())
	{
	cv::namedWindow("image_sensor");
        cv::namedWindow("depth_sensor");
	}
        cv::namedWindow("image_segmented");

        Vector4 camParam = cameraIntrinsicParameters.getValue();

        rgbIntrinsicMatrix(0,0) = camParam[0];
        rgbIntrinsicMatrix(1,1) = camParam[1];
        rgbIntrinsicMatrix(0,2) = camParam[2];
        rgbIntrinsicMatrix(1,2) = camParam[3];

        initsegmentation = true;        

}


bool g_bDisplay = true;

bool destroy=false;
cv::Rect box;
bool drawing_box = false;

void draw_box(cv::Mat _img, cv::Rect rect)
{
  cv::rectangle(_img, cvPoint(box.x, box.y), cvPoint(box.x+box.width,box.y+box.height),
              cvScalar(0,0,255) ,2);

  cv::Rect rect2=cv::Rect(box.x,box.y,box.width,box.height);
  //cvSetImageROI(image, rect2);   //here I wanted to set the drawn rect as ROI
}

// Implement mouse callback
void my_mouse_callback( int event, int x, int y, int flags, void* param )
{
  cv::Mat* frame = (cv::Mat*) param;
  
  switch( event )
  {
      case CV_EVENT_MOUSEMOVE:
      {
          if( drawing_box )
          {
              box.width = x-box.x;
              box.height = y-box.y;
          }
      }
      break;

      case CV_EVENT_LBUTTONDOWN:
      {
          drawing_box = true;
          box = cvRect( x, y, 0, 0 );
      }
      break;

      case CV_EVENT_LBUTTONUP:
      {
          drawing_box = false;
          if( box.width < 0 )
          {
              box.x += box.width;
              box.width *= -1;
          }

          if( box.height < 0 )
          {
              box.y += box.height;
              box.height *= -1;
          }

          draw_box(*frame, box);
      }
      break;

      case CV_EVENT_RBUTTONUP:
      {
          destroy=true;
      }
      break;

      default:
      break;
   }

}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::initSegmentation()
{
	
    cv::Mat mask,maskimg,mask0,roimask,mask1; // segmentation result (4 possible values)
    cv::Mat bgModel,fgModel; // the models (internally used)
	
    cv::Mat downsampledbox,downsampled;

    int scaleSeg = scaleSegmentation.getValue();
    if (scaleSeg>1)
    cv::resize(color, downsampledbox, cv::Size(color.cols/scaleSeg, color.rows/scaleSeg));
    else downsampledbox = color.clone();
	
    //cv::imwrite("colorinit.png", color);
	
    cv::Mat temp;
    cv::Mat tempgs;
    cv::Mat imgs,imgklt;
    temp = downsampledbox.clone();
    cv::Mat temp1 = temp.clone();
		
    const char* name = "image";

    cv::namedWindow(name);
    box = cvRect(0,0,1,1);

    //cv::imshow("image",	color);
    //tempm.resize(image.step1());
	  
    // Set up the callback
    cv::setMouseCallback(name, my_mouse_callback, (void*) &temp);
	  
    std::cout << " Time " << (int)this->getContext()->getTime() << std::endl;

    // Main loop
    while(1)
    {
      if (destroy)
      {
        cv::destroyWindow(name); break;
      }
	  temp1 = temp.clone();

      if (drawing_box)
          draw_box(temp1, box);
		  
      cv::moveWindow(name, 200, 100);
      cv::imshow(name, temp1);
      //tempm.resize(image.step1());
      int key=waitKey(10);
      if ((char)key == 27) break;
	
    }
   // delete temp1;
    temp1.release();
    cv::setMouseCallback(name, NULL, NULL);
	
	    //cvDestroyWindow(name);
    cv::Rect rectangle(69,47,198,171);

    rectangle = box;
    seg.setRectangle(rectangle);
    
    if (scaleSeg>1)
    cv::resize(color, downsampled, cv::Size(color.cols/scaleSeg, color.rows/scaleSeg));
    else downsampled = color.clone();

    foregroundS = cv::Mat(downsampled.size(),CV_8UC3,cv::Scalar(255,255,255));

    seg.segmentationFromRect(downsampled,foregroundS);

    cv::resize(foregroundS, foreground, color.size());

    // draw rectangle on original image
    //cv::rectangle(image, rectangle, cv::Scalar(255,255,255),1);
    
    // display result
    if (displaySegmentation.getValue()){
    cv::imshow("image_segmented",foregroundS);
    cv::waitKey(1);
    }

    if (waitKey(20) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
   {
        cout << "esc key is pressed by user" << endl;
   }
   
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::segment()
{
    cv::Mat downsampled,downsampled1;
    double timef = (double)getTickCount();
    int scaleSeg = scaleSegmentation.getValue();
    if (scaleSeg>1)
    cv::resize(color, downsampled, cv::Size(color.cols/scaleSeg, color.rows/scaleSeg));
    else downsampled = color.clone();

    seg.updateMask(foregroundS);
    //cv::GaussianBlur( downsampled, downsampled1, cv::Size( 3, 3), 0, 0 );
    //cv::imwrite("downsampled.png", downsampled);
    seg.updateSegmentation(downsampled,foregroundS);
    cv::resize(foregroundS, foreground, color.size(), INTER_NEAREST);
    if(useContour.getValue())
    {
    cv::resize(seg.dotImage, dotimage, color.size(), INTER_NEAREST);
    cv::resize(seg.distImage, distimage, color.size(), INTER_NEAREST);
    }
    //foreground = foregroundS.clone();
    timeSegmentation = ((double)getTickCount() - timef)/getTickFrequency();
    std::cout << "TIME SEGMENTATION " << timeSegmentation << std::endl;

    //seg.updateSegmentationCrop(downsampled,foreground);

    // display result
    if (displaySegmentation.getValue()){
    cv::imshow("image_segmented",foregroundS);
    cv::waitKey(1);
    }
}

template<class DataTypes>
void RGBDDataProcessing<DataTypes>::segmentSynth()
{

	cv::Mat downsampled;
	//cv::resize(color, downsampled, cv::Size(color.cols/2, color.rows/2));
	downsampled = color;
	
	foreground = cv::Mat(downsampled.size(),CV_8UC4,cv::Scalar(255,255,255,0));
	
	//cv::imwrite("downsampled0.png",downsampled);
	
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
				
                        if (color.at<cv::Vec3b>(i,j)[0] == 255 && color.at<cv::Vec3b>(i,j)[1] == 255 && color.at<cv::Vec3b>(i,j)[2] == 255)
			{
				foreground.at<cv::Vec4b>(i,j)[3] = 0;
			}
			else foreground.at<cv::Vec4b>(i,j)[3] = 255;
		}

	cv::Mat distanceMap(foreground.size(),CV_8U,cv::Scalar(255));
	cv::Mat dot(foreground.size(),CV_8U,cv::Scalar(0));
	
		//std::cout << " ok seg synth 1 " << std::endl; 

	seg.filter(foreground,distanceMap,dot);	
	
        //cv::imwrite("downsampled1.png",foreground);

	//distImage = distanceMap;
	//dotImage = dot;
	
	/*cv::namedWindow("Image");
    cv::imshow("Image",downsampled);*/

    // display result
    cv::imshow("image_segmented",foreground);

}

template <class DataTypes>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr RGBDDataProcessing<DataTypes>::PCDFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	//pcl::PointCloud<pcl::PointXYZRGB> pointcloud; 
	outputPointcloud->points.resize(0);  
		
	int sample;
	int offsetx;
	
    switch (sensorType.getValue())
    {
    // FLOAT ONE CHANNEL
    case 0:
	sample = 2;
	offsetx = 0;
	break;
	case 1:
	
	sample = samplePCD.getValue();//3
	offsetx = offsetX.getValue();//3;
	break;
        }

	float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
        float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy
	pcl::PointXYZRGB newPoint;
	int offsety = offsetY.getValue();

	for (int i=0;i<(int)(depthImage.rows-offsety)/sample;i++)
	{
		for (int j=0;j<(int)(depthImage.cols-offsetx)/sample;j++)
		{

                        float depthValue = (float)depthImage.at<float>(sample*(i+offsety),sample*(j+offsetx))*0.819;
			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)rgbImage.at<Vec4b>(sample*i,sample*j)[3];
			if (avalue > 0 && depthValue>0)                // if depthValue is not NaN
			{
				// Find 3D position respect to rgb frame:
                                newPoint.z = depthValue;
				newPoint.x = (sample*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
				newPoint.y = (sample*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
				newPoint.r = rgbImage.at<cv::Vec4b>(sample*i,sample*j)[2];
				newPoint.g = rgbImage.at<cv::Vec4b>(sample*i,sample*j)[1];
				newPoint.b = rgbImage.at<cv::Vec4b>(sample*i,sample*j)[0];
				outputPointcloud->points.push_back(newPoint);
				
			}

		}
	}

	if (useGroundTruth.getValue())
	{
	int sample1 = 2;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud1(new pcl::PointCloud<pcl::PointXYZRGB>);

	for (int i=0;i<(int)(depthImage.rows-offsety)/sample1;i++)
	{
		for (int j=0;j<(int)(depthImage.cols-offsetx)/sample1;j++)
		{
                        float depthValue = (float)depthImage.at<float>(sample1*(i+offsety),sample1*(j+offsetx))*0.819;
			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)rgbImage.at<Vec4b>(sample1*i,sample1*j)[3];
			if (avalue > 0 && depthValue>0)                // if depthValue is not NaN
			{
				// Find 3D position respect to rgb frame:
				newPoint.z = depthValue;
				newPoint.x = (sample1*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
				newPoint.y = (sample1*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
				newPoint.r = rgbImage.at<cv::Vec4b>(sample1*i,sample1*j)[2];
				newPoint.g = rgbImage.at<cv::Vec4b>(sample1*i,sample1*j)[1];
				newPoint.b = rgbImage.at<cv::Vec4b>(sample1*i,sample1*j)[0];
				outputPointcloud1->points.push_back(newPoint);
				
				
			}

		}
	}
		
	targetPointCloud =outputPointcloud1;
		
        }

        if (useCurvature.getValue())
        {
        int sample1 = samplePCD.getValue();//3;
        pcl::PointCloud<pcl::PointXYZ>::Ptr outputPointcloud1(new pcl::PointCloud<pcl::PointXYZ>);
        outputPointcloud1->points.resize(0);

        pcl::PointXYZ newPoint1;

        for (int i=0;i<(int)(depthImage.rows-offsety)/sample1;i++)
        {
                for (int j=0;j<(int)(depthImage.cols-offsetx)/sample1;j++)
                {
                        float depthValue = (float)depthImage.at<float>(sample1*(i+offsety),sample1*(j+offsetx))*0.819;
                        //depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
                        int avalue = (int)rgbImage.at<Vec4b>(sample1*i,sample1*j)[3];
                        if (avalue > 0 && depthValue>0)                // if depthValue is not NaN
                        {
                                // Find 3D position respect to rgb frame:
                                newPoint1.z = depthValue;
                                newPoint1.x = (sample1*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
                                newPoint1.y = (sample1*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
                                outputPointcloud1->points.push_back(newPoint1);


                        }

                }
        }


        // Compute the normals
          pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;
          normalEstimation.setInputCloud (outputPointcloud1);
          pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
          normalEstimation.setSearchMethod (tree);
          pcl::PointCloud<pcl::Normal>::Ptr cloudWithNormals (new pcl::PointCloud<pcl::Normal>);
          normalEstimation.setRadiusSearch (0.02);
          normalEstimation.compute (*cloudWithNormals);


          /*vector<pcl::PointXYZ> normals_for_gpu(cloudWithNormals->points.size());
          std::transform(cloudWithNormals->points.begin(), cloudWithNormals->points.end(), normals_for_gpu.begin(), pcl::gpu::DataSource::Normal2PointXYZ());
          pcl::gpu::PrincipalCurvaturesEstimation::PointCloud cloud_gpu;
              cloud_gpu.upload(outputPointcloud1->points);
              pcl::gpu::PrincipalCurvaturesEstimation::Normals normals_gpu;
              normals_gpu.upload(normals_for_gpu);
              pcl::gpu::DeviceArray<pcl::PrincipalCurvatures> pc_features;
              pcl::gpu::PrincipalCurvaturesEstimation pc_gpu;
              pc_gpu.setInputCloud(cloud_gpu);
              pc_gpu.setInputNormals(normals_gpu);
              pc_gpu.setRadiusSearch(0.03, 50);
              pc_gpu.compute(pc_features);
              vector<pcl::PrincipalCurvatures> downloaded;
              pc_features.download(downloaded);
              pcl::PrincipalCurvaturesEstimation<pcl::PointXYZ, pcl::Normal, pcl::PrincipalCurvatures> fe;
              fe.setInputCloud (outputPointcloud1);
              fe.setInputNormals (cloudWithNormals);
              fe.setRadiusSearch(1);
              pcl::PointCloud<pcl::PrincipalCurvatures> pc;
              fe.compute (pc);
              for(size_t i = 0; i < downloaded.size(); ++i)
              {
                  pcl::PrincipalCurvatures& gpu = downloaded[i];
                  pcl::PrincipalCurvatures& cpu = pc.points[i];
              }*/

          // Setup the principal curvatures computation
          pcl::PrincipalCurvaturesEstimation<pcl::PointXYZ, pcl::Normal, pcl::PrincipalCurvatures> principalCurvaturesEstimation;

          // Provide the original point cloud (without normals)
          principalCurvaturesEstimation.setInputCloud (outputPointcloud1);

          // Provide the point cloud with normals
          principalCurvaturesEstimation.setInputNormals(cloudWithNormals);

          // Use the same KdTree from the normal estimation
          principalCurvaturesEstimation.setSearchMethod (tree);
          principalCurvaturesEstimation.setRadiusSearch(0.03);

          // Actually compute the principal curvatures
          pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principalCurvatures (new pcl::PointCloud<pcl::PrincipalCurvatures> ());
          principalCurvaturesEstimation.compute (*principalCurvatures);

          std::cout << "output points.size (): " << principalCurvatures->points.size () << std::endl;

          // Display and retrieve the shape context descriptor vector for the 0th point.
          pcl::PrincipalCurvatures descriptor = principalCurvatures->points[0];

          std::vector<double> curvs;

          for (int k = 0; k < principalCurvatures->points.size (); k++)
          {
              pcl::PrincipalCurvatures descriptor0 = principalCurvatures->points[k];

              double curv = abs(descriptor0.pc1*descriptor0.pc1);
              curvs.push_back(curv);
          }

          curvatures.setValue(curvs);

}

        if (useSIFT3D.getValue())
        {

        int sample1 = samplePCD.getValue();//3;
        pcl::PointCloud<pcl::PointXYZ>::Ptr outputPointcloud1(new pcl::PointCloud<pcl::PointXYZ>);
        outputPointcloud1->points.resize(0);

        pcl::PointXYZ newPoint1;

        for (int i=0;i<(int)(depthImage.rows-offsety)/sample1;i++)
        {
                for (int j=0;j<(int)(depthImage.cols-offsetx)/sample1;j++)
                {
                        float depthValue = (float)depthImage.at<float>(sample1*(i+offsety),sample1*(j+offsetx))*0.819;
                        //depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
                        int avalue = (int)rgbImage.at<Vec4b>(sample1*i,sample1*j)[3];
                        if (avalue > 0 && depthValue>0)                // if depthValue is not NaN
                        {
                                // Find 3D position respect to rgb frame:
                                newPoint1.z = depthValue;
                                newPoint1.x = (sample1*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
                                newPoint1.y = (sample1*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
                                outputPointcloud1->points.push_back(newPoint1);


                        }

                }
        }

        // Compute the normals
          pcl::NormalEstimation<pcl::PointXYZ, pcl::PointNormal> normalEstimation;
          normalEstimation.setInputCloud (outputPointcloud1);
          pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
          normalEstimation.setSearchMethod (tree);
          pcl::PointCloud<pcl::PointNormal>::Ptr cloudWithNormals (new pcl::PointCloud<pcl::PointNormal>);
          normalEstimation.setRadiusSearch (0.02);
          normalEstimation.compute (*cloudWithNormals);

  // Parameters for sift computation
  const float min_scale = 0.01f;
  const int n_octaves = 3;
  const int n_scales_per_octave = 4;
  const float min_contrast = 0.001f;

  // Copy the xyz info from cloud_xyz and add it to cloud_normals as the xyz field in PointNormals estimation is zero
  for(size_t i = 0; i<cloudWithNormals->points.size(); ++i)
  {
    cloudWithNormals->points[i].x = outputPointcloud1->points[i].x;
    cloudWithNormals->points[i].y = outputPointcloud1->points[i].y;
    cloudWithNormals->points[i].z = outputPointcloud1->points[i].z;
  }

  // Estimate the sift interest points using normals values from xyz as the Intensity variants
  pcl::SIFTKeypoint<pcl::PointNormal, pcl::PointWithScale> sift;
  pcl::PointCloud<pcl::PointWithScale> result;
  pcl::search::KdTree<pcl::PointNormal>::Ptr treenormal(new pcl::search::KdTree<pcl::PointNormal> ());
  sift.setSearchMethod(treenormal);
  sift.setScales(min_scale, n_octaves, n_scales_per_octave);
  sift.setMinimumContrast(min_contrast);
  sift.setInputCloud(cloudWithNormals);
  sift.compute(result);

  pcl::Feature<pcl::PointXYZRGB, pcl::FPFHSignature33>::Ptr feature_extractor (new pcl::FPFHEstimationOMP<pcl::PointXYZRGB, pcl::Normal, pcl::FPFHSignature33>);
  feature_extractor->setSearchMethod (pcl::search::Search<pcl::PointXYZRGB>::Ptr (new pcl::search::KdTree<pcl::PointXYZRGB>));
  feature_extractor->setRadiusSearch (0.05);

  typename pcl::PointCloud<pcl::PointXYZRGB>::Ptr kpts(new pcl::PointCloud<pcl::PointXYZRGB>);
  kpts->points.resize(result.points.size());

  pcl::copyPointCloud(result, *kpts);

  feature_extractor->setSearchSurface(outputPointcloud);
  feature_extractor->setInputCloud(kpts);

  /*typename pcl::PointCloud<pcl::FPFHSignature33>::Ptr features(new pcl::PointCloud<pcl::FPFHSignature33>);

  typename pcl::FeatureFromNormals<pcl::PointXYZRGB, pcl::Normal, FeatureType>::Ptr feature_from_normals = boost::dynamic_pointer_cast<pcl::FeatureFromNormals<pcl::PointXYZRGB, pcl::Normal, FeatureType> > (feature_extractor_);

  cout << "descriptor extraction..." << std::flush;
  feature_extractor->compute (*features);
  cout << "OK" << endl;*/

  std::cout << "No of SIFT points in the result are " << result.points.size () << std::endl;
        }

		
	return outputPointcloud;
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::computeCenter(vpImage<unsigned char> &Itemp, vpImagePoint &cog,double &angle, int &surface)
{
vpImagePoint ip;
int it;
it=0;
/*for (int n=0; n < Itemp.getHeight() ; n++)
{
for (int m=0 ; m < Itemp.getWidth(); m++)
  {
//std::cout << "Dot : " << (int)Itemp[n][m] << std::endl;
	if((int)Itemp[n][m]!=16)
    {Itemp[n][m]=255;
    if(it==1500){
	ip.set_i(n);
	ip.set_j(m);}
    it++;
    }
	else{
		Itemp[n][m]=0;
	}
  }
}*/
double mu11,mu20,mu02,mup11,mup20,mup02,m00,m01,m10,m11,m02,m20;
surface = 0;
angle = 0;
m00=0;
m01 = 0;
m10 =0 ;
m11 = 0;
m20 = 0;
m02 = 0;
mu11 = 0;
mu20=0;
mu02=0;

for (int n=0; n < Itemp.getHeight()-0 ; n++)
{
for (int m=0 ; m < Itemp.getWidth()-0; m++)
  {
	if(Itemp[n][m] == 255)
	{
m00 ++;
m01 += n;
m10 += m;
m11 += m*n;
m02 += n*n;
m20 += m*m;
	}
  }
}
cog.set_u(m10/m00);
cog.set_v(m01/m00);

surface = m00;

mu11 = m11-m10*m01/m00;
mu20 = m20 - m10*m10/m00;
mu02 = m02 - m01*m01/m00;

mup11 = mu11/m00;
mup02 = mu02/m00;
mup20 = mu20/m00;

angle =0.5*atan2(2*mup11,(mup20-mup02));

std::cout << "Dot : " << " COG : " << cog.get_u() << " " << cog.get_v() << " Angle : " << angle << " Surface : " << surface << std::endl;
}

template <class DataTypes>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr RGBDDataProcessing<DataTypes>::PCDContourFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage, cv::Mat& distImage, cv::Mat& dotImage)
{
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    //pcl::PointCloud<pcl::PointXYZRGB> pointcloud;
	
    cv::Mat frgd;
    cv::Mat distimg,dotimg;
    distimg = distImage;
    dotimg = dotImage;

    helper::vector<bool> targetborder;
    targetborder.resize(0);

    helper::vector<double> targetweights;
    targetweights.resize(0);
    int sample;
    int offsetx;

    frgd = rgbImage;
    sample = samplePCD.getValue();//3;
    offsetx = offsetX.getValue();//0;

    float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);// 1/fx
    float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);// 1/fy
    pcl::PointXYZRGB newPoint;
    ntargetcontours = 0;
    int jj = 0;
    int offsety = offsetY.getValue();
    double totalweights = 0;

	for (int i=0;i<(int)(depthImage.rows-offsety)/sample;i++)
	{
            for (int j=0;j<(int)(depthImage.cols-offsetx)/sample;j++)
            {
                float depthValue = (float)depthImage.at<float>(sample*(i+offsety),sample*(j+offsetx))*0.819;
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

                    targetweights.push_back((double)exp(-bvalue/sigmaWeight.getValue()));
                    totalweights += targetweights[jj];

                    jj++;

                    if (avalue > 0 && bvalue < borderThdPCD.getValue() /*4*/)
                    {
                        targetborder.push_back(true);
                        ntargetcontours++;
                    }
                    else targetborder.push_back(false);
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
	
        for (int i=0; i < targetweights.size();i++)
	{
            targetweights[i]*=((double)targetweights.size()/totalweights);
            //std::cout << " weights " << totalweights << " " << (double)targetweights[i] << std::endl;
	}

        targetWeights.setValue(targetweights);
        targetBorder.setValue(targetborder);

        if (useCurvature.getValue())
        {
        int sample1 = samplePCD.getValue();//3;
        pcl::PointCloud<pcl::PointXYZ>::Ptr outputPointcloud1(new pcl::PointCloud<pcl::PointXYZ>);
        outputPointcloud1->points.resize(0);

        pcl::PointXYZ newPoint1;

        for (int i=0;i<(int)(depthImage.rows-offsety)/sample1;i++)
        {
                for (int j=0;j<(int)(depthImage.cols-offsetx)/sample1;j++)
                {
                        float depthValue = (float)depthImage.at<float>(sample1*(i+offsety),sample1*(j+offsetx))*0.819;
                        //depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
                        int avalue = (int)rgbImage.at<Vec4b>(sample1*i,sample1*j)[3];
                        if (avalue > 0 && depthValue>0)                // if depthValue is not NaN
                        {
                                // Find 3D position respect to rgb frame:
                                newPoint1.z = depthValue;
                                newPoint1.x = (sample1*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
                                newPoint1.y = (sample1*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
                                outputPointcloud1->points.push_back(newPoint1);


                        }

                }
        }


        // Compute the normals
          pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimation;
          normalEstimation.setInputCloud (outputPointcloud1);

          pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
          normalEstimation.setSearchMethod (tree);

          pcl::PointCloud<pcl::Normal>::Ptr cloudWithNormals (new pcl::PointCloud<pcl::Normal>);

          normalEstimation.setRadiusSearch (0.02);

          normalEstimation.compute (*cloudWithNormals);


          /*vector<pcl::PointXYZ> normals_for_gpu(cloudWithNormals->points.size());
          std::transform(cloudWithNormals->points.begin(), cloudWithNormals->points.end(), normals_for_gpu.begin(), pcl::gpu::DataSource::Normal2PointXYZ());


          pcl::gpu::PrincipalCurvaturesEstimation::PointCloud cloud_gpu;
              cloud_gpu.upload(outputPointcloud1->points);

              pcl::gpu::PrincipalCurvaturesEstimation::Normals normals_gpu;
              normals_gpu.upload(normals_for_gpu);

              pcl::gpu::DeviceArray<pcl::PrincipalCurvatures> pc_features;

              pcl::gpu::PrincipalCurvaturesEstimation pc_gpu;
              pc_gpu.setInputCloud(cloud_gpu);
              pc_gpu.setInputNormals(normals_gpu);
              pc_gpu.setRadiusSearch(0.03, 50);
              pc_gpu.compute(pc_features);

              vector<pcl::PrincipalCurvatures> downloaded;
              pc_features.download(downloaded);

              pcl::PrincipalCurvaturesEstimation<pcl::PointXYZ, pcl::Normal, pcl::PrincipalCurvatures> fe;
              fe.setInputCloud (outputPointcloud1);
              fe.setInputNormals (cloudWithNormals);
              fe.setRadiusSearch(1);

              pcl::PointCloud<pcl::PrincipalCurvatures> pc;
              fe.compute (pc);

              for(size_t i = 0; i < downloaded.size(); ++i)
              {
                  pcl::PrincipalCurvatures& gpu = downloaded[i];
                  pcl::PrincipalCurvatures& cpu = pc.points[i];

              }*/

          // Setup the principal curvatures computation
          pcl::PrincipalCurvaturesEstimation<pcl::PointXYZ, pcl::Normal, pcl::PrincipalCurvatures> principalCurvaturesEstimation;

          // Provide the original point cloud (without normals)
          principalCurvaturesEstimation.setInputCloud (outputPointcloud1);

          // Provide the point cloud with normals
          principalCurvaturesEstimation.setInputNormals(cloudWithNormals);

          // Use the same KdTree from the normal estimation
          principalCurvaturesEstimation.setSearchMethod (tree);
          principalCurvaturesEstimation.setRadiusSearch(0.03);

          // Actually compute the principal curvatures
          pcl::PointCloud<pcl::PrincipalCurvatures>::Ptr principalCurvatures (new pcl::PointCloud<pcl::PrincipalCurvatures> ());
          principalCurvaturesEstimation.compute (*principalCurvatures);

          std::cout << "output points.size (): " << principalCurvatures->points.size () << std::endl;

          // Display and retrieve the shape context descriptor vector for the 0th point.
          pcl::PrincipalCurvatures descriptor = principalCurvatures->points[0];

          std::vector<double> curvs;

          for (int k = 0; k < principalCurvatures->points.size (); k++)
          {
              pcl::PrincipalCurvatures descriptor0 = principalCurvatures->points[k];

              double curv = abs(descriptor0.pc1*descriptor0.pc1);
              curvs.push_back(curv);
          }

          curvatures.setValue(curvs);
          std::cout << " curvature " << descriptor << std::endl;

        }

        if (useSIFT3D.getValue())
        {

        int sample1 = samplePCD.getValue();//3;
        pcl::PointCloud<pcl::PointXYZ>::Ptr outputPointcloud1(new pcl::PointCloud<pcl::PointXYZ>);
        outputPointcloud1->points.resize(0);

        pcl::PointXYZ newPoint1;

        for (int i=0;i<(int)(depthImage.rows-offsety)/sample1;i++)
        {
                for (int j=0;j<(int)(depthImage.cols-offsetx)/sample1;j++)
                {
                        float depthValue = (float)depthImage.at<float>(sample1*(i+offsety),sample1*(j+offsetx))*0.819;
                        //depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
                        int avalue = (int)rgbImage.at<Vec4b>(sample1*i,sample1*j)[3];
                        if (avalue > 0 && depthValue>0)                // if depthValue is not NaN
                        {
                                // Find 3D position respect to rgb frame:
                                newPoint1.z = depthValue;
                                newPoint1.x = (sample1*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
                                newPoint1.y = (sample1*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
                                outputPointcloud1->points.push_back(newPoint1);


                        }

                }
        }

        // Compute the normals
          pcl::NormalEstimation<pcl::PointXYZ, pcl::PointNormal> normalEstimation;
          normalEstimation.setInputCloud (outputPointcloud1);
          pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
          normalEstimation.setSearchMethod (tree);
          pcl::PointCloud<pcl::PointNormal>::Ptr cloudWithNormals (new pcl::PointCloud<pcl::PointNormal>);
          normalEstimation.setRadiusSearch (0.02);
          normalEstimation.compute (*cloudWithNormals);

  // Parameters for sift computation
  const float min_scale = 0.01f;
  const int n_octaves = 3;
  const int n_scales_per_octave = 4;
  const float min_contrast = 0.001f;

  // Copy the xyz info from cloud_xyz and add it to cloud_normals as the xyz field in PointNormals estimation is zero
  for(size_t i = 0; i<cloudWithNormals->points.size(); ++i)
  {
    cloudWithNormals->points[i].x = outputPointcloud1->points[i].x;
    cloudWithNormals->points[i].y = outputPointcloud1->points[i].y;
    cloudWithNormals->points[i].z = outputPointcloud1->points[i].z;
  }

  // Estimate the sift interest points using normals values from xyz as the Intensity variants
  pcl::SIFTKeypoint<pcl::PointNormal, pcl::PointWithScale> sift;
  pcl::PointCloud<pcl::PointWithScale> result;
  pcl::search::KdTree<pcl::PointNormal>::Ptr treenormal(new pcl::search::KdTree<pcl::PointNormal> ());
  sift.setSearchMethod(treenormal);
  sift.setScales(min_scale, n_octaves, n_scales_per_octave);
  sift.setMinimumContrast(min_contrast);
  sift.setInputCloud(cloudWithNormals);
  sift.compute(result);



std::cout << "No of SIFT points in the result are " << result.points.size () << std::endl;
        }

	/*const std::string file = "test_pcdf%03d.pcd";
    char buf[FILENAME_MAX];
    sprintf(buf, file.c_str(), frame);
    std::string filename(buf);*/
	//pcl::io::savePCDFileASCII (filename, outputPointcloud);
	//std::cout << "Saved " << outputPointcloud.points.size () << " data points to test_pcd.pcd." << std::endl;
	
	return outputPointcloud;
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::ContourFromRGBSynth(cv::Mat& rgbImage, cv::Mat& distImage, cv::Mat& dotImage)
{
	
	ntargetcontours = 0;	
	cv::Mat frgd;
	cv::Mat distimg,dotimg;
	
	//cv::resize(rgbImage, frgd, cv::Size(foreground.cols/2, foreground.rows/2));
	distimg = distImage.clone();
	dotimg = dotImage.clone();
	frgd = rgbImage.clone();
	
        helper::vector<bool> targetborder;
        targetborder.resize(0);
        helper::vector<double> targetweights;
        targetweights.resize(0);
	const VecCoord& targetp = targetPositions.getValue();
	
	int nvisible = 0;
	double totalweights = 0;

	for (int k = 0; k < targetp.size(); k++)
	{
	int x_u = (int)(targetp[k][0]*rgbIntrinsicMatrix(0,0)/targetp[k][2] + rgbIntrinsicMatrix(0,2));
        int x_v = (int)(targetp[k][1]*rgbIntrinsicMatrix(1,1)/targetp[k][2] + rgbIntrinsicMatrix(1,2));

	int avalue = (int)frgd.at<Vec4b>(x_v,x_u)[3];
	int bvalue = (int)distimg.at<uchar>(x_v,x_u);
	int dvalue = (int)dotimg.at<uchar>(x_v,x_u);
		{
								
                                targetweights.push_back((double)exp(-bvalue/7));
                                totalweights += targetweights[k];
								
                                if (bvalue < borderThdPCD.getValue() /*4*/) {targetborder.push_back(true);ntargetcontours++;}
                                        else targetborder.push_back(false);
		}
	
}
	
for (int i=0; i < targetweights.size();i++)
	{
                targetweights[i]*=((double)targetweights.size()/totalweights);
		//std::cout << " weights " << (double)targetWeights[i] << std::endl;
	}
	
VecCoord targetContourpos;
targetContourpos.resize(ntargetcontours);
std::cout << " ntargetcontours " << ntargetcontours << std::endl;
int kk = 0;
for (int i = 0; i < targetborder.size(); i++)
	{
                if (targetborder[i]){
            targetContourpos[kk]=targetp[i];
			kk++;
		}
	}
    const VecCoord&  p1 = targetContourpos;
    targetContourPositions.setValue(p1);
    targetWeights.setValue(targetweights);
    targetBorder.setValue(targetborder);
		
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::extractTargetPCD()
{
        targetP.reset(new pcl::PointCloud<pcl::PointXYZRGB>); 
	targetP = PCDFromRGBD(depth,foreground);
	VecCoord targetpos;
	
	if (targetP->size() > 10)
	{
	target.reset(new pcl::PointCloud<pcl::PointXYZRGB>);	
	target = targetP;
	targetpos.resize(target->size());

            Vector3 pos;
            Vector3 col;

	for (unsigned int i=0; i<target->size(); i++)
	{
            pos[0] = (double)target->points[i].x;
            pos[1] = (double)target->points[i].y;
            pos[2] = (double)target->points[i].z;
            targetpos[i]=pos;
            //std::cout << " target " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	} 
    const VecCoord&  p = targetpos;
    targetPositions.setValue(p);
	
}

}


template <class DataTypes>
void RGBDDataProcessing<DataTypes>::extractTargetPCDContour()
{
	
    targetP.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
		
    double cannyTh1 = 150;
    double cannyTh2 = 80;
    cv::Mat contour,dist,dist0;
	
    //cv::imwrite("depthmap.png", seg.distImage);
    cv::Canny( dotimage, contour, cannyTh1, cannyTh2, 3);
    contour = cv::Scalar::all(255) - contour;
    cv::distanceTransform(contour, dist, CV_DIST_L2, 3);
    dist.convertTo(dist0, CV_8U, 1, 0);	
    seg.distImage = dist0.clone();
    targetP = PCDContourFromRGBD(depth,foreground, distimage, dotimage);
    VecCoord targetpos;
	
	if (targetP->size() > 10)
	{
            target.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
            target = targetP;
            targetpos.resize(target->size());

            Vector3 pos;
                for (unsigned int i=0; i<target->size(); i++)
                {
                    pos[0] = (double)target->points[i].x;
                    pos[1] = (double)target->points[i].y;
                    pos[2] = (double)target->points[i].z;
                    targetpos[i]=pos;
                }
            VecCoord targetContourpos;
            targetContourpos.resize(ntargetcontours);
            int kk = 0;
                for (unsigned int i=0; i<target->size(); i++)
                {
                    if (targetBorder.getValue()[i])
                    {
                        pos[0] = (double)target->points[i].x;
                        pos[1] = (double)target->points[i].y;
                        pos[2] = (double)target->points[i].z;
                        targetContourpos[kk]=pos;
                        kk++;
                    }
                }

            const VecCoord&  p0 = targetpos;
            targetPositions.setValue(p0);
            const VecCoord&  p1 = targetContourpos;
            targetContourPositions.setValue(p1);
            std::cout << " target contour " << p1.size() << std::endl;
	}

}

template<class DataTypes>
void RGBDDataProcessing<DataTypes>::setCameraPose()
{
    pcl::PointCloud<pcl::PointXYZRGB>& point_cloud = *target;
    if (target->size() > 0)
    {
    Vec3 cameraposition;
    Quat cameraorientation;
    cameraposition[0] = point_cloud.sensor_origin_[0];
    cameraposition[1] = point_cloud.sensor_origin_[1];
    cameraposition[2] = point_cloud.sensor_origin_[2];
    cameraorientation[0] = point_cloud.sensor_orientation_.w ();
    cameraorientation[1] = point_cloud.sensor_orientation_.x ();
    cameraorientation[2] = point_cloud.sensor_orientation_.y ();
    cameraorientation[3] = point_cloud.sensor_orientation_.z ();

    cameraPosition.setValue(cameraposition);
    cameraOrientation.setValue(cameraorientation);
    cameraChanged.setValue(true);
    }
}


template <class DataTypes>
void RGBDDataProcessing<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{
        if (dynamic_cast<simulation::AnimateBeginEvent*>(event))
	{
	double timeT = (double)getTickCount();
        double timeAcq0 = (double)getTickCount();

        sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
        typename sofa::core::objectmodel::ImageConverter<DataTypes,DepthTypes>::SPtr imconv;
        root->get(imconv);
	
        typename sofa::core::objectmodel::DataIO<DataTypes>::SPtr dataio;
	root->get(dataio);

        bool okimages =false;
        bool newimages = false;
	
	if (useRealData.getValue())
	{
        if (useSensor.getValue()){
                //color_1 = color.clone();
                //depth_1 = depth.clone();
            if(!((imconv->depth).empty()) && !((imconv->color).empty()))
            {
		if (scaleImages.getValue() > 1)
		{	
                cv::resize(imconv->depth, depth, cv::Size(imconv->depth.cols/scaleImages.getValue(), imconv->depth.rows/scaleImages.getValue()), 0, 0);
                cv::resize(imconv->color, color, cv::Size(imconv->color.cols/scaleImages.getValue(), imconv->color.rows/scaleImages.getValue()), 0, 0);
		}
		else
		{
		color = imconv->color;
		depth = imconv->depth;
                }
                okimages = true;
                newimages=imconv->newImages.getValue();
            }
                //cv::imwrite("depth22.png", depth);
	}
        else {
		color = dataio->color;
		depth = dataio->depth;
                color_1 = dataio->color_1;
                okimages = true;
                newimages=dataio->newImages.getValue();
	 }

        double timeAcq1 = (double)getTickCount();
        cout <<"TIME GET IMAGES " << (timeAcq1 - timeAcq0)/getTickFrequency() << endl;

        std::cout << "newimages " << newimages << std::endl;

        imagewidth.setValue(color.cols);
        imageheight.setValue(color.rows);


        if (displayImages.getValue() && displayDownScale.getValue() > 0 && !depth.empty() && !color.empty())
        {
        int scale = displayDownScale.getValue();
        cv::Mat colorS, depthS;
        cv::resize(depth, depthS, cv::Size(imagewidth.getValue()/scale, imageheight.getValue()/scale), 0, 0);
        cv::resize(color, colorS, cv::Size(imagewidth.getValue()/scale, imageheight.getValue()/scale), 0, 0);

        /*cv::Mat depthmat1;
        depthS.convertTo(depthmat1, CV_8UC1, 255);
        cv::imwrite("depthS0.png", depthmat1);*/
        cv::imshow("image_sensor",colorS);
        cv::waitKey(1);
        /*cv::imshow("depth_sensor",depthS);
        cv::waitKey(1);*/
        }

        if (saveImages.getValue())
	{
                cv::Mat* imgl = new cv::Mat;
                *imgl = color.clone();
                cv::Mat* depthl = new cv::Mat;
                *depthl = depth.clone();

                dataio->listimg.push_back(imgl);
                dataio->listdepth.push_back(depthl);
	}
	}
        else
        {
            dataio->readData();
            okimages=true;
        }


        if (okimages && newimages)
        {
        if (initsegmentation)
	{

	if (useRealData.getValue())
	{
		initSegmentation();
                extractTargetPCD();
        }
        setCameraPose();
        initsegmentation = false;
	}
	else
        {
            if(useRealData.getValue() && !stopatinit.getValue())
            {
            segment() ;
            timePCD = (double)getTickCount();

            if(!useContour.getValue())
            extractTargetPCD();
            else extractTargetPCDContour();

            timePCD = ((double)getTickCount() - timePCD)/getTickFrequency();
            std::cout << "TIME PCD " << timePCD << std::endl;

            }
            else if (!stopatinit.getValue()){
            segmentSynth();
            if(useContour.getValue())
            ContourFromRGBSynth(foreground, distimage,dotimage);
            }

            if (saveImages.getValue())
            {
            cv::Mat* imgseg = new cv::Mat;
            *imgseg = foreground.clone();
            dataio->listimgseg.push_back(imgseg);
            }
            cameraChanged.setValue(false);
        }
        }

        std::cout << "TIME RGBDDATAPROCESSING " << ((double)getTickCount() - timeT)/getTickFrequency() << std::endl;


}
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

    ReadAccessor< Data< VecCoord > > xtarget(targetPositions);
    vparams->drawTool()->saveLastState();

    if (displayBackgroundImage.getValue())
    {
    GLfloat projectionMatrixData[16];
    glGetFloatv(GL_PROJECTION_MATRIX, projectionMatrixData);
    GLfloat modelviewMatrixData[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, modelviewMatrixData);

    cv::Mat colorrgb = color.clone();
    if (!color.empty())
    cv::cvtColor(color, colorrgb, CV_RGB2BGR);

    std::stringstream imageString;
    imageString.write((const char*)colorrgb.data, colorrgb.total()*colorrgb.elemSize());
    // PERSPECTIVE

    glMatrixMode(GL_PROJECTION);	//init the projection matrix
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, -1, 1);  // orthogonal view
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // BACKGROUND TEXTURING
    //glDepthMask (GL_FALSE);		// disable the writing of zBuffer
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);	// enable the texture
    glDisable(GL_LIGHTING);		// disable the light

    glBindTexture(GL_TEXTURE_2D, 0);  // texture bind
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, colorrgb.cols, colorrgb.rows, 0, GL_RGB, GL_UNSIGNED_BYTE, imageString.str().c_str());

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);	// Linear Filtering

                                                                        // BACKGROUND DRAWING
                                                                        //glEnable(GL_DEPTH_TEST);

    glBegin(GL_QUADS); //we draw a quad on the entire screen (0,0 1,0 1,1 0,1)
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glTexCoord2f(0, 1);		glVertex2f(0, 0);
    glTexCoord2f(1, 1);		glVertex2f(1, 0);
    glTexCoord2f(1, 0);		glVertex2f(1, 1);
    glTexCoord2f(0, 0);		glVertex2f(0, 1);
    glEnd();

    //glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);		// enable light
    glDisable(GL_TEXTURE_2D);	// disable texture 2D
    glEnable(GL_DEPTH_TEST);
    //glDepthMask (GL_TRUE);		// enable zBuffer

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);


    vparams->drawTool()->restoreLastState();
    }

    if (drawPointCloud.getValue() && xtarget.size() > 0){

                    std::vector< sofa::defaulttype::Vector3 > points;
                    sofa::defaulttype::Vector3 point;

      for (unsigned int i=0; i< xtarget.size(); i++)
        {
          points.resize(0);
          point = DataTypes::getCPos(xtarget[i]);
          points.push_back(point);
         // std::cout << curvatures.getValue()[i] << std::endl;
          //if (targetWeights.getValue().size()>0) vparams->drawTool()->drawPoints(points, 10, sofa::defaulttype::Vec<4,float>(0.5*targetWeights.getValue()[i],0,0,1));
         vparams->drawTool()->drawPoints(points, 10, sofa::defaulttype::Vec<4,float>(1,0.5,0.5,1));
        }

    }

}



}
}
} // namespace sofa

