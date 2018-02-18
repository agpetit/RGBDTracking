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

#include <pcl/common/common_headers.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/impl/point_types.hpp>

#include <sofa/helper/gl/Color.h>
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseVisual/InteractiveCamera.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/simulation/Simulation.h>

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
	, niterations(initData(&niterations,3,"niterations","Number of iterations in the tracking process"))
	, nimages(initData(&nimages,1500,"nimages","Number of images to read"))
	, samplePCD(initData(&samplePCD,4,"samplePCD","Sample step for the point cloud"))
	, offsetX(initData(&offsetX,3,"offsetX","offset along x for the point cloud"))
	, offsetY(initData(&offsetY,0,"offsetY","offset along y for the point cloud"))
	, borderThdPCD(initData(&borderThdPCD,4,"borderThdPCD","border threshold on the target silhouette"))
	, inputPath(initData(&inputPath,"inputPath","Path for data readings",false))
	, outputPath(initData(&outputPath,"outputPath","Path for data writings",false))
	, useDistContourNormal(initData(&useDistContourNormal,false,"outputPath","Path for data writings"))
	, windowKLT(initData(&windowKLT,1500,"nimages","Number of images to read"))
 {
	this->f_listening.setValue(true); 
	iter_im = 0;
    // initialize paramters
    //

    // Tracker parameters
    tracker.setTrackerId(1);
    //tracker.setOnMeasureFeature(&modifyFeature);
    tracker.setMaxFeatures(200);
    tracker.setWindowSize(10);
    tracker.setQuality(0.01);
    tracker.setMinDistance(10);
    tracker.setHarrisFreeParameter(0.04);
    tracker.setBlockSize(9);
    tracker.setUseHarris(1);
    tracker.setPyramidLevels(3); 
	
	tracker1.setTrackerId(1);
    //tracker.setOnMeasureFeature(&modifyFeature);
    tracker1.setMaxFeatures(200);
    tracker1.setWindowSize(10);
    tracker1.setQuality(0.01);
    tracker1.setMinDistance(10);
    tracker1.setHarrisFreeParameter(0.04);
    tracker1.setBlockSize(9);
    tracker1.setUseHarris(1);
    tracker1.setPyramidLevels(3);
	
	bool g_bQATest = false;
	int  g_nDevice = 0;
}

template <class DataTypes>
RGBDDataProcessing<DataTypes>::~RGBDDataProcessing()
{
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::init()
{

	std::string configFile;	
    bool use_cuda = true;
		
	if (configFile.empty())
	configFile = vpIoTools::path("param/pizza.lua");

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glEnable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        glDepthMask(GL_TRUE);

		
    seg.init(configFile);
	
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

	Vector4 camParam = cameraIntrinsicParameters.getValue();
	
        rgbIntrinsicMatrix(0,0) = camParam[0];
	rgbIntrinsicMatrix(1,1) = camParam[1];
	rgbIntrinsicMatrix(0,2) = camParam[2];
	rgbIntrinsicMatrix(1,2) = camParam[3];
	
	cv::Mat mask,maskimg,mask0,roimask,mask1; // segmentation result (4 possible values)
    cv::Mat bgModel,fgModel; // the models (internally used)

	//cv::imwrite("color.png",color);
	
	cv::Mat downsampledbox,downsampled;
	
	//cv::pyrDown(color, downsampledbox, cv::Size(color.cols/2, color.rows/2));
	downsampledbox = color;
	
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
	
		/*for (int i = 0; i<color.rows; i++)
		  for (int j = 0; j<color.cols; j++)
	  std::cout << (int)color.at<Vec3b>(i,j)[0] << std::endl;*/
	  
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
      int key=cvWaitKey(10);
      if ((char)key == 27) break;
	
    }
   // delete temp1;
    temp1.release();
    cv::setMouseCallback(name, NULL, NULL);
	
	    //cvDestroyWindow(name);
    cv::Rect rectangle(69,47,198,171);
    cv::Rect recttemp = rectangle;

    rectangle = box;
    seg.setRectangle(rectangle);

    int width = downsampled.cols;
    int height = downsampled.rows;

    foreground = cv::Mat(downsampled.size(),CV_8UC3,cv::Scalar(255,255,255));
		//cv::pyrDown(color, downsampled, cv::Size(color.cols/2, color.rows/2));
	downsampled = color;

    seg.segmentationFromRect(downsampled,foreground);

    // draw rectangle on original image
    //cv::rectangle(image, rectangle, cv::Scalar(255,255,255),1);
    
    // display result
    cv::namedWindow("Segmented Image");
    cv::imshow("Segmented Image",foreground);

    if (waitKey(20) == 27) //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
   {
        cout << "esc key is pressed by user" << endl;
   }
   
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::segment()
{

        cv::Mat downsampled;
        downsampled = color.clone();
		
	int width = downsampled.cols;
        int height = downsampled.rows;
		
	seg.updateMask(foreground);
        seg.updateSegmentation(downsampled,foreground);
        //seg.updateSegmentationCrop(downsampled,foreground);
	
	//foreground = foreground0.clone();
	
	/*cv::namedWindow("Image");
        cv::imshow("Image",downsampled);*/

        // display result
        cv::namedWindow("image_segmented");
        cv::imshow("image_segmented",foreground);

        cv::waitKey(1);
}

template<class DataTypes>
void RGBDDataProcessing<DataTypes>::segmentSynth()
{

	cv::Mat downsampled;
	//cv::pyrDown(color, downsampled, cv::Size(color.cols/2, color.rows/2));
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
	
  cv::imwrite("downsampled1.png",foreground);

	//distImage = distanceMap;
	//dotImage = dot;
	
	/*cv::namedWindow("Image");
    cv::imshow("Image",downsampled);*/

    // display result
    cv::namedWindow("Segmented Image");
    cv::imshow("Segmented Image",foreground);

}

template <class DataTypes>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr RGBDDataProcessing<DataTypes>::PCDFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	//pcl::PointCloud<pcl::PointXYZRGB> pointcloud;
	
	cv::Mat frgd;
	
	int sample;
	int offsetx;
	
	switch (sensorType.getValue())
    {
    // FLOAT ONE CHANNEL
    case 0:
	frgd = rgbImage;
	sample = 2;
	offsetx = 0;
	break;
	case 1:
	frgd = rgbImage;
	
	//157.8843587, 113.3869574, 260.0061005, 262.0832823
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
			float depthValue = (float)depthImage.at<float>(sample1*(i+offsety),sample1*(j+offsetx));
			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)frgd.at<Vec4b>(sample1*i,sample1*j)[3];
			if (avalue > 0 && depthValue>0)                // if depthValue is not NaN
			{
				// Find 3D position respect to rgb frame:
				newPoint.z = depthValue;
				newPoint.x = (sample1*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
				newPoint.y = (sample1*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
				newPoint.r = frgd.at<cv::Vec4b>(sample1*i,sample1*j)[2];
				newPoint.g = frgd.at<cv::Vec4b>(sample1*i,sample1*j)[1];
				newPoint.b = frgd.at<cv::Vec4b>(sample1*i,sample1*j)[0];
				outputPointcloud1->points.push_back(newPoint);
				
				
			}

		}
	}
		
		targetPointCloud =outputPointcloud1;
		
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
	
	targetBorder.resize(0);
	targetWeights.resize(0);
	int sample;
	int offsetx;
	
	//std::cout << " ok " << std::endl;
	
    switch (sensorType.getValue())
    {
    // FLOAT ONE CHANNEL
    case 0:
	//cv::pyrDown(rgbImage, frgd, cv::Size(foreground.cols/2, foreground.rows/2));
	frgd = rgbImage;

	sample = samplePCD.getValue();//2
	offsetx = offsetX.getValue();//0;
	break;
	case 1:
	frgd = rgbImage;  
	sample = samplePCD.getValue();//3;
	offsetx = offsetX.getValue();//0;
	break;
	}
	
	//std::cout << " ok0 " << std::endl;

	float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
	float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy
	pcl::PointXYZRGB newPoint;
	ntargetcontours = 0;
	int jj = 0;
		int offsety = offsetY.getValue();
	double totalweights = 0;
	for (int i=0;i<(int)(depthImage.rows-offsety)/sample;i++)
	{
		for (int j=0;j<(int)(depthImage.cols-offsetx)/sample;j++)
		{
			float depthValue = (float)depthImage.at<float>(sample*(i+offsety),sample*(j+offsetx));
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
				
				targetWeights.push_back((double)exp(-bvalue/sigmaWeight.getValue()));

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
				
				if (avalue > 0 && bvalue < borderThdPCD.getValue() /*4*/) {targetBorder.push_back(true);ntargetcontours++;}
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
		
	//std::cout << " ok 1" << std::endl;	
	
	return outputPointcloud;
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::ContourFromRGBSynth(cv::Mat& rgbImage, cv::Mat& distImage, cv::Mat& dotImage)
{
	
	ntargetcontours = 0;	
	cv::Mat frgd;
	cv::Mat distimg,dotimg;
	
	//cv::pyrDown(rgbImage, frgd, cv::Size(foreground.cols/2, foreground.rows/2));
	distimg = distImage.clone();
	dotimg = dotImage.clone();
	frgd = rgbImage.clone();
	
	targetBorder.resize(0);
	targetWeights.resize(0);
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
								
				targetWeights.push_back((double)exp(-bvalue/7));
				totalweights += targetWeights[k];
								
				if (bvalue < borderThdPCD.getValue() /*4*/) {targetBorder.push_back(true);ntargetcontours++;}
					else targetBorder.push_back(false);
		}
	
}
	
for (int i=0; i < targetWeights.size();i++)
	{
		targetWeights[i]*=((double)targetWeights.size()/totalweights);
		//std::cout << " weights " << (double)targetWeights[i] << std::endl;
	}
	
std::cout << " " << ntargetcontours << std::endl;

VecCoord targetContourpos;
targetContourpos.resize(ntargetcontours);
std::cout << " ntargetcontours " << ntargetcontours << std::endl;
int kk = 0;
for (int i = 0; i < targetBorder.size(); i++)
	{
		if (targetBorder[i]){
            targetContourpos[kk]=targetp[i];
			kk++;
		}
	}
    const VecCoord&  p1 = targetContourpos;
	targetContourPositions.setValue(p1);
		
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
    //targetKdTree.build(p);

	//targetKdTree = sourceKdTree;
	
    // detect border
    //if(targetBorder.size()!=p.size()) { targetBorder.resize(p.size()); detectBorder(targetBorder,targetTriangles.getValue()); }

}
template <class DataTypes>
void RGBDDataProcessing<DataTypes>::extractTargetPCDContour()
{
	
	targetP.reset(new pcl::PointCloud<pcl::PointXYZRGB>); 
		
	double cannyTh1 = 150;
	double cannyTh2 = 80;
	cv::Mat contour,dist,dist0;
	
	//cv::imwrite("depthmap.png", seg.distImage);
	cv::Canny( seg.dotImage, contour, cannyTh1, cannyTh2, 3);
    contour = cv::Scalar::all(255) - contour;

	cv::distanceTransform(contour, dist, CV_DIST_L2, 3);

    dist.convertTo(dist0, CV_8U, 1, 0);
	
	seg.distImage = dist0.clone();
	targetP = PCDContourFromRGBD(depth,foreground, seg.distImage,seg.dotImage);
		
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

	}

//std::cout << " target contour " << ntargetcontours << std::endl;
VecCoord targetContourpos;
targetContourpos.resize(ntargetcontours);
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
	}

	const VecCoord&  p0 = targetpos;
	targetPositions.setValue(p0);
    const VecCoord&  p1 = targetContourpos;
	targetContourPositions.setValue(p1);
	}

}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::KLTPointsTo3D()
{
	
  /*  int outside = 0;
	
	float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
	float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy
		
	VecCoord targetpos;
	targetpos.resize(tracker.getMaxFeatures());
	
	double znear;// = currentCamera->getZNear();
	double zfar;// = currentCamera->getZFar();
	
	 znear = 0.0716081;
	 zfar  = 72.8184;
            Vector3 pos;
            Vector3 col;
float xp, yp;
int id;

	 for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
                Mat3x3d m,mt;
				double xt, yt;
			     Vector3 xim0,xim1,xim2;
				 const int kk = k;
	tracker.getFeature(kk, id, xp, yp);
				int n0 = (int)yp;
				 int m0 = (int)xp;
				 float depthValue;
							if (!useRealData.getValue())
				 			depthValue = (float)depth.at<float>(2*yp,2*xp);
							else depthValue = (float)depth.at<float>(yp,xp);

			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)color.at<Vec4b>(yp,xp)[3];
			if ( depthValue>0 && depthValue < 1)                // if depthValue is not NaN
			{				
				double clip_z = (depthValue - 0.5) * 2.0;
                if (!useRealData.getValue()) pos[2] = -2*znear*zfar/(clip_z*(zfar-znear)-(zfar+znear));
				else pos[2] = depthValue;
				pos[0] = (xp - rgbIntrinsicMatrix(0,2)) * pos[2] * rgbFocalInvertedX;
				pos[1] = (yp - rgbIntrinsicMatrix(1,2)) * pos[2] * rgbFocalInvertedY;
				targetpos[id]=pos;
				
			}
		
            }
    const VecCoord&  p = targetpos;
	targetKLTPositions.setValue(p);*/

}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{
        if (dynamic_cast<simulation::AnimateBeginEvent*>(event))
	{
	
	int t = (int)this->getContext()->getTime();
	
	double timeT = (double)getTickCount();

        double timeAcq0 = (double)getTickCount();
	
	typename sofa::core::objectmodel::DataIO<DataTypes>::SPtr dataio;// = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();
	sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
	root->get(dataio);
	
	if (useRealData.getValue())
	{
	if (useSensor.getValue()){
		typename sofa::core::objectmodel::ImageConverter<DataTypes,DepthTypes>::SPtr imconv;// = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();

		root->get(imconv);
		color = imconv->color;
                //cv::imwrite("color02.png", color);
		depth = imconv->depth;
        	
		cv::Mat* imgl = new cv::Mat;
		*imgl = color.clone();
		cv::Mat* depthl = new cv::Mat;
		*depthl = depth.clone();
	
		dataio->listimg.push_back(imgl);
		dataio->listdepth.push_back(depthl);              
		
		depth00 = depth.clone();
                //cv::imwrite("depth02.png", depth00);
	}
	 else {
		color = dataio->color;
		depth = dataio->depth;
		depth00 = depth.clone();
		color_1 = dataio->color_1;
		color_5 = dataio->color_5.clone();
		color_4 = dataio->color_4.clone();
		color_3 = dataio->color_3.clone();
		color_2 = dataio->color_2.clone();
	 }
	}
	else dataio->readData();

        double timeAcq1 = (double)getTickCount();
        cout <<"time acq 0 " << (timeAcq1 - timeAcq0)/getTickFrequency() << endl;

        cv::namedWindow("image_sensor");
        cv::imshow("image_sensor",color);

        cv::namedWindow("depth_sensor");
        cv::imshow("depth_sensor",depth);

        cv::waitKey(1);

	if (t == 0)
	{

	if (useRealData.getValue())
	{
		initSegmentation();
		extractTargetPCD();
	}
	
	}
	else
        {
                if (t > 0 && t%niterations.getValue() == 0){
		
                if(useRealData.getValue())
		{	
                segment() ;

		if(!useContour.getValue())
		extractTargetPCD();
		else extractTargetPCDContour();
		
	    }
		else{
		segmentSynth();
		if(useContour.getValue())
		{
		cv::Mat distimg = seg.distImage;
		cv::Mat dotimg = seg.dotImage;	
		ContourFromRGBSynth(foreground, distimg,dotimg);
                }
                }

                cv::Mat* imgseg = new cv::Mat;
                *imgseg = foreground.clone();
                dataio->listimgseg.push_back(imgseg);
	
	}
	
    }


}
}

/*template <class DataTypes>
void RGBDDataProcessing<DataTypes>::CCDPointsTo3D(std::vector<pointCCD> pointsCCD)
{
	
    int outside = 0;
	
	float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
	float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy
		
	VecCoord targetpos;
	targetpos.resize(pointsCCD.size());
	
	double znear = currentCamera->getZNear();
	double zfar = currentCamera->getZFar();

            Vector3 pos;
            Vector3 col;
float xp, yp;
int id;

	 for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
                Mat3x3d m,mt;
				double xt, yt;
			     Vector3 xim0,xim1,xim2;
				 const int kk = k;
	            xp = pointsCCD[k].xu;
				yp = pointsCCD[k].xv;
				int n0 = (int)yp;
				 int m0 = (int)xp;
				 float depthValue;
							if (!useRealData.getValue())
				 			depthValue = (float)depth.at<float>(2*yp,2*xp);
							else depthValue = (float)depth.at<float>(yp,xp);

			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)color.at<Vec4b>(yp,xp)[3];
			if ( depthValue>0 && depthValue < 1)                // if depthValue is not NaN
			{				
				double clip_z = (depthValue - 0.5) * 2.0;
                if (!useRealData.getValue()) pos[2] = -2*znear*zfar/(clip_z*(zfar-znear)-(zfar+znear));
				else pos[2] = depthValue;
				pos[0] = (xp - rgbIntrinsicMatrix(0,2)) * pos[2] * rgbFocalInvertedX;
				pos[1] = (yp - rgbIntrinsicMatrix(1,2)) * pos[2] * rgbFocalInvertedY;
				targetpos[id]=pos;
				
			}
		
            }
    const VecCoord&  p = targetpos;
	targetCCDPositions.setValue(p);
}*/

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

}



}
}
} // namespace sofa

