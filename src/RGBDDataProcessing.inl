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

#pragma once
#include <limits>
#include <iterator>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/io.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/principal_curvatures.h>
#include <pcl/gpu/features/features.hpp>
#include <pcl/search/impl/search.hpp>


//#include "DataSource.hpp"

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

#include <algorithm>

#include "RGBDDataProcessing.h"

using std::cerr;
using std::endl;

namespace sofa {

namespace rgbdtracking {

//using namespace sofa::defaulttype;
//using namespace helper;

template <class DataTypes>
RGBDDataProcessing<DataTypes>::RGBDDataProcessing( )
 : Inherit()
    , l_imconv(initLink("imageconverter", "Link to image converter"))
    , l_dataio(initLink("dataio", "Link to dataio component"))
    , cameraIntrinsicParameters(initData(&cameraIntrinsicParameters,Vector4(),"cameraIntrinsicParameters","camera parameters"))

    , useContour(initData(&useContour,false,"useContour","Emphasize forces close to the target contours"))
    , useSensor(initData(&useSensor,true,"useSensor","Use real data"))
    , sensorType(initData(&sensorType, 0,"sensorType","Type of the sensor"))

    , niterations(initData(&niterations,1,"niterations","Number of iterations in the tracking process"))

    , samplePCD(initData(&samplePCD,4,"samplePCD","Sample step for the point cloud"))
    , sigmaWeight(initData(&sigmaWeight,(Real)4,"sigmaWeight","sigma weights"))
    , borderThdPCD(initData(&borderThdPCD,4,"borderThdPCD","border threshold on the target silhouette"))

    //segmentation params
    , segNghb(initData(&segNghb,8,"segnghb","Neighbourhood for segmentation"))
    , segImpl(initData(&segImpl,1,"segimpl","Implementation mode for segmentation (CUDA, OpenCV)"))
    , segMsk(initData(&segMsk,1,"segmsk","Mask type for segmentation"))

    , scaleSegmentation(initData(&scaleSegmentation,1,"downscalesegmentation","Down scaling factor on the RGB image for segmentation"))

    //display params
    , scaleImages(initData(&scaleImages,1,"downscaleimages","Down scaling factor on the RGB and depth images"))
    , displayImages(initData(&displayImages,true,"displayimages","Option to display RGB and Depth images"))
    , displayDownScale(initData(&displayDownScale,1,"downscaledisplay","Down scaling factor for the RGB and Depth images to be displayed"))
    , displaySegmentation(initData(&displaySegmentation,true,"displaySegmentation","Option to display the segmented image"))
    , drawPointCloud(initData(&drawPointCloud,false,"drawPointCloud"," "))
    , displayBackgroundImage(initData(&displayBackgroundImage,false,"displayBackgroundImage"," "))

    , saveImages(initData(&saveImages,false,"saveimages","Option to save RGB and Depth images on disk"))
    , useCurvature(initData(&useCurvature,false,"useCurvature"," "))

    //output
    , targetPositions(initData(&targetPositions,"targetPositions","Points of the target point cloud."))
    , targetNormals(initData(&targetNormals,"targetNormals","normals of the target point cloud."))
    , targetContourPositions(initData(&targetContourPositions,"targetContourPositions","Contour points of the target point cloud."))
    , targetBorder(initData(&targetBorder,"targetBorder","Border of the target point cloud."))
    , targetWeights(initData(&targetWeights,"targetWeights","weigths of the target point cloud."))
    , curvatures(initData(&curvatures, "curvatures", "output curvatures"))

    , cameraPosition(initData(&cameraPosition,"cameraPosition","Position of the camera w.r.t the point cloud"))
    , cameraOrientation(initData(&cameraOrientation,"cameraOrientation","Orientation of the camera w.r.t the point cloud"))
    , cameraChanged(initData(&cameraChanged,false,"cameraChanged","If the camera has changed or not"))

    , stopatinit(initData(&stopatinit,false,"stopatinit","stopatinit."))
    , safeModeSeg(initData(&safeModeSeg,false,"safeModeSeg","safe mode when segmentation fails"))
    , segTolerance(initData(&segTolerance,0.5,"segTolerance","tolerance or segmentation"))
{
    this->f_listening.setValue(true);
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

    if(displayImages.getValue()) {
        cv::namedWindow("image_sensor");
        cv::namedWindow("depth_sensor");
    }
    cv::namedWindow("image_segmented");

    Vector4 camParam = cameraIntrinsicParameters.getValue();

    rgbIntrinsicMatrix(0,0) = camParam[0];
    rgbIntrinsicMatrix(1,1) = camParam[1];
    rgbIntrinsicMatrix(0,2) = camParam[2];
    rgbIntrinsicMatrix(1,2) = camParam[3];

}

// Mouse Handler parameters
bool MouseEventHandler::drawing_box = false ;
bool MouseEventHandler::destroy=false ;
cv::Rect MouseEventHandler::box = cvRect(0,0,1,1);
template <class DataTypes>
void RGBDDataProcessing<DataTypes>::initSegmentation()
{
    std::cout << "initseg" << std::endl ;
    cv::Mat downsampledbox,downsampled;

    int scaleSeg = scaleSegmentation.getValue();
    if (scaleSeg>1) {
        cv::resize(color, downsampledbox, cv::Size(color.cols/scaleSeg, color.rows/scaleSeg));
    } else {
        downsampledbox = color.clone();
    }
	
    //cv::imwrite("colorinit.png", color);
	

    const char* name = "image";

    cv::namedWindow(name);

    //cv::imshow("image",	color);
    //tempm.resize(image.step1());

    // Set up the callback

    cv::Mat temp = downsampledbox.clone();
    cv::setMouseCallback(name, MouseEventHandler::my_mouse_callback, (void*) &temp);

    // Main loop
    cv::Mat temp1 ;
    while(true) {
        if (MouseEventHandler::destroy) {
            cv::destroyWindow(name);
            break;
        }
        temp1 = temp.clone();

        if (MouseEventHandler::drawing_box) {
            MouseEventHandler::draw_box(temp1, MouseEventHandler::box);
        }

        cv::moveWindow(name, 200, 100);
        cv::imshow(name, temp1);

        int key = waitKey(10);
        if ((char)key == 27) {
            break;
        }
    }

    // delete temp1;
    temp1.release();
    cv::setMouseCallback(name, NULL, NULL);

    //cvDestroyWindow(name);
    cv::Rect rectangle(69,47,198,171);

    rectangle = MouseEventHandler::box;
    seg.setRectangle(rectangle);
    
    if (scaleSeg>1) {
        cv::resize(
            color,
            downsampled,
            cv::Size(color.cols/scaleSeg, color.rows/scaleSeg));
    } else {
        downsampled = color.clone();
    }

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

    if (waitKey(20) == 27) {
    //wait for 'esc' key press for 30ms. If 'esc' key is pressed, break loop
        std::cout << "esc key is pressed by user" << std::endl;
    }
   
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::segment()
{
    helper::AdvancedTimer::stepBegin("Segmentation") ;

    int scaleSeg = scaleSegmentation.getValue();
    cv::Mat downsampled ;
    if (scaleSeg>1) {
        cv::resize(color, downsampled, cv::Size(color.cols/scaleSeg, color.rows/scaleSeg));
    } else {
        downsampled = color.clone();
    }

    seg.updateMask(foregroundS);
    //cv::GaussianBlur( downsampled, downsampled1, cv::Size( 3, 3), 0, 0 );
    //cv::imwrite("downsampled.png", downsampled);
    seg.updateSegmentation(downsampled,foregroundS);
    cv::resize(foregroundS, foreground, color.size(), INTER_NEAREST);

    if(useContour.getValue()) {
        cv::resize(seg.dotImage, dotimage, color.size(), INTER_NEAREST);
        cv::resize(seg.distImage, distimage, color.size(), INTER_NEAREST);
    }
    //foreground = foregroundS.clone();
    helper::AdvancedTimer::stepEnd("Segmentation") ;

    //seg.updateSegmentationCrop(downsampled,foreground);

    // display result
    if (displaySegmentation.getValue()){
        cv::imshow("image_segmented",foregroundS);
        cv::waitKey(1);
    }
}

template <class DataTypes>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr RGBDDataProcessing<DataTypes>::PCDFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage)
{
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr outputPointcloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    //pcl::PointCloud<pcl::PointXYZRGB> pointcloud;
    outputPointcloud->points.resize(0);
		
    int sample;
	
    switch (sensorType.getValue()) {
        // FLOAT ONE CHANNEL
        case 0:
            sample = 2;
            break;
        case 1:
            sample = samplePCD.getValue();
            break;
        default :
            break ;
    }

    float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
    float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy
    pcl::PointXYZRGB newPoint;
    std::cout << "inrgbpcl" <<std::endl ;
    for (int i=0;i<(int)depthImage.rows/sample;i++) {
        for (int j=0;j<(int)depthImage.cols/sample;j++) {

            float depthValue = (float)depthImage.at<float>(sample*i,sample*j);//*0.819;
            int avalue = (int)rgbImage.at<Vec4b>(sample*i,sample*j)[3];

            if (avalue > 0 && depthValue>0) {
                std::cout << "IN" << std::endl ;
                // if depthValue is not NaN
                // Find 3D position respect to rgb frame:
                newPoint.z = depthValue;
                newPoint.x = (sample*j - rgbIntrinsicMatrix(0,2)) * newPoint.z * rgbFocalInvertedX;
                newPoint.y = (sample*i - rgbIntrinsicMatrix(1,2)) * newPoint.z * rgbFocalInvertedY;
                newPoint.r = rgbImage.at<cv::Vec4b>(sample*i,sample*j)[2];
                newPoint.g = rgbImage.at<cv::Vec4b>(sample*i,sample*j)[1];
                newPoint.b = rgbImage.at<cv::Vec4b>(sample*i,sample*j)[0];
                outputPointcloud->points.push_back(newPoint);
                std::cout << "OUT" << std::endl ;
            }
        }
    }
    std::cout << "outrgbpcl" << outputPointcloud->size() <<std::endl ;

    if (useCurvature.getValue()){
        int sample1 = samplePCD.getValue();
        pcl::PointCloud<pcl::PointXYZ>::Ptr outputPointcloud1(new pcl::PointCloud<pcl::PointXYZ>);
        outputPointcloud1->points.resize(0);

        pcl::PointXYZ newPoint1;

        for (int i=0;i<(int)depthImage.rows/sample1;i++) {
            for (int j=0;j<(int)depthImage.cols/sample1;j++) {
                float depthValue = (float)depthImage.at<float>(sample1*i,sample1*j);//*0.819;
                int avalue = (int)rgbImage.at<Vec4b>(sample1*i,sample1*j)[3];
                if (avalue > 0 && depthValue>0) {
                    // if depthValue is not NaN
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
        std::vector<double> curvs;
        for (int k = 0; k < principalCurvatures->points.size (); k++) {
            pcl::PrincipalCurvatures descriptor0 = principalCurvatures->points[k];
            double curv = abs(descriptor0.pc1*descriptor0.pc1);
            curvs.push_back(curv);
        }

        curvatures.setValue(curvs);
    }
    return outputPointcloud;
}

template <class DataTypes>
pcl::PointCloud<pcl::PointXYZRGB>::Ptr RGBDDataProcessing<DataTypes>::PCDContourFromRGBD(cv::Mat& depthImage, cv::Mat& rgbImage, cv::Mat& distImage, cv::Mat& dotImage) {
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

    frgd = rgbImage;
    sample = samplePCD.getValue();

    float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);// 1/fx
    float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);// 1/fy
    ntargetcontours = 0;
    int jj = 0;
    double totalweights = 0;
    pcl::PointXYZRGB newPoint;

    for (int i=0;i<(int)depthImage.rows/sample;i++) {
        for (int j=0;j<(int)depthImage.cols/sample;j++) {
            float depthValue = (float)depthImage.at<float>(sample*i,sample*j);//*0.819;
            int avalue = (int)frgd.at<Vec4b>(sample*i,sample*j)[3];
            int bvalue = (int)distimg.at<uchar>(sample*i,sample*(j));
            int dvalue = (int)dotimg.at<uchar>(sample*i,sample*(j));

            if (dvalue == 0 && depthValue>0) {
                // if depthValue is not NaN
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

                if (avalue > 0 && bvalue < borderThdPCD.getValue()) {
                    targetborder.push_back(true);
                    ntargetcontours++;
                } else {
                    targetborder.push_back(false);
                }
            }
        }
    }
	
    for (int i=0; i < targetweights.size();i++) {
        targetweights[i]*=((double)targetweights.size()/totalweights);
        //std::cout << " weights " << totalweights << " " << (double)targetweights[i] << std::endl;
    }

    targetWeights.setValue(targetweights);
    targetBorder.setValue(targetborder);

    if (useCurvature.getValue()) {
        int sample1 = samplePCD.getValue();
        pcl::PointCloud<pcl::PointXYZ>::Ptr outputPointcloud1(new pcl::PointCloud<pcl::PointXYZ>);
        outputPointcloud1->points.resize(0);

        pcl::PointXYZ newPoint1;

        for (int i=0;i<(int)depthImage.rows/sample1;i++) {
            for (int j=0;j<(int)depthImage.cols/sample1;j++) {
                float depthValue = (float)depthImage.at<float>(sample1*i,sample1*j);//*0.819;
                int avalue = (int)rgbImage.at<Vec4b>(sample1*i,sample1*j)[3];
                if (avalue > 0 && depthValue>0) {
                    // if depthValue is not NaN
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

        for (int k = 0; k < principalCurvatures->points.size (); k++) {
            pcl::PrincipalCurvatures descriptor0 = principalCurvatures->points[k];
            double curv = abs(descriptor0.pc1*descriptor0.pc1);
            curvs.push_back(curv);
        }

        curvatures.setValue(curvs);
        std::cout << " curvature " << descriptor << std::endl;
    }

    return outputPointcloud;
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::extractTargetPCD() {

    int t = (int)this->getContext()->getTime();
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetP ;

    targetP.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
    targetP = PCDFromRGBD(depth,foreground);
    std::cout << "111" << std::endl ;
    VecCoord targetpos;
	
    if (targetP->size() <= 10) {
        return ;
    }

    target.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
    target = targetP;
    targetpos.resize(target->size());

    Vector3 pos;
    for (unsigned int i=0; i<target->size(); i++) {
        pos[0] = (double)target->points[i].x;
        pos[1] = (double)target->points[i].y;
        pos[2] = (double)target->points[i].z;
        targetpos[i]=pos;
        //std::cout << " target " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    }
    const VecCoord&  p = targetpos;

    if (safeModeSeg.getValue()) {
        bool guard =
            abs((double)p.size() - (double)sizeinit)/(double)sizeinit
            < segTolerance.getValue() ;
        if (t<20*niterations.getValue()) {
            sizeinit = p.size();
            targetPositions.setValue(p);
        } else if (guard) {
            targetPositions.setValue(p);
        }
    } else {
        targetPositions.setValue(p);
    }
	
}


template <class DataTypes>
void RGBDDataProcessing<DataTypes>::extractTargetPCDContour() {
#define CANNY_TH1 150
#define CANNY_TH2 80
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr targetP ;

    targetP.reset(new pcl::PointCloud<pcl::PointXYZRGB>);


    cv::Mat contour, dist, dist0;
	
    //cv::imwrite("depthmap.png", seg.distImage);
    cv::Canny(dotimage, contour, CANNY_TH1, CANNY_TH2, 3);
    contour = cv::Scalar::all(255) - contour;
    cv::distanceTransform(contour, dist, CV_DIST_L2, 3);
    dist.convertTo(dist0, CV_8U, 1, 0);
    seg.distImage = dist0.clone();
    targetP = PCDContourFromRGBD(depth,foreground, distimage, dotimage);
    VecCoord targetpos;
	
    if (targetP->size() <= 10) {
        return ;
    }
    target.reset(new pcl::PointCloud<pcl::PointXYZRGB>);
    target = targetP;
    targetpos.resize(target->size());

    for (unsigned int i=0; i<target->size(); i++) {
        Vector3 pos (
            (double)target->points[i].x,
            (double)target->points[i].y,
            (double)target->points[i].z
        ) ;
        targetpos[i]=pos;
    }

    VecCoord targetContourpos;
    targetContourpos.resize(ntargetcontours);
    int kk = 0;
    for (unsigned int i=0; i<target->size(); i++) {
        if (targetBorder.getValue()[i]) {
            Vector3 pos (
                (double)target->points[i].x,
                (double)target->points[i].y,
                (double)target->points[i].z
            );
            targetContourpos[kk++]=pos;
        }
    }

    const VecCoord&  p0 = targetpos;
    targetPositions.setValue(p0);
    const VecCoord&  p1 = targetContourpos;
    targetContourPositions.setValue(p1);
//    std::cout << " target contour " << p1.size() << std::endl;

}

template<class DataTypes>
void RGBDDataProcessing<DataTypes>::setCameraPose()
{
    pcl::PointCloud<pcl::PointXYZRGB>& point_cloud = *target;
    if (point_cloud.size() > 0) {
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
bool RGBDDataProcessing<DataTypes>::loadFromDataIO () {
    color = l_dataio->color;
    depth = l_dataio->depth;
    return l_dataio->newImages.getValue();
}

template <class DataTypes>
bool RGBDDataProcessing<DataTypes>::loadFromImageConverter () {
    if (!l_imconv->depth.empty() && !l_imconv->color.empty()) {
        // this is for scaling images
        int scale = (scaleImages.getValue() < 1) ? 1 : scaleImages.getValue() ;
        cv::resize(
            l_imconv->depth, depth,
            cv::Size(
                l_imconv->depth.cols/scale,
                l_imconv->depth.rows/scale),
            0, 0);
        cv::resize(
            l_imconv->color, color,
            cv::Size(
                l_imconv->color.cols/scale,
                l_imconv->color.rows/scale),
            0, 0);
        return l_imconv->newImages.getValue();
    } else {
        return false ;
    }
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::saveToDataIO () {
    cv::Mat* imgl = new cv::Mat;
    *imgl = color.clone();
    cv::Mat* depthl = new cv::Mat;
    *depthl = depth.clone();

    l_dataio->listimg.push_back(imgl);
    l_dataio->listdepth.push_back(depthl);
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::displayDownScaledImage () {
    int scale = displayDownScale.getValue();
    cv::Mat colorS, depthS;
    cv::resize(depth, depthS, cv::Size(depth.cols/scale, depth.rows/scale), 0, 0);
    cv::resize(color, colorS, cv::Size(color.cols/scale, color.rows/scale), 0, 0);

    cv::imshow("image_sensor",colorS);
    cv::waitKey(1);
    /*cv::imshow("depth_sensor",depthS);
    cv::waitKey(1);*/
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::saveSegmentationToDataIO () {
    cv::Mat *imgseg = new cv::Mat;
    *imgseg = foreground.clone();
    l_dataio->listimgseg.push_back(imgseg);
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event) {
    if (dynamic_cast<simulation::AnimateBeginEvent*>(event)) {
        if (! l_dataio && ! l_imconv) {
            std::cerr << "(RGBDDataProcessing) Can't find link to dataio and imageconverter" << std::endl ;
            return ;
        }
        helper::AdvancedTimer::stepBegin("RGBDDataProcessing") ;
        helper::AdvancedTimer::stepBegin("RGBDDataProcessingImageLoading") ;
        int t = (int)this->getContext()->getTime();

        bool newimages = false;

        if (l_imconv) {
        //load from image converter
            newimages = loadFromImageConverter() ;
        } else if (l_dataio) {
        // load from dataio
            newimages = loadFromDataIO() ;
        }
        /// end
        helper::AdvancedTimer::stepEnd("RGBDDataProcessingImageLoading") ;
        //std::cout << "newimages " << newimages << std::endl;

        /// dispay image with cv::imshow
        if (displayImages.getValue() &&
            displayDownScale.getValue() > 0 &&
            !depth.empty() &&
            !color.empty()
        ) {
            displayDownScaledImage();
        }

        /// save images if desired
        if (saveImages.getValue() && t % niterations.getValue() == 0 ) {
            saveToDataIO();
        }

        //// if image loading successful
        static bool initsegmentation = true ;
        if (newimages) {
            if (initsegmentation) {
            // 1st segmentation
                initSegmentation();
                extractTargetPCD();
                setCameraPose();
                initsegmentation = false;
                std::cout << "done" << std::endl ;
            } else if(!stopatinit.getValue()) {
            // nth iteration
                segment() ;
                helper::AdvancedTimer::stepBegin("PCDExtraction") ;
                if(!useContour.getValue()) {
                    extractTargetPCD();
                } else {
                    extractTargetPCDContour();
                }
                helper::AdvancedTimer::stepEnd("PCDExtraction") ;
            }

            if (l_dataio &&
                saveImages.getValue() &&
                t%niterations.getValue() == 0
            ) {
                saveSegmentationToDataIO();
            }
            cameraChanged.setValue(false);
        }
        helper::AdvancedTimer::stepEnd("RGBDDataProcessing") ;
    }
}

template <class DataTypes>
void RGBDDataProcessing<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

    helper::ReadAccessor< Data< VecCoord > > xtarget(targetPositions);
    vparams->drawTool()->saveLastState();

    if (displayBackgroundImage.getValue()) {
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

        for (unsigned int i=0; i< xtarget.size(); i++) {
            points.resize(0);
            point = DataTypes::getCPos(xtarget[i]);
            points.push_back(point);
            // std::cout << curvatures.getValue()[i] << std::endl;
            //if (targetWeights.getValue().size()>0) vparams->drawTool()->drawPoints(points, 10, sofa::defaulttype::Vec<4,float>(0.5*targetWeights.getValue()[i],0,0,1));
            vparams->drawTool()->drawPoints(points, 10, sofa::defaulttype::Vec<4,float>(1,0.5,0.5,1));
        }

    }

}

} // rgbdtracking

} // namespace sofa

