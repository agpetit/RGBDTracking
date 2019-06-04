/*
 * segmentation.h
 *
 *  Created on: May 2, 2014
 *      Author: apetit
 */

#pragma once
/*#ifndef min
#define min(a,b) ((a < b) ? a:b)
#endif
#ifndef max
#define max(a,b) ((a > b) ? a:b)
#endif*/

#include <opencv/cv.h>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <RGBDTracking/config.h>
#include <RGBDTracking/src/Segmentor/CVSegmentation.h>
#include <RGBDTracking/src/Segmentor/CUDASegmentation.h>

#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>

#include <GL/glew.h>
#include <GL/freeglut.h>

using namespace cv;
using namespace std;

#define MAX_IMAGE_SIZE 1024.0f

class segmentation {

public :
    typedef enum {
        CVGRAPHCUT,
        CUDAGRAPHCUT,
        APGRAPHCUT
    } implementation;

    CVSegmentation cvseg ;
    CUDASegmentation cudaseg ;

    int crop_width, crop_height;

    int neighborhood;

    //struct cudaGraphicsResource *pbo_resource;

    //segmentationParameters
    implementation type;
    masktype mskt;

    cv::Rect rectangle;
    cv::Mat mask ;
    cv::Mat distImage, dotImage;

    segmentation();
    virtual ~segmentation();

    void init(int nghb, int impl, int msk);
    void setRectangle(cv::Rect _rectangle){rectangle = _rectangle;}
    void segmentationFromRect(cv::Mat &image, cv::Mat &foreground);
    void updateMask(cv::Mat &foreground);
    void updateSegmentation(cv::Mat &image, cv::Mat &foreground);
    void filter(cv::Mat &out,cv::Mat &dt,cv::Mat &dot);
    void updateSegmentationCrop(cv::Mat &image, cv::Mat &foreground);
    void clear();

    //void setSegmentationParameters(segmentationParameters &_segParam){segParam = _segParam;}
    #ifdef HAVECUDA // && (CUDART_VERSION == 7000)
    FIBITMAP* convertCVFree(cv::Mat &in);
    void convertFreeCV(FIBITMAP* h_Image,cv::Mat &out);
    void convertFreeCV8(FIBITMAP* h_Image,cv::Mat &out);
    #endif
    void maskFromDt(cv::Mat &_dt, cv::Mat &mask_);
    void trimapFromDt(cv::Mat &_dt,cv::Mat &dot);
    inline int cudaDeviceInit();

protected:
    void updateMaskBBox(cv::Mat &foreground);
    void updateMaskContour(cv::Mat &foreground);
};
