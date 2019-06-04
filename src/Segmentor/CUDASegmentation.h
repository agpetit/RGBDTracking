/*
 * segmentation.h
 *
 *  Created on: May 2, 2014
 *      Author: apetit
 */

#pragma once

#include <RGBDTracking/config.h>


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

#ifdef HAVECUDA // && (CUDART_VERSION == 7000)
#include <cuda_runtime.h>
#include <npp.h>
#include <nppi.h>

#include <cuda_gl_interop.h>
#include <helper_cuda.h>
#include <helper_string.h>

#include "FreeImage.h"
#include "cudaSegmentation.h"
#endif

#include <opencv/cv.h>
#include <opencv2/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

using namespace cv;
using namespace std;

#define MAX_IMAGE_SIZE 1024.0f

typedef enum {
    BBOX,
    CONTOUR
} masktype;

class CUDASegmentation {
public :
    void segmentationFromRect(cv::Mat &image, cv::Mat &foreground) ;
    void updateMaskBBox(cv::Mat &foreground, std::vector<cv::Point> & ptfgd) ;
    void updateSegmentation(cv::Mat &foreground, cv::Mat &image) ;
    void updateSegmentationCrop(cv::Mat &foreground, cv::Mat &image) ;
    void clean () ;
    CUDASegmentation () ;
    ~CUDASegmentation () ;
protected :
    void saveResult(const char *filename) ;
    void getResult(cv::Mat &out) ;
    void getResultCrop(cv::Mat &out) ;
    bool verifyResult(const char *filename) ;

#ifdef HAVECUDA
    FIBITMAP* convertCVFree(cv::Mat &in) ;
#endif

    cv::Mat mask ;
    masktype mskt ;

#ifdef HAVECUDA
    int width, height;
    NppiRect rect;
    uchar4 *d_image;
    size_t image_pitch;

    unsigned char *d_trimap;
    size_t trimap_pitch;
    unsigned char *d_dt;

    NppiRect crop_rect;
    uchar4 *d_crop_image;
    size_t crop_image_pitch;

    unsigned char *d_crop_trimap;
    size_t crop_trimap_pitch;
    unsigned char *d_crop_dt;

    cudaSegmentation *cudaseg;
#endif

} ;
