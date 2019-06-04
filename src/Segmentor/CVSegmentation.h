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

class CVSegmentation {
public :
    CVSegmentation () ;
    ~CVSegmentation() ;
    void segmentationFromRect (cv::Mat &image, cv::Mat &foreground) ;
    void updateMaskBBox(std::vector<cv::Point> & ptfgd) ;
    void updateSegmentation(cv::Mat &foreground, cv::Mat &image) ;
    void updateSegmentationCrop(cv::Mat &image, cv::Mat &foreground) ;

    cv::Mat mask, maskimg ;
    cv::Rect rectangle;
} ;
