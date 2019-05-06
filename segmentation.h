/*
 * segmentation.h
 *
 *  Created on: May 2, 2014
 *      Author: apetit
 */

#ifndef SEGMENTATION_H_
#define SEGMENTATION_H_

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

using namespace cv;
using namespace std;

#define MAX_IMAGE_SIZE 1024.0f

class segmentation {

public :
	typedef enum
	{
	CVGRAPHCUT,
	CUDAGRAPHCUT,
	APGRAPHCUT
	}implementation;

	typedef enum
	{
	BBOX,
	CONTOUR
	}masktype;

int width, height;
int crop_width, crop_height;

#ifdef HAVECUDA
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


int neighborhood;

//struct cudaGraphicsResource *pbo_resource;

implementation type;
masktype mskt;
//segmentationParameters segParam;

cv::Rect rectangle;
cv::Mat mask, maskimg;
cv::Mat distImage, dotImage;



segmentation();
virtual ~segmentation();

void init(int nghb, int impl, int msk);
void setRectangle(cv::Rect _rectangle){rectangle = _rectangle;}
void segmentationFromRect(cv::Mat &image, cv::Mat &foreground);
void clear();
//void setSegmentationParameters(segmentationParameters &_segParam){segParam = _segParam;}
void updateMask(cv::Mat &foreground);
void updateSegmentation(cv::Mat &image, cv::Mat &foreground);
void updateSegmentationCrop(cv::Mat &image, cv::Mat &foreground);
void saveResult(const char *filename);
void getResult(cv::Mat &out);
void getResultCrop(cv::Mat &out);
bool verifyResult(const char *filename);
void filter(cv::Mat &out,cv::Mat &dt,cv::Mat &dot);
#ifdef HAVECUDA // && (CUDART_VERSION == 7000)
FIBITMAP* convertCVFree(cv::Mat &in);
void convertFreeCV(FIBITMAP* h_Image,cv::Mat &out);
void convertFreeCV8(FIBITMAP* h_Image,cv::Mat &out);
#endif
void maskFromDt(cv::Mat &_dt, cv::Mat &mask_);
void trimapFromDt(cv::Mat &_dt,cv::Mat &dot);
inline int cudaDeviceInit();
void clean();
};

#endif /* SEGMENTATION_H_ */
