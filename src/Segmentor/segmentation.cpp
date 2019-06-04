/*
 * segmentation.cpp
 *
 *  Created on: May 2, 2014
 *      Author: apetit
 */

#include "segmentation.h"

#include <sofa/helper/AdvancedTimer.h>

segmentation::segmentation() {
    // TODO Auto-generated constructor stub
    neighborhood = 8;
    mskt = BBOX;
    type = CVGRAPHCUT;

    cvseg = CVSegmentation () ;
    cudaseg = CUDASegmentation () ;
}

segmentation::~segmentation() {
}

void segmentation::init(int nghb, int impl, int msk)
{
    neighborhood = nghb;
    switch (impl) {
    case 0:
        type = CVGRAPHCUT;
        break;
    case 1:
        type = CUDAGRAPHCUT;
        break;
    }
    switch (msk) {
    case 0:
        mskt = BBOX;
        break;
    case 1:
        mskt = CONTOUR;
        break;
    }

    if (type == CVGRAPHCUT) {
        mskt = BBOX;
    }
}

void segmentation::segmentationFromRect(cv::Mat &image, cv::Mat &foreground)
{

    switch(type){
        case CVGRAPHCUT:
            cvseg.segmentationFromRect(image, foreground);
            break;
        case CUDAGRAPHCUT:
            cudaseg.segmentationFromRect(image, foreground);
            break;
    }
}

void segmentation::updateMaskBBox(cv::Mat &foreground) {

    std::vector<cv::Point> ptfgd;
    ptfgd.resize(0);

    switch(type){
        case CVGRAPHCUT:
            /// instantiante class please
            cvseg.updateMaskBBox(ptfgd);
            mask = cvseg.mask ; // update mask
            break;
        case CUDAGRAPHCUT:
            cudaseg.updateMaskBBox(foreground, ptfgd);
            break;
    }

    rectangle = cv::boundingRect(ptfgd);
    //std::cout << " ptfgd " << frame_count << " rect1 " << rectangle.x << " " << rectangle.y << " rect2 " << rectangle.width << " " << rectangle.height << std::endl;

    rectangle.x -= 10;
    rectangle.y -= 10;
    rectangle.height += 20;
    rectangle.width += 20;

//    std::cout << " rect1 " << rectangle.x << " " << rectangle.y << std::endl;

    for(int x = 0; x<mask.cols; x++) {
        for(int y = 0; y<mask.rows; y++) {
            if (x < rectangle.x ||
                x > rectangle.x + rectangle.width ||
                y < rectangle.y ||
                y > rectangle.y + rectangle.height
            ) {
                    mask.at<uchar>(y,x) = 0;
            }
        }
    }
}

void segmentation::updateMaskContour(cv::Mat &foreground)
{
    cv::Mat distanceMap(foreground.size(),CV_8U,cv::Scalar(255));
    cv::Mat dot(foreground.size(),CV_8U,cv::Scalar(0));
    filter(foreground,distanceMap,dot);
    //maskFromDt(distanceMap,mask1);
    trimapFromDt(distanceMap,dot);

    //std::cout << " rect1 " << rectangle.x << " " << rectangle.y << std::endl;

    distImage = distanceMap;
    dotImage = dot;
}

void segmentation::updateMask(cv::Mat &foreground) {
    switch(mskt){
        case BBOX:
            updateMaskBBox(foreground);
            break ;
        case CONTOUR:
            updateMaskContour(foreground);
            break;
    }
}

void segmentation::updateSegmentation(cv::Mat &image,cv::Mat &foreground) {

    switch(type){
        case CVGRAPHCUT:
            cvseg.updateSegmentation(foreground, image);
            break;

        case CUDAGRAPHCUT:
            cudaseg.updateSegmentation(foreground, image);
            break;
    }

}

void segmentation::updateSegmentationCrop(cv::Mat &image,cv::Mat &foreground)
{
    switch(type){
        case CVGRAPHCUT:
            cvseg.updateSegmentationCrop(image, foreground);
            break;

        case CUDAGRAPHCUT:
            cudaseg.updateSegmentationCrop(foreground, image);
            break;
    }
}

void segmentation::filter(cv::Mat &out,cv::Mat &dt,cv::Mat &dot) {
    cv::Mat dt0;

    //cv::imwrite("out.png",out);
    for(int x = 0; x<out.cols; x++) {
        for(int y = 0; y<out.rows; y++) {
            if (out.at<cv::Vec4b>(y,x)[0] > 0 ||
                out.at<cv::Vec4b>(y,x)[1] > 0 ||
                out.at<cv::Vec4b>(y,x)[2] > 0
            ) {
                dot.at<uchar>(y,x) = 0;
            } else {
                dot.at<uchar>(y,x) = 255;
            }
        }
    }

    cv::Mat edges;

    /*cv::namedWindow("Image");
    cv::imshow("Image",out);*/

    cv::Canny(dot,edges,10,350,3);
    edges = cv::Scalar::all(255) - edges;

    //cv::imwrite("edges.png",edges);

    cv::distanceTransform(edges, dt0, CV_DIST_L2, 3);

    //cv::imwrite("dt0.png",dt0);
    //dt0 *= 5000;
    //pow(dt0, 0.5, dt0);
    dt0.convertTo(dt, CV_8U, 1, 0);

    distImage = dt;
    dotImage = dot;

    //cv::imwrite("dt.png",dt);

    /*cv::namedWindow("Image");
    cv::imshow("Image",dot);*/

    /*cv::namedWindow("Canny");
    cv::imshow("Canny",edges);*/

    //normalize(dt, dt, 0.0, 1.0, NORM_MINMAX);
    //imshow("normalized", dt);

    //waitKey(0);
    //getchar();
}

void segmentation::maskFromDt(cv::Mat &_dt,cv::Mat &mask_)
{
    cv::Mat mask(_dt.size(),CV_8U,cv::Scalar(0));

    for(int x = 0; x<_dt.cols; x++) {
        for(int y = 0; y<_dt.rows; y++) {
            if (_dt.at<uchar>(y,x) > 10) {
                mask.at<uchar>(y,x) = 0;
            } else {
                mask.at<uchar>(y,x) = 255;
            }
        }
    }

    mask_ = mask;

    /*cv::namedWindow("Mask");
    cv::imshow("Mask",mask_);*/

    //waitKey(0);
    //getchar();
}

void segmentation::trimapFromDt(cv::Mat &_dt,cv::Mat &dot) {

    mask = cv::Mat::zeros(_dt.size(),CV_8U);
    std::vector<cv::Point> ptfgd;
    ptfgd.resize(0);
    cv::Mat strip = dot.clone();
    int band = 10;
    for(int x = 0; x<_dt.cols; x++) {
        for(int y = 0; y<_dt.rows; y++) {
            cv::Point pt;
            if (_dt.at<uchar>(y,x) > band && dot.at<uchar>(y,x) == 255) {
                mask.at<uchar>(y,x) = 0;
                strip.at<uchar>(y,x) = 0;
            } else if (_dt.at<uchar>(y,x) > band && dot.at<uchar>(y,x) == 0) {
                mask.at<uchar>(y,x) = 2;
                pt.x = x;
                pt.y = y;
                ptfgd.push_back(pt);
                strip.at<uchar>(y,x) = 255;
            } else {
                mask.at<uchar>(y,x) = 1;
                pt.x = x;
                pt.y = y;
                ptfgd.push_back(pt);
                strip.at<uchar>(y,x) = 127;
            }
        }
    }

    std::cout << " ok bd rect 0 " << ptfgd.size() <<  std::endl;
    if (ptfgd.size() >= 200) {
        rectangle = cv::boundingRect(ptfgd);
    }

    //cv::imwrite("strip.png",strip);
    //cv::imwrite("dot.png",dot);

    //std::cout << " ok bd rect 1 " << rectangle.x << " " << rectangle.y << std::endl;
    if (ptfgd.size() < 200) {
        std::cout << " fail " << std::endl;
        for(int x = 0; x<_dt.cols; x++) {
            for(int y = 0; y<_dt.rows; y++) {
                if (_dt.at<uchar>(y,x) > band && dot.at<uchar>(y,x) == 255) {
                    mask.at<uchar>(y,x) = 1;
                } else if (_dt.at<uchar>(y,x) > band && dot.at<uchar>(y,x) == 0) {
                    mask.at<uchar>(y,x) = 2;
                } else {
                    mask.at<uchar>(y,x) = 1;
                }
            }
        }
    } else {
        rectangle.x -= 0;
        rectangle.y -= 0;
        rectangle.height += 0;
        rectangle.width += 0;
    }

    //std::cout << " ptfgd " << frame_count << " rect1 " << rectangle.x << " " << rectangle.y << " rect2 " << rectangle.width << " " << rectangle.height << std::endl;
    //std::cout << " rect1 " << rectangle.x << " " << rectangle.y << std::endl;

    /*cv::namedWindow("Mask");
    cv::imshow("Mask",mask);*/
}
