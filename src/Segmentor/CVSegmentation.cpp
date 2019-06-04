/*
 * segmentation.cpp
 *
 *  Created on: May 2, 2014
 *      Author: apetit
 */

#include "CVSegmentation.h"

#include <sofa/helper/AdvancedTimer.h>

CVSegmentation::CVSegmentation () {
}
CVSegmentation::~CVSegmentation () {
}

void CVSegmentation::segmentationFromRect(cv::Mat &image, cv::Mat &foreground)
{
    sofa::helper::AdvancedTimer::stepBegin("CVGraphCut") ;
    cv::Mat bgModel,fgModel;
    cv::grabCut(
        image,                  // input image
        mask,                   // segmentation result
        rectangle,              // rectangle containing foreground
        bgModel,fgModel,        // models
        1,                      // number of iterations
        cv::GC_INIT_WITH_RECT); // use rectangle
    // Get the pixels marked as likely foreground
    sofa::helper::AdvancedTimer::stepEnd("CVGraphCut") ;

    cv::compare(mask,cv::GC_PR_FGD,maskimg,cv::CMP_EQ);
    // Generate output image
    image.copyTo(foreground,maskimg);
}

void CVSegmentation::updateMaskBBox(std::vector<cv::Point> & ptfgd)
{
    for(int x = 0; x<mask.cols; x++) {
        for(int y = 0; y<mask.rows; y++) {
            if (mask.at<uchar>(y,x) == 1 || mask.at<uchar>(y,x) == 3) {
                    Point pt ;
                    pt.x = x;
                    pt.y = y;
                    ptfgd.push_back(pt);
            } else if(mask.at<uchar>(y,x) == 0) {
                    mask.at<uchar>(y,x) = 2;
            }
        }
    }
}

void CVSegmentation::updateSegmentation(cv::Mat &foreground, cv::Mat &image)
{
    sofa::helper::AdvancedTimer::stepBegin("CVGraphCut") ;

    cv::Mat bgModel,fgModel;

    cv::Mat _image = image.clone();
    cv::Mat _mask = mask.clone();
    cv::Mat foreground_;

    cv::grabCut(_image,_mask,rectangle,bgModel,fgModel,5,cv::GC_INIT_WITH_MASK);

    sofa::helper::AdvancedTimer::stepEnd("CVGraphCut") ;

    //std::cout << " mask " << mask << std::endl;
    cv::compare(_mask,cv::GC_PR_FGD,maskimg,cv::CMP_EQ);
    // Generate output image
    //cv::Mat foreground(image.size(),CV_8UC3,cv::Scalar(255,255,255));
    _image.copyTo(foreground_,maskimg); // bg pixels not copied
    cv::Mat alpha(image.size(),CV_8UC1,Scalar(0));

    for(int x = 0; x<image.cols; x++) {
        for(int y = 0; y<image.rows; y++){
            if (foreground_.at<cv::Vec3b>(y,x)[0] == 0 &&
                foreground_.at<cv::Vec3b>(y,x)[1] == 0 &&
                foreground_.at<cv::Vec3b>(y,x)[2] == 0
            ) {
                  alpha.at<uchar>(y,x) = 0;
            } else {
                alpha.at<uchar>(y,x) = 255;
            }
        }
    }

    cv::Mat rgb[4];
    cv::split(foreground_,rgb);

    cv::Mat rgba[4]={rgb[0],rgb[1],rgb[2],alpha};
    cv::merge(rgba,4,foreground_);
    foreground = foreground_.clone();
}

void CVSegmentation::updateSegmentationCrop(cv::Mat &image, cv::Mat &foreground) {
    sofa::helper::AdvancedTimer::stepBegin("CVGraphCut") ;

    cv::Mat bgModel,fgModel;
    cv::grabCut(image,mask,rectangle,bgModel,fgModel,1,cv::GC_INIT_WITH_MASK);

    sofa::helper::AdvancedTimer::stepBegin("CVGraphCut") ;

    //std::cout << " mask " << mask << std::endl;
    cv::compare(mask,cv::GC_PR_FGD,maskimg,cv::CMP_EQ);
    // Generate output image
    image.copyTo(foreground,maskimg);
}

