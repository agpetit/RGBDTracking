/*
 * segmentation.cpp
 *
 *  Created on: May 2, 2014
 *      Author: apetit
 */

#include "segmentation.h"

segmentation::segmentation() {
        // TODO Auto-generated constructor stub
neighborhood = 8;
mskt = BBOX;
type = CVGRAPHCUT;
}

segmentation::~segmentation() {
        // TODO Auto-generated destructor stub

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

        if (type == CVGRAPHCUT) mskt = BBOX;
}

void segmentation::segmentationFromRect(cv::Mat &image, cv::Mat &foreground)
{

switch(type){
        case CVGRAPHCUT:
        {
                cv::Mat bgModel,fgModel;

            cv::grabCut(image,    // input image
            mask,   // segmentation result
            rectangle,// rectangle containing foreground
            bgModel,fgModel, // models
            1,        // number of iterations
            cv::GC_INIT_WITH_RECT); // use rectangle
            // Get the pixels marked as likely foreground

            double t = ((double)getTickCount() - t)/getTickFrequency();
            cout << "Times passed in seconds: " << t << endl;

            //std::cout << " mask " << mask << std::endl;
            cv::compare(mask,cv::GC_PR_FGD,maskimg,cv::CMP_EQ);
            // Generate output image
            //cv::Mat foreground(image.size(),CV_8UC3,cv::Scalar(255,255,255));
            image.copyTo(foreground,maskimg); // bg pixels not copied
                break;
        }
        case CUDAGRAPHCUT:
        {
            break;
        }
        }
}

void segmentation::updateMask(cv::Mat &foreground)
{
        switch(mskt){
        case BBOX:
        {
            Point pt;
                std::vector<cv::Point> ptfgd;
                ptfgd.resize(0);

                switch(type){
                case CVGRAPHCUT:
                {
            for(int x = 0; x<mask.cols; x++)
                for(int y = 0; y<mask.rows; y++)
                        if (mask.at<uchar>(y,x) == 1 || mask.at<uchar>(y,x) == 3)
                        {
                                pt.x = x;
                                pt.y = y;
                                ptfgd.push_back(pt);
                        }
                        else if(mask.at<uchar>(y,x) == 0)
                                mask.at<uchar>(y,x) = 2;
            break;
                }
                case CUDAGRAPHCUT:
                {
            break;
                }
                }

            rectangle = cv::boundingRect(ptfgd);
            //std::cout << " ptfgd " << frame_count << " rect1 " << rectangle.x << " " << rectangle.y << " rect2 " << rectangle.width << " " << rectangle.height << std::endl;

            rectangle.x -= 10;
            rectangle.y -= 10;
            rectangle.height += 20;
            rectangle.width += 20;

            std::cout << " rect1 " << rectangle.x << " " << rectangle.y << std::endl;

            for(int x = 0; x<mask.cols; x++)
                for(int y = 0; y<mask.rows; y++)
                        if (x < rectangle.x || x > rectangle.x + rectangle.width || y < rectangle.y || y > rectangle.y + rectangle.height)
                        {
                                mask.at<uchar>(y,x) = 0;
                        }

                break;
        }
        case CONTOUR:
        {
            cv::Mat distanceMap(foreground.size(),CV_8U,cv::Scalar(255));
            cv::Mat dot(foreground.size(),CV_8U,cv::Scalar(0));
            filter(foreground,distanceMap,dot);
        //maskFromDt(distanceMap,mask1);
            trimapFromDt(distanceMap,dot);

                //std::cout << " rect1 " << rectangle.x << " " << rectangle.y << std::endl;

                distImage = distanceMap;
                dotImage = dot;
            break;
        }
        }
}

void segmentation::updateSegmentation(cv::Mat &image,cv::Mat &foreground)
{

        switch(type){
        case CVGRAPHCUT:
        {
    double t = (double)getTickCount();

    cv::Mat bgModel,fgModel;

    cv::Mat _image = image.clone();
    cv::Mat _mask = mask.clone();
    cv::Mat foreground_;

    cv::grabCut(_image,_mask,rectangle,bgModel,fgModel,5,cv::GC_INIT_WITH_MASK);

    // do something ...
    t = ((double)getTickCount() - t)/getTickFrequency();
    cout << "Times passed in seconds: " << t << endl;

    //std::cout << " mask " << mask << std::endl;
    cv::compare(_mask,cv::GC_PR_FGD,maskimg,cv::CMP_EQ);
    // Generate output image
    //cv::Mat foreground(image.size(),CV_8UC3,cv::Scalar(255,255,255));
    _image.copyTo(foreground_,maskimg); // bg pixels not copied
    cv::Mat alpha(image.size(),CV_8UC1,Scalar(0));

    for(int x = 0; x<image.cols; x++)
        for(int y = 0; y<image.rows; y++){
                if (foreground_.at<cv::Vec3b>(y,x)[0] == 0 && foreground_.at<cv::Vec3b>(y,x)[1] == 0 && foreground_.at<cv::Vec3b>(y,x)[2] == 0)
                      alpha.at<uchar>(y,x) = 0;
                else alpha.at<uchar>(y,x) = 255;
        }

    cv::Mat rgb[4];
    cv::split(foreground_,rgb);

    cv::Mat rgba[4]={rgb[0],rgb[1],rgb[2],alpha};
    cv::merge(rgba,4,foreground_);
    foreground = foreground_.clone();

        break;
        }
        case CUDAGRAPHCUT:
        {
            break;
                }
        }

}

void segmentation::updateSegmentationCrop(cv::Mat &image,cv::Mat &foreground)
{

        switch(type){
        case CVGRAPHCUT:
        {
    double t = (double)getTickCount();

    cv::Mat bgModel,fgModel;

    cv::Mat roiimage = image(rectangle);
    cv::Mat roimask = mask(rectangle);

    cv::grabCut(image,mask,rectangle,bgModel,fgModel,1,cv::GC_INIT_WITH_MASK);

    // do something ...
    t = ((double)getTickCount() - t)/getTickFrequency();
    cout << "Times passed in seconds: " << t << endl;

    //std::cout << " mask " << mask << std::endl;
    cv::compare(mask,cv::GC_PR_FGD,maskimg,cv::CMP_EQ);
    // Generate output image
    //cv::Mat foreground(image.size(),CV_8UC3,cv::Scalar(255,255,255));
    image.copyTo(foreground,maskimg); // bg pixels not copied
        break;
        }
        case CUDAGRAPHCUT:
        {

            break;
                }
        }

}

void segmentation::saveResult(const char *filename)
{

}


void segmentation::getResult(cv::Mat &out)
{
}

void segmentation::getResultCrop(cv::Mat &out)
{
}


void segmentation::filter(cv::Mat &out,cv::Mat &dt,cv::Mat &dot)
{

cv::Mat dt0;

//cv::imwrite("out.png",out);


for(int x = 0; x<out.cols; x++)
        for(int y = 0; y<out.rows; y++)
                if (out.at<cv::Vec4b>(y,x)[0] > 0 || out.at<cv::Vec4b>(y,x)[1] > 0 || out.at<cv::Vec4b>(y,x)[2] > 0)
                {
                dot.at<uchar>(y,x) = 0;
                }
                else dot.at<uchar>(y,x) = 255;

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

for(int x = 0; x<_dt.cols; x++)
        for(int y = 0; y<_dt.rows; y++)
        {
                if (_dt.at<uchar>(y,x) > 10)
                {
                mask.at<uchar>(y,x) = 0;
                }
                else
                {
                mask.at<uchar>(y,x) = 255;
                }
        }

mask_ = mask;

/*cv::namedWindow("Mask");
cv::imshow("Mask",mask_);*/

//waitKey(0);
//getchar();

}

void segmentation::trimapFromDt(cv::Mat &_dt,cv::Mat &dot)
{

mask = cv::Mat::zeros(_dt.size(),CV_8U);
Point pt;
std::vector<cv::Point> ptfgd;
ptfgd.resize(0);
cv::Mat strip = dot.clone();
int band = 10;
for(int x = 0; x<_dt.cols; x++)
        for(int y = 0; y<_dt.rows; y++)
        {
                if (_dt.at<uchar>(y,x) > band && dot.at<uchar>(y,x) == 255)
                {
                mask.at<uchar>(y,x) = 0;
                strip.at<uchar>(y,x) = 0;
                }
                else if (_dt.at<uchar>(y,x) > band && dot.at<uchar>(y,x) == 0)
                {
                mask.at<uchar>(y,x) = 2;
                pt.x = x;
                pt.y = y;
                ptfgd.push_back(pt);
                strip.at<uchar>(y,x) = 255;

                }
                else
                {
                mask.at<uchar>(y,x) = 1;
                pt.x = x;
                pt.y = y;
                ptfgd.push_back(pt);
                strip.at<uchar>(y,x) = 127;
                }
        }

std::cout << " ok bd rect 0 " << ptfgd.size() <<  std::endl;
if (ptfgd.size() >= 200)
rectangle = cv::boundingRect(ptfgd);

//cv::imwrite("strip.png",strip);
//cv::imwrite("dot.png",dot);



//std::cout << " ok bd rect 1 " << rectangle.x << " " << rectangle.y << std::endl;


if (ptfgd.size() < 200)
{
        std::cout << " fail " << std::endl;
        for(int x = 0; x<_dt.cols; x++)
                for(int y = 0; y<_dt.rows; y++)
                {
                        if (_dt.at<uchar>(y,x) > band && dot.at<uchar>(y,x) == 255)
                        {
                        mask.at<uchar>(y,x) = 1;
                        }
                        else if (_dt.at<uchar>(y,x) > band && dot.at<uchar>(y,x) == 0)
                        {
                        mask.at<uchar>(y,x) = 2;
                        }
                        else
                        {
                        mask.at<uchar>(y,x) = 1;
                        }
                }
}
else {
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

bool segmentation::verifyResult(const char *filename)
{

        return false;
}

FIBITMAP* segmentation::convertCVFree(cv::Mat &in)
{

        FIBITMAP* out = NULL;
        width  = in.size().width;
        height = in.size().height;

        switch(in.type())
        {
        case CV_8U  :{out = FreeImage_AllocateT(FIT_BITMAP,width, height, 8) ;}break;  // 8  bit grayscale
        case CV_8UC3:{out = FreeImage_AllocateT(FIT_BITMAP,width, height, 24);}break;  // 24 bit RGB
        case CV_16U :{out = FreeImage_AllocateT(FIT_UINT16,width, height, 16);}break;  // 16 bit grayscale
        case CV_16S :{out = FreeImage_AllocateT(FIT_INT16 ,width, height, 16);}break;
        case CV_32S :{out = FreeImage_AllocateT(FIT_INT32 ,width, height, 32);}break;
        case CV_32F :{out = FreeImage_AllocateT(FIT_FLOAT ,width, height, 32);}break;
        case CV_64F :{out = FreeImage_AllocateT(FIT_DOUBLE,width, height, 32);}break;
        }

        //if(out==NULL)
                //return FALSE;

        int srcRowBytes = width  * in.elemSize();

        for (int ih=0;ih<height;ih++)
        {
                BYTE* ptr2Line = FreeImage_GetScanLine(out,(height-1)-ih);
                memcpy(ptr2Line,in.ptr(ih),srcRowBytes);
        }

        return out;

}

void segmentation::clean()
{
}
