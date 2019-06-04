/*
 * segmentation.cpp
 *
 *  Created on: May 2, 2014
 *      Author: apetit
 */

#include "CUDASegmentation.h"

#include <sofa/helper/AdvancedTimer.h>

#ifdef HAVECUDA
// Functions from GrabcutUtil.cu
cudaError_t TrimapFromRect(Npp8u *alpha, int alpha_pitch, NppiRect rect, int width, int height);
cudaError_t TrimapFromMask(Npp8u *alpha, int alpha_pitch, Npp8u *dt, int width, int height);
cudaError_t ApplyMatte(int mode, uchar4 *result, int result_pitch, const uchar4 *image, int image_pitch, const unsigned char *matte, int matte_pitch, int width, int height);
#endif

CUDASegmentation::CUDASegmentation() {
    // TODO Auto-generated constructor stub
    mskt = BBOX;
}

CUDASegmentation::~CUDASegmentation() {
        // TODO Auto-generated destructor stub
#ifdef HAVECUDA
checkCudaErrors(cudaFree(d_image));
checkCudaErrors(cudaFree(d_trimap));
checkCudaErrors(cudaFree(d_dt));
checkCudaErrors(cudaFree(d_crop_image));
checkCudaErrors(cudaFree(d_crop_trimap));
checkCudaErrors(cudaFree(d_crop_dt));

delete cudaseg;
#endif
//delete pbo_resource;
}

void CUDASegmentation::segmentationFromRect(cv::Mat &image, cv::Mat &foreground)
{
#ifdef HAVECUDA

   FIBITMAP* _dib = NULL;
   FIBITMAP* _dib4 = NULL;

   if(_dib)  {
   // get rid of the current dib.
       FreeImage_Unload(_dib);
   }

   //_dib = convertCVFree(image);

   width  = image.size().width;
   height = image.size().height;

   switch(image.type()) {
       case CV_8U  :
           _dib = FreeImage_AllocateT(FIT_BITMAP,width, height, 8) ;
           break;  // 8  bit grayscale
       case CV_8UC3:
           _dib = FreeImage_AllocateT(FIT_BITMAP,width, height, 24);
           break;  // 24 bit RGB
       case CV_16U :
           _dib = FreeImage_AllocateT(FIT_UINT16,width, height, 16);
           break;  // 16 bit grayscale
       case CV_16S :
           _dib = FreeImage_AllocateT(FIT_INT16 ,width, height, 16);
           break;
       case CV_32S :
           _dib = FreeImage_AllocateT(FIT_INT32 ,width, height, 32);
           break;
       case CV_32F :
           _dib = FreeImage_AllocateT(FIT_FLOAT ,width, height, 32);
           break;
       case CV_64F :
           _dib = FreeImage_AllocateT(FIT_DOUBLE,width, height, 32);
           break;
   }

   //if(out==NULL)
       //return FALSE;

   int srcRowBytes = width  * image.elemSize();

   for (int ih=0;ih<height;ih++) {
       BYTE* ptr2Line = FreeImage_GetScanLine(_dib,(height-1)-ih);
       memcpy(ptr2Line,image.ptr(ih),srcRowBytes);
   }

   width = FreeImage_GetWidth(_dib);
   height = FreeImage_GetHeight(_dib);

   if (width > MAX_IMAGE_SIZE || height > MAX_IMAGE_SIZE) {

       float scale_factor = min(MAX_IMAGE_SIZE / width, MAX_IMAGE_SIZE / height);

       FIBITMAP *pResampled = FreeImage_Rescale(
           _dib,
           (int)(scale_factor * width),
           (int)(scale_factor * height),
           FILTER_BICUBIC);
       width = FreeImage_GetWidth(pResampled);
       height = FreeImage_GetHeight(pResampled);

       _dib4 = FreeImage_ConvertTo32Bits(pResampled);

       FreeImage_Unload(pResampled);

   } else {
       _dib4 = FreeImage_ConvertTo32Bits(_dib);
   }

   /*if (!g_bQATest && g_bDisplay)
   {
       initGL(&argc, argv, width, height);
   }*/

   std::cout << " width " << width << " height " << height << std::endl;

   checkCudaErrors(cudaMallocPitch(&d_image, &image_pitch, width * sizeof(uchar4), height));
   checkCudaErrors(cudaMemcpy2D(d_image, image_pitch, FreeImage_GetBits(_dib4) , FreeImage_GetPitch(_dib4), width * sizeof(uchar4), height, cudaMemcpyHostToDevice));

   FreeImage_Unload(_dib);
   FreeImage_Unload(_dib4);

   checkCudaErrors(cudaMallocPitch(&d_trimap, &trimap_pitch, width, height));

   // Setup GrabCut
   cudaseg = new cudaSegmentation(d_image, (int) image_pitch, d_trimap, (int) trimap_pitch, width, height);

   // Default selection rectangle
   rect.x = (int) ceil(width * 0.1);
   rect.y = (int) ceil(height * 0.1);
   rect.width = width - 2 * rect.x;
   rect.height = height - 2 * rect.y;

   rect.x = rectangle.x;
   rect.y = rectangle.y;
   rect.width = rectangle.width;
   rect.height = rectangle.height;

   checkCudaErrors(TrimapFromRect(d_trimap, (int) trimap_pitch, rect, width, height));

   cudaseg->computeSegmentationFromTrimap();

   /*if (!g_bQATest && g_bDisplay)
   {
       glutMainLoop();
   }*/

   int qaStatus = EXIT_SUCCESS;

   /*if (g_bQATest) {
       qaStatus = verifyResult(sdkFindFilePath((char *)sReferenceFile.c_str(), argv[0])) ? EXIT_SUCCESS : EXIT_FAILURE;
   }*/

   getResult(foreground);

   // Cleanup
   //delete cudaseg;

   /*checkCudaErrors(cudaFree(d_image));
   checkCudaErrors(cudaFree(d_trimap));*/
   #endif
}

void CUDASegmentation::updateMaskBBox(cv::Mat &foreground, std::vector<cv::Point> & ptfgd)
{
    for(int x = 0; x<mask.cols; x++) {
        for(int y = 0; y<mask.rows; y++) {
            if (foreground.at<cv::Vec4b>(y,x)[0] > 0 ||
                foreground.at<cv::Vec4b>(y,x)[1] > 0 ||
                foreground.at<cv::Vec4b>(y,x)[2] > 0
            ) {
                    cv::Point pt (x, y);
                    //pt.x = x;
                    //pt.y = y;
                    ptfgd.push_back(pt);
            }
        }
    }
}

void CUDASegmentation::updateSegmentation(cv::Mat &foreground, cv::Mat &image)
{
   #ifdef HAVECUDA
// this code might be duplicate, function should be refactored

    sofa::helper::AdvancedTimer::stepBegin("CUDAGraphCut") ;

    FIBITMAP* _dib = NULL;
    FIBITMAP* _dib4 = NULL;

    if(_dib)  // get rid of the current dib.
       FreeImage_Unload(_dib);

    //_dib = convertCVFree(image);

    width  = image.size().width;
    height = image.size().height;

    switch(image.type()) {
       case CV_8U  :
           _dib = FreeImage_AllocateT(FIT_BITMAP,width, height, 8) ;
           break;  // 8  bit grayscale
       case CV_8UC3:
           _dib = FreeImage_AllocateT(FIT_BITMAP,width, height, 24);
           break;  // 24 bit RGB
       case CV_16U :
           _dib = FreeImage_AllocateT(FIT_UINT16,width, height, 16);
           break;  // 16 bit grayscale
       case CV_16S :
           _dib = FreeImage_AllocateT(FIT_INT16 ,width, height, 16);
           break;
       case CV_32S :
           _dib = FreeImage_AllocateT(FIT_INT32 ,width, height, 32);
           break;
       case CV_32F :
           _dib = FreeImage_AllocateT(FIT_FLOAT ,width, height, 32);
           break;
       case CV_64F :
           _dib = FreeImage_AllocateT(FIT_DOUBLE,width, height, 32);
           break;
    }


    //if(out==NULL)
       //return FALSE;

    int srcRowBytes = width  * image.elemSize();

    for (int ih=0;ih<height;ih++) {
       BYTE* ptr2Line = FreeImage_GetScanLine(_dib,(height-1)-ih);
       memcpy(ptr2Line,image.ptr(ih),srcRowBytes);
    }

    width = FreeImage_GetWidth(_dib);
    height = FreeImage_GetHeight(_dib);


    if (width > MAX_IMAGE_SIZE || height > MAX_IMAGE_SIZE) {

       float scale_factor = min(MAX_IMAGE_SIZE / width, MAX_IMAGE_SIZE / height);

       FIBITMAP *pResampled = FreeImage_Rescale(_dib, (int)(scale_factor * width), (int)(scale_factor * height), FILTER_BICUBIC);
       width = FreeImage_GetWidth(pResampled);
       height = FreeImage_GetHeight(pResampled);

       _dib4 = FreeImage_ConvertTo32Bits(pResampled);

       FreeImage_Unload(pResampled);
    } else {
       _dib4 = FreeImage_ConvertTo32Bits(_dib);
    }

    /*if (!g_bQATest && g_bDisplay)
    {
       initGL(&argc, argv, width, height);
    }*/

    //std::cout << " width " << width << " height " << height << std::endl;

    checkCudaErrors(cudaMallocPitch(&d_image, &image_pitch, width * sizeof(uchar4), height));
    checkCudaErrors(cudaMemcpy2D(d_image, image_pitch, FreeImage_GetBits(_dib4) , FreeImage_GetPitch(_dib4), width * sizeof(uchar4), height, cudaMemcpyHostToDevice));

    FreeImage_Unload(_dib);
    FreeImage_Unload(_dib4);

    checkCudaErrors(cudaMallocPitch(&d_trimap, &trimap_pitch, width, height));

    switch(mskt){
    case BBOX:
       rect.x = rectangle.x;
       rect.y = rectangle.y;
       rect.width = rectangle.width;
       rect.height = rectangle.height;

       //checkCudaErrors(TrimapFromMask(d_trimap, (int) trimap_pitch, d_dt, width, height));
       checkCudaErrors(TrimapFromRect(d_trimap, (int) trimap_pitch, rect, width, height));
       cudaseg->updateTrimap(d_trimap);
       break;
    case CONTOUR:
       //FIBITMAP *h_tri = FreeImage_Allocate(width, height, 8);

       //checkCudaErrors(cudaMemcpy2D(FreeImage_GetBits(h_tri) , FreeImage_GetPitch(h_tri), d_trimap, trimap_pitch, width, height, cudaMemcpyDeviceToHost));

       //checkCudaErrors(cudaMemcpy2D(d_trimap, trimap_pitch, FreeImage_GetBits(tri) , FreeImage_GetPitch(tri), width * sizeof(uchar), height, cudaMemcpyDeviceToHost));

       FIBITMAP* _distanceMap = NULL;

       if(_distanceMap)  {
       // get rid of the current dib.
           FreeImage_Unload(_distanceMap);
       }

       _distanceMap = convertCVFree(mask);

       size_t dt_pitch;

       checkCudaErrors(cudaMallocPitch(&d_dt, &dt_pitch, width * sizeof(uchar), height));
       checkCudaErrors(cudaMemcpy2D(d_dt, dt_pitch, FreeImage_GetBits(_distanceMap) , FreeImage_GetPitch(_distanceMap), width * sizeof(uchar), height, cudaMemcpyHostToDevice));

       cudaseg->updateTrimap(d_dt);

       //FreeImage_Unload(h_tri);

       FreeImage_Unload(_distanceMap);
       break;
    }

    checkCudaErrors(cudaFree(d_dt));

    //cudaseg->computeSegmentationFromTrimap();
    cudaseg->updateImage(d_image);

    // Setup GrabCut
    //cudaseg = new cudaSegmentation(d_image, (int) image_pitch, d_trimap, (int) trimap_pitch, width, height);

    //cudaseg->computeSegmentationFromTrimap();
    cudaseg->updateSegmentation();

    /*if (!g_bQATest && g_bDisplay) {
       glutMainLoop();
    }*/

    int qaStatus = EXIT_SUCCESS;

    /*if (g_bQATest) {
       qaStatus = verifyResult(sdkFindFilePath((char *)sReferenceFile.c_str(), argv[0])) ? EXIT_SUCCESS : EXIT_FAILURE;
    }*/
    getResult(foreground);

    //cv::imwrite("foreground01.png", foreground.clone());

    //waitKey(0);
    //getchar();

    // Cleanup
    //delete cudaseg;

    checkCudaErrors(cudaFree(d_image));
    checkCudaErrors(cudaFree(d_trimap));

    sofa::helper::AdvancedTimer::stepEnd("CUDAGaphCut") ;

    #endif
}

void CUDASegmentation::updateSegmentationCrop(cv::Mat &foreground, cv::Mat &image)
{
   #ifdef HAVECUDA
   cv::Mat crop_image = image(rectangle);
   cv::Mat crop_mask = mask(rectangle);

   sofa::helper::AdvancedTimer::stepBegin("CUDAGraphCut") ;

   FIBITMAP* _dib = NULL;
   FIBITMAP* _dib4 = NULL;

   FIBITMAP* _crop_dib = NULL;
   FIBITMAP* _crop_dib4 = NULL;

   if(_dib)   {
       // get rid of the current dib.
       FreeImage_Unload(_dib);
       FreeImage_Unload(_crop_dib);
   }

   //_dib = convertCVFree(image);

   width  = image.size().width;
   height = image.size().height;

   crop_width  = crop_image.size().width;
   crop_height = crop_image.size().height;

   //std::cout << " crop_w " << crop_width << " " << crop_height << std::endl;

   switch(image.type()) {
       case CV_8U  :
           _dib = FreeImage_AllocateT(FIT_BITMAP,width, height, 8) ;
           _crop_dib = FreeImage_AllocateT(FIT_BITMAP,crop_width, crop_height, 8) ;
           break;  // 8  bit grayscale
       case CV_8UC3:
           _dib = FreeImage_AllocateT(FIT_BITMAP,width, height, 24);
           _crop_dib = FreeImage_AllocateT(FIT_BITMAP,crop_width, crop_height, 24) ;
           break;  // 24 bit RGB
       case CV_16U :
           _dib = FreeImage_AllocateT(FIT_UINT16,width, height, 16);
           _crop_dib = FreeImage_AllocateT(FIT_UINT16,crop_width, crop_height, 16);
           break;  // 16 bit grayscale
       case CV_16S :
           _dib = FreeImage_AllocateT(FIT_INT16 ,width, height, 16);
           _crop_dib = FreeImage_AllocateT(FIT_INT16 ,crop_width, crop_height, 16);
           break;
       case CV_32S :
           _dib = FreeImage_AllocateT(FIT_INT32 ,width, height, 32);
           _crop_dib = FreeImage_AllocateT(FIT_INT32 ,crop_width, crop_height, 32);
           break;
       case CV_32F :
           _dib = FreeImage_AllocateT(FIT_FLOAT ,width, height, 32);
           _crop_dib = FreeImage_AllocateT(FIT_FLOAT ,crop_width, crop_height, 32);
           break;
       case CV_64F :
           _dib = FreeImage_AllocateT(FIT_DOUBLE,width, height, 32);
           _crop_dib = FreeImage_AllocateT(FIT_DOUBLE,crop_width, crop_height, 32);
           break;
   }

   //if(out==NULL)
       //return FALSE;

   int srcRowBytes = width  * image.elemSize();

   int crop_srcRowBytes = crop_width  * crop_image.elemSize();

   for (int ih=0;ih<height;ih++) {
       BYTE* ptr2Line = FreeImage_GetScanLine(_dib,(height-1)-ih);
       memcpy(ptr2Line,image.ptr(ih),srcRowBytes);
   }

   for (int ih=0;ih<crop_height;ih++) {
       BYTE* crop_ptr2Line = FreeImage_GetScanLine(_crop_dib,(crop_height-1)-ih);
       memcpy(crop_ptr2Line,crop_image.ptr(ih),crop_srcRowBytes);
   }


   if (width > MAX_IMAGE_SIZE || height > MAX_IMAGE_SIZE) {

       float scale_factor = min(MAX_IMAGE_SIZE / width, MAX_IMAGE_SIZE / height);

       FIBITMAP *pResampled = FreeImage_Rescale(_dib, (int)(scale_factor * width), (int)(scale_factor * height), FILTER_BICUBIC);
       width = FreeImage_GetWidth(pResampled);
       height = FreeImage_GetHeight(pResampled);

       _dib4 = FreeImage_ConvertTo32Bits(pResampled);

       FreeImage_Unload(pResampled);
   } else {
       _dib4 = FreeImage_ConvertTo32Bits(_dib);
   }

   if (crop_width > MAX_IMAGE_SIZE || crop_height > MAX_IMAGE_SIZE)
   {

       float scale_factor = min(MAX_IMAGE_SIZE / crop_width, MAX_IMAGE_SIZE / crop_height);

       FIBITMAP *crop_pResampled = FreeImage_Rescale(_crop_dib, (int)(scale_factor * crop_width), (int)(scale_factor * crop_height), FILTER_BICUBIC);
       crop_width = FreeImage_GetWidth(crop_pResampled);
       crop_height = FreeImage_GetHeight(crop_pResampled);

       _crop_dib4 = FreeImage_ConvertTo32Bits(crop_pResampled);

       FreeImage_Unload(crop_pResampled);

   } else {
       _crop_dib4 = FreeImage_ConvertTo32Bits(_crop_dib);
   }

   /*if (!g_bQATest && g_bDisplay) {
       initGL(&argc, argv, width, height);
   }*/

   std::cout << " width " << width << " height " << height << std::endl;

   checkCudaErrors(cudaMallocPitch(&d_image, &image_pitch, width * sizeof(uchar4), height));
   checkCudaErrors(cudaMemcpy2D(d_image, image_pitch, FreeImage_GetBits(_dib4) , FreeImage_GetPitch(_dib4), width * sizeof(uchar4), height, cudaMemcpyHostToDevice));

   checkCudaErrors(cudaMallocPitch(&d_crop_image, &crop_image_pitch, crop_width * sizeof(uchar4), crop_height));
   checkCudaErrors(cudaMemcpy2D(d_crop_image, crop_image_pitch, FreeImage_GetBits(_crop_dib4) , FreeImage_GetPitch(_crop_dib4), crop_width * sizeof(uchar4), crop_height, cudaMemcpyHostToDevice));

   //FreeImage_Save(FIF_PNG, _crop_dib4, "img_crop.png", 0);

   FreeImage_Unload(_dib);
   FreeImage_Unload(_dib4);

   FreeImage_Unload(_crop_dib);
   FreeImage_Unload(_crop_dib4);

   //checkCudaErrors(cudaMallocPitch(&d_trimap, &trimap_pitch, width, height));
   //checkCudaErrors(cudaMallocPitch(&d_crop_trimap, &crop_trimap_pitch, crop_width, crop_height));

   switch(mskt){
       case BBOX:
           rect.x = rectangle.x;
           rect.y = rectangle.y;
           rect.width = rectangle.width;
           rect.height = rectangle.height;

           //checkCudaErrors(TrimapFromMask(d_trimap, (int) trimap_pitch, d_dt, width, height));
           checkCudaErrors(TrimapFromRect(d_trimap, (int) trimap_pitch, rect, width, height));
           cudaseg->updateTrimap(d_trimap);
           break;

       case CONTOUR:
           //FIBITMAP *h_tri = FreeImage_Allocate(width, height, 8);
           //checkCudaErrors(cudaMemcpy2D(FreeImage_GetBits(h_tri) , FreeImage_GetPitch(h_tri), d_trimap, trimap_pitch, width, height, cudaMemcpyDeviceToHost));
           //checkCudaErrors(cudaMemcpy2D(d_trimap, trimap_pitch, FreeImage_GetBits(tri) , FreeImage_GetPitch(tri), width * sizeof(uchar), height, cudaMemcpyDeviceToHost));

           FIBITMAP* _distanceMap = NULL;
           FIBITMAP* _crop_distanceMap = NULL;

           if(_distanceMap)  {
           // get rid of the current dib.
               FreeImage_Unload(_distanceMap);
           }

           if(_crop_distanceMap) {
           // get rid of the current dib.
               FreeImage_Unload(_crop_distanceMap);
           }

           _distanceMap = convertCVFree(mask);
           _crop_distanceMap = convertCVFree(crop_mask);

           size_t dt_pitch;
           size_t dt_crop_pitch;

           //FreeImage_Save(FIF_PNG, _crop_distanceMap, "crop_dmap.png", 0);
           //FreeImage_Save(FIF_PNG, _distanceMap, "dmap.png", 0);

           checkCudaErrors(cudaMallocPitch(&d_dt, &dt_pitch, width * sizeof(uchar), height));
           checkCudaErrors(cudaMemcpy2D(d_dt, dt_pitch, FreeImage_GetBits(_distanceMap) , FreeImage_GetPitch(_distanceMap), width * sizeof(uchar), height, cudaMemcpyHostToDevice));

           checkCudaErrors(cudaMallocPitch(&d_crop_dt, &dt_crop_pitch, crop_width * sizeof(uchar), crop_height));
           checkCudaErrors(cudaMemcpy2D(d_crop_dt, dt_crop_pitch, FreeImage_GetBits(_crop_distanceMap) , FreeImage_GetPitch(_crop_distanceMap), crop_width * sizeof(uchar), crop_height, cudaMemcpyHostToDevice));

           cudaseg->updateTrimapCrop(d_dt,d_crop_dt,(int)dt_crop_pitch);

           //FreeImage_Unload(h_tri);

           FreeImage_Unload(_distanceMap);
           FreeImage_Unload(_crop_distanceMap);
           break;
   }

   checkCudaErrors(cudaFree(d_dt));
   checkCudaErrors(cudaFree(d_crop_dt));

   //cudaseg->computeSegmentationFromTrimap();
   cudaseg->updateImage(d_image);
   cudaseg->updateImageCrop(d_crop_image, (int)crop_image_pitch, crop_width, crop_height);

   // Setup GrabCut
   //cudaseg = new cudaSegmentation(d_image, (int) image_pitch, d_trimap, (int) trimap_pitch, width, height);

   //cudaseg->computeSegmentationFromTrimap();
   cudaseg->updateSegmentationCrop();

   /*if (!g_bQATest && g_bDisplay) {
       glutMainLoop();
   }*/

   int qaStatus = EXIT_SUCCESS;

   /*if (g_bQATest) {
       qaStatus = verifyResult(sdkFindFilePath((char *)sReferenceFile.c_str(), argv[0])) ? EXIT_SUCCESS : EXIT_FAILURE;
   }*/

   cv::Mat crop_foreground,foreground0;
   getResultCrop(crop_foreground);
   foreground0 = cv::Mat::zeros(image.size(),CV_8UC4);
   cv::Mat roiforeground = foreground0(rectangle);
   crop_foreground.copyTo(roiforeground); // bg pixels not copied
   foreground = foreground0;

   //waitKey(0);
   //getchar();

   // Cleanup
   //delete cudaseg;

   checkCudaErrors(cudaFree(d_image));
   //checkCudaErrors(cudaFree(d_trimap));

       //checkCudaErrors(cudaFree(d_crop_image));
   //checkCudaErrors(cudaFree(d_crop_trimap));

   sofa::helper::AdvancedTimer::stepEnd("CUDAGraphCut") ;
   #endif
}

void CUDASegmentation::saveResult(const char *filename) {
#ifdef HAVECUDA
    uchar4 *d_result;
    size_t result_pitch;

    checkCudaErrors(cudaMallocPitch(&d_result, &result_pitch, width*4, height));

    std::cout << " wi " << width << " he " << height << std::endl;

    ApplyMatte(2, d_result, (int) result_pitch, d_image, (int) image_pitch, cudaseg->getAlpha(), cudaseg->getAlphaPitch(), width, height);

    FIBITMAP *h_Image = FreeImage_Allocate(width, height, 32);

    checkCudaErrors(cudaMemcpy2D(FreeImage_GetBits(h_Image) , FreeImage_GetPitch(h_Image), d_result, result_pitch, width * 4, height, cudaMemcpyDeviceToHost));

    FreeImage_Save(FIF_PNG, h_Image, filename, 0);

    FreeImage_Unload(h_Image);

    checkCudaErrors(cudaFree(d_result));

    printf("Saved result as %s\n", filename);
#endif
}

void CUDASegmentation::getResult(cv::Mat &out)
{
#ifdef HAVECUDA
    uchar4 *d_result;
    size_t result_pitch;

    checkCudaErrors(cudaMallocPitch(&d_result, &result_pitch, width*4, height));

    ApplyMatte(2, d_result, (int) result_pitch, d_image, (int) image_pitch, cudaseg->getAlpha(), cudaseg->getAlphaPitch(), width, height);

    FIBITMAP *h_Image = FreeImage_Allocate(width, height, 32);

    checkCudaErrors(cudaMemcpy2D(FreeImage_GetBits(h_Image) , FreeImage_GetPitch(h_Image), d_result, result_pitch, width * 4, height, cudaMemcpyDeviceToHost));

    //convertFreeCV(h_Image,out);

    IplImage *l_pCvImg = NULL;
    CvSize l_size;
    int x, y;
    unsigned int l_pitch;
    FREE_IMAGE_TYPE l_image_type;

    // FreeImage is bottom up ordered but not OpenCV
    FreeImage_FlipVertical(h_Image);

    // Get image property
    l_size.width = FreeImage_GetWidth(h_Image);
    l_size.height = FreeImage_GetHeight(h_Image);
    l_pitch = FreeImage_GetPitch(h_Image);
    l_image_type = FreeImage_GetImageType(h_Image);

    //std::cout << " width " << l_size.width << " height " << l_size.height << std::endl;

    // Test image type
    if((l_image_type == FIT_BITMAP) && (FreeImage_GetBPP(h_Image) == 32)) {
        // Create OpenCV image
        l_pCvImg = cvCreateImage(l_size, IPL_DEPTH_8U, 4);

        BYTE *l_pLine = (BYTE*)FreeImage_GetBits(h_Image);
        for(y = 0; y < l_size.height; y++) {
            BYTE *l_pPixel = (BYTE*)l_pLine;
            for(x = 0; x < l_size.width; x++) {
                // Get pixel addresses
                uchar* b = &((uchar*)(l_pCvImg->imageData + l_pCvImg->widthStep*y))[x*4];
                uchar* g = &((uchar*)(l_pCvImg->imageData + l_pCvImg->widthStep*y))[x*4+1];
                uchar* r = &((uchar*)(l_pCvImg->imageData + l_pCvImg->widthStep*y))[x*4+2];
                uchar* a = &((uchar*)(l_pCvImg->imageData + l_pCvImg->widthStep*y))[x*4+3];

                // Copy pixel value
                *b = l_pPixel[FI_RGBA_BLUE];
                *r = l_pPixel[FI_RGBA_RED];
                *g = l_pPixel[FI_RGBA_GREEN];
                *a = l_pPixel[FI_RGBA_ALPHA];

                // Next pixel in the line
                l_pPixel += 4;
            }
            // Next line
            l_pLine += l_pitch;
        }
    }

    out = cv::cvarrToMat(l_pCvImg);

    FreeImage_Unload(h_Image);

    checkCudaErrors(cudaFree(d_result));
#endif
}

void CUDASegmentation::getResultCrop(cv::Mat &out)
{
#ifdef HAVECUDA
    uchar4 *d_crop_result;
    size_t crop_result_pitch;

    checkCudaErrors(cudaMallocPitch(&d_crop_result, &crop_result_pitch, crop_width*4, crop_height));

    ApplyMatte(2, d_crop_result, (int) crop_result_pitch, d_crop_image, (int) crop_image_pitch, cudaseg->getAlphaCrop(), cudaseg->getAlphaPitchCrop(), crop_width, crop_height);

    FIBITMAP *h_crop_Image = FreeImage_Allocate(crop_width, crop_height, 32);

    checkCudaErrors(cudaMemcpy2D(FreeImage_GetBits(h_crop_Image) , FreeImage_GetPitch(h_crop_Image), d_crop_result, crop_result_pitch, crop_width * 4, crop_height, cudaMemcpyDeviceToHost));

    //convertFreeCV(h_Image,out);

    IplImage *l_pCvImg = NULL;
    CvSize l_size;
    int x, y;
    unsigned int l_pitch;
    FREE_IMAGE_TYPE l_image_type;

    // FreeImage is bottom up ordered but not OpenCV
    FreeImage_FlipVertical(h_crop_Image);

    // Get image property
    l_size.width = FreeImage_GetWidth(h_crop_Image);
    l_size.height = FreeImage_GetHeight(h_crop_Image);
    l_pitch = FreeImage_GetPitch(h_crop_Image);
    l_image_type = FreeImage_GetImageType(h_crop_Image);

    //std::cout << " width " << l_size.width << " height " << l_size.height << std::endl;

    // Test image type
    if((l_image_type == FIT_BITMAP) && (FreeImage_GetBPP(h_crop_Image) == 32)) {
        // Create OpenCV image
        l_pCvImg = cvCreateImage(l_size, IPL_DEPTH_8U, 4);

        BYTE *l_pLine = (BYTE*)FreeImage_GetBits(h_crop_Image);
        for(y = 0; y < l_size.height; y++) {
            BYTE *l_pPixel = (BYTE*)l_pLine;
            for(x = 0; x < l_size.width; x++) {
                // Get pixel addresses
                uchar* b = &((uchar*)(l_pCvImg->imageData + l_pCvImg->widthStep*y))[x*4];
                uchar* g = &((uchar*)(l_pCvImg->imageData + l_pCvImg->widthStep*y))[x*4+1];
                uchar* r = &((uchar*)(l_pCvImg->imageData + l_pCvImg->widthStep*y))[x*4+2];
                uchar* a = &((uchar*)(l_pCvImg->imageData + l_pCvImg->widthStep*y))[x*4+3];

                // Copy pixel value
                *b = l_pPixel[FI_RGBA_BLUE];
                *r = l_pPixel[FI_RGBA_RED];
                *g = l_pPixel[FI_RGBA_GREEN];
                *a = l_pPixel[FI_RGBA_ALPHA];


                // Next pixel in the line
                l_pPixel += 4;
            }
            // Next line
            l_pLine += l_pitch;
        }
    }

    cv::Mat out0 = cv::cvarrToMat(l_pCvImg);
    out = out0;

    FreeImage_Unload(h_crop_Image);

    checkCudaErrors(cudaFree(d_crop_result));
#endif
}

bool CUDASegmentation::verifyResult(const char *filename)
{
#ifdef HAVECUDA
    uchar4 *d_result;
    size_t result_pitch;

    FREE_IMAGE_FORMAT eFormat = FreeImage_GetFileType(filename);
    FIBITMAP *goldImage = FreeImage_Load(eFormat, filename);

    if (goldImage == 0) {
        printf("Could not open gold image: %s\n", filename);
        return false;
    }

    if (FreeImage_GetHeight(goldImage) != (unsigned int)height ||
        FreeImage_GetWidth(goldImage) != (unsigned int)width
    ) {
        printf("Gold image size != result image size\n");
        FreeImage_Unload(goldImage);
        return false;
    }

    FreeImage_Save(FIF_PNG, goldImage, "gold.png", 0);

    checkCudaErrors(cudaMallocPitch(&d_result, &result_pitch, width*4, height));

    ApplyMatte(2, d_result, (int) result_pitch, d_image, (int) image_pitch, cudaseg->getAlpha(), cudaseg->getAlphaPitch(), width, height);

    FIBITMAP *h_Image = FreeImage_Allocate(width, height, 32);
    checkCudaErrors(cudaMemcpy2D(FreeImage_GetBits(h_Image) , FreeImage_GetPitch(h_Image), d_result, result_pitch, width * 4, height, cudaMemcpyDeviceToHost));


    bool result = true;

    int bytespp = FreeImage_GetLine(h_Image) / FreeImage_GetWidth(h_Image);

    for (int y = 0; y < height; y++) {
        BYTE *goldBits = FreeImage_GetScanLine(goldImage, y);
        BYTE *resultBits = FreeImage_GetScanLine(h_Image, y);

        for (int x = 0; x < width * bytespp; x++) {
            if (goldBits[x] != resultBits[x]) {
                result = false;
            }
        }
    }

    printf("Checking grabcut results with reference file <%s>\n", filename);
    printf("Images %s\n", result ? "Match!" : "Mismatched!");

    FreeImage_Unload(h_Image);
    FreeImage_Unload(goldImage);
    checkCudaErrors(cudaFree(d_result));

    return result;
#endif
}
#ifdef HAVECUDA
FIBITMAP* CUDASegmentation::convertCVFree(cv::Mat &in)
{

    FIBITMAP* out = NULL;
    width  = in.size().width;
    height = in.size().height;

    switch(in.type()) {
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

    for (int ih=0;ih<height;ih++) {
        BYTE* ptr2Line = FreeImage_GetScanLine(out,(height-1)-ih);
        memcpy(ptr2Line,in.ptr(ih),srcRowBytes);
    }

    return out;
}
#endif

void CUDASegmentation::clean()
{
#ifdef HAVECUDA
delete cudaseg;

checkCudaErrors(cudaFree(d_image));
checkCudaErrors(cudaFree(d_trimap));
checkCudaErrors(cudaFree(d_dt));
#endif
}

