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

#define SOFA_RGBDTRACKING_IMAGECONVERTER_INL

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/Mapping.inl>
#include <sofa/simulation/Simulation.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/gui/BaseGUI.h>
#include <sofa/gui/BaseViewer.h>
#include <sofa/gui/GUIManager.h>

#ifdef USING_OMP_PRAGMAS
#include <omp.h>
#endif

#include <limits>
#include <iterator>
#include <sofa/helper/gl/Color.h>

#ifdef Success
  #undef Success
#endif

#include <pcl/common/common_headers.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/registration/icp.h>
#include <pcl/common/transforms.h>

#include "ImageConverter.h"

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

template <class DataTypes, class DepthTypes>
ImageConverter<DataTypes, DepthTypes>::ImageConverter()
    : Inherit()
    , depthImage(initData(&depthImage,DepthTypes(),"depthImage","depth map"))
    //, depthTransform(initData(&depthTransform, TransformType(), "depthTransform" , ""))
    , image(initData(&image,ImageTypes(),"image","image"))
	//, transform(initData(&transform, TransformType(), "transform" , ""))
	, useRealData(initData(&useRealData,true,"useRealData","Use real data"))
	, useSensor(initData(&useSensor,false,"useSensor","Use the sensor"))
	, sensorType(initData(&sensorType, 0,"sensorType","Type of the sensor"))
	, niterations(initData(&niterations,3,"niterations","Number of iterations in the tracking process"))
{
	//softk.init();
	this->f_listening.setValue(true); 
	this->addAlias(&depthImage, "depthImage");
	depthImage.setGroup("depthImage");
	depthImage.setReadOnly(true); 
	this->addAlias(&image, "image");
	image.setGroup("image");
	image.setReadOnly(true);
  	//niterations.setValue(2);

}

template <class DataTypes, class DepthTypes>
ImageConverter<DataTypes, DepthTypes>::~ImageConverter()
{
}

template <class DataTypes, class DepthTypes>
void ImageConverter<DataTypes, DepthTypes>::init()
{
    this->Inherit::init();
    core::objectmodel::BaseContext* context = this->getContext();

	mstate = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());

		std::cout << " k init " << std::endl;
		
}

template<class DataTypes, class DepthTypes>
void ImageConverter<DataTypes, DepthTypes>::getImages()
{    
	//int t = (int)this->getContext()->getTime();
		cv::Rect ROI(160, 120, 320, 240);

		//if (t%niterations.getValue() == 0)
		{	
	raDepth rdepth(this->depthImage);
	if( rdepth->isEmpty() || !mstate)  return;

	const CImg<dT>& depthimg =rdepth->getCImg(0); 
    cv::Mat depth0,depthuc; 	
	int width, height;		
	switch (sensorType.getValue())
    {
    // FLOAT ONE CHANNEL
    case 0:
	height = 240;
	width = 320;
	break;
	case 1:
	height = 480;
	width = 640;
	break;
	}
	cv::Mat depth_single = cv::Mat::zeros(height,width,CV_32FC1); 
	memcpy(depth_single.data, (float*)depthimg.data(), height*width*sizeof(float));
		   
		   	switch (sensorType.getValue())
			{
		// FLOAT ONE CHANNEL
		case 0:
		depth = depth_single;
		break;
		case 1:
		depth = depth_single(ROI);
		break;
			}
			
		cv::namedWindow("depth_softkinetic");
		cv::imshow("depth_softkinetic",depth);
		

	raImage rimg(this->image);
	if( rimg->isEmpty() || !mstate)  return;
	const CImg<T>& img =rimg->getCImg(0);

	cv::Mat color0; 
	color0 = cv::Mat::zeros(img.height(),img.width(), CV_8UC3); 
	if(img.spectrum()==3)  // deinterlace
		            {
                        unsigned char* rgb = (unsigned char*)color0.data;
                        const unsigned char *ptr_r = img.data(0,0,0,2), *ptr_g = img.data(0,0,0,1), *ptr_b = img.data(0,0,0,0);
                        for ( int siz = 0 ; siz<img.width()*img.height(); siz++)    {*(rgb++) = *(ptr_r++) ; *(rgb++) = *(ptr_g++); *(rgb++) = *(ptr_b++); }
                    }
	
	        color_1 = color;
			//color_2 = color; 
	
		    switch (sensorType.getValue())
			{
		case 0:		
		cv::pyrDown(color0, color, cv::Size(color0.cols/2, color0.rows/2));
		break;
		case 1:
		color = color0(ROI);
		break;
			}
					std::cout << " k init 0 " << std::endl;

	cv::namedWindow("image_softkinetic");
    cv::imshow("image_softkinetic",color);
							//cv::imwrite("color01.png", color);
		}
}

template <class DataTypes, class DepthTypes>
void ImageConverter<DataTypes, DepthTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{

	int t = (int)this->getContext()->getTime();
			
	//if (useRealData.getValue())
	//if (useSensor.getValue())
	 getImages();
		
}


template<class DataTypes, class DepthTypes>
void ImageConverter<DataTypes, DepthTypes>::draw(const core::visual::VisualParams* vparams)
{	


}

}
}
} // namespace sofa



