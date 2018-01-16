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

#define SOFA_RGBDTRACKING_RENDERTEXTUREAR_INL

#include <SofaGeneralEngine/NormalsFromPoints.h>
#include <limits>
#include <set>
#include <iterator>
#include <sofa/helper/gl/Color.h>
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseVisual/InteractiveCamera.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/simulation/Simulation.h>

#include <pcl/common/common_headers.h>
#include <pcl/point_cloud.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/point_types.h>
#include <pcl/impl/point_types.hpp>

#include <algorithm>    // std::max
#include "RenderTextureAR.h"


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

template <class DataTypes>
RenderTextureAR<DataTypes>::RenderTextureAR()
 : Inherit(){

}


template <class DataTypes>
RenderTextureAR<DataTypes>::~RenderTextureAR()
{
}

template <class DataTypes>
void RenderTextureAR<DataTypes>::init()
{
    sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
    root->get(renderingmanager);
}

template<class DataTypes>
void RenderTextureAR<DataTypes>::renderToTexture(cv::Mat &_rtt)
{
    renderingmanager->getTexture(_rtt);
}

template<class DataTypes>
void RenderTextureAR<DataTypes>::renderToTextureD(cv::Mat &_rtt, cv::Mat &color_1)
{
    renderingmanager->getTexture(_rtt);
    int hght = _rtt.rows;
    int wdth = _rtt.cols;
    cv::Mat _rttd,_dd,_rttd_;
    _rttd.create(hght,wdth, CV_8UC3);
    //cv::imwrite("rtt.png",_rtt);

    //_rttd = _rtt;

    cv::resize(_rtt, _rttd_,_rttd.size());
    cv::flip( _rttd_,_rttd,0);

    cv::Mat _rtd, _rtd0, depthmat;
    _rtd.create(hght, wdth, CV_32F);
    _rtd0.create(hght,wdth, CV_8UC1);

    renderingmanager->getTexture(depthmat);

        for (int j = 0; j < wdth; j++)
                for (int i = 0; i< hght; i++)
                {
                        if ((double)depthmat.at<uchar>(i,j)< 1)
                        {
                        _rtd0.at<uchar>(i,j) = 255;
                        }
                        else _rtd0.at<uchar>(i,j) = 0;
                }
        cv::Mat col1 = color_1.clone();

        for (int j = 0; j < wdth; j++)
                for (int i = 0; i< hght; i++)
                {
                        if (color_1.at<Vec3b>(i,j)[0] == 0 && color_1.at<Vec3b>(i,j)[1] == 0 && color_1.at<Vec3b>(i,j)[2] == 0)
                        {
                                /*col1.at<Vec3b>(i,j)[0] = 255;
                                col1.at<Vec3b>(i,j)[2] = 255;
                                col1.at<Vec3b>(i,j)[1] = 255;*/
                                /*_rttd.at<Vec3b>(i,j)[0] = 255;
                                _rttd.at<Vec3b>(i,j)[2] = 255;
                                _rttd.at<Vec3b>(i,j)[1] = 255;*/
                        }


                        if (/*_rtd0.at<uchar>(i,j) == 0 ||*/ (_rttd.at<Vec3b>(i,j)[0] > 0 && _rttd.at<Vec3b>(i,j)[1] > 0 && _rttd.at<Vec3b>(i,j)[2] > 0))
                        {
                                _rttd.at<Vec3b>(i,j)[0] = col1.at<Vec3b>(i,j)[2];
                                _rttd.at<Vec3b>(i,j)[2] = col1.at<Vec3b>(i,j)[0];
                                _rttd.at<Vec3b>(i,j)[1] = col1.at<Vec3b>(i,j)[1];
                        }



                }

        cv::namedWindow("dd");
        cv::imshow("dd",_rttd);
        _rtt = _rttd.clone();
        cv::waitKey(1);
}


template<class DataTypes>
void RenderTextureAR<DataTypes>::renderToTextureDepth(cv::Mat &_rtt, cv::Mat &_rttdepth)
{
    renderingmanager->getTexture(_rtt);
    int hght = _rtt.rows;
    int wdth = _rtt.cols;

    cv::Mat _dd;
	
    _rttdepth.create(hght,wdth,CV_32F);
    _dd.create(hght,wdth,CV_32F);
	
    renderingmanager->getDepths(_rttdepth);

    double znear = renderingmanager->getZNear();
    double zfar = renderingmanager->getZFar();
	
	for (int j = 0; j < wdth; j++)
		for (int i = 0; i< hght; i++)
		{
                        if (_rttdepth.at<float>(hght-i-1,j)<1)
			{
                                double clip_z = (_rttdepth.at<float>(hght-i-1,j) - 0.5) * 2.0;
                        _dd.at<float>(i,j) = -2*znear*zfar/(clip_z*(zfar-znear)-(zfar+znear));
			//double clip_z = (depths1[j-rectRtt.x+(i-rectRtt.y)*(rectRtt.width)] - 0.5) * 2.0;
                        //std::cout << " depth " << znear << " " << zfar << " " << -2*znear*zfar/(clip_z*(zfar-znear)-(zfar+znear)) << std::endl;

			}


		}
        _rttdepth = _dd.clone();
}


}
}
} // namespace sofa

