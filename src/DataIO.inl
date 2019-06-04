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

#define SOFA_RGBDTRACKING_DATAIO_INL
#include <limits>
#include <iterator>
#include <sofa/core/behavior/ForceField.inl>

#include <pcl/common/common_headers.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/registration/icp.h>
#include <pcl/common/transforms.h>

#include <sofa/helper/AdvancedTimer.h>

#include "DataIO.h"

using std::cerr;
using std::endl;

namespace sofa
{

namespace rgbdtracking {
	
using namespace sofa::defaulttype;
using namespace helper;


template <class DataTypes>
DataIO<DataTypes>::DataIO()
    : Inherit()
    , inputPath(initData(&inputPath,"inputPath","Path for data readings",false))
    , outputPath(initData(&outputPath,"outputPath","Path for data writings",false))
    , nimages(initData(&nimages,"nimages","Number of images",false))
    , startimage(initData(&startimage,1,"startimage","Number of images"))
    , niterations(initData(&niterations,1,"niterations","Number of images"))
{

    std::cout << " init data " << std::endl;
    this->f_listening.setValue(true);

}

template <class DataTypes>
DataIO<DataTypes>::~DataIO()
{
}

template<class DataTypes>
void DataIO<DataTypes>::readImages()
{
    int t = (int)this->getContext()->getTime();
    if (t%(npasses) != 0 ) {
        return ;
    } else if (t % niterations.getValue() != 0) {
        newImages.setValue(false) ;
    } else {
        newImages.setValue(true);
        RGBDFileSystemIO mydumper (inputPath.getValue(), iter_im) ;
        iter_im++;

        depth = mydumper.read_depth_file() ;
        color = mydumper.read_image() ;

        color_1 = color.clone();
    }
}


template <class DataTypes>
void DataIO<DataTypes>::init()
{
    iter_im = startimage.getValue();
    //iter_im = 1; //crocodile disk
   //iter_im = 20; // cube disk liver
    //iter_im = 440; //pig 2018
   //iter_im = 1;
   //iter_im = 300; //patient1liver1
   //iter_im = 50; //patient1liver2
   //iter_im = 1; //patient2liver1
   //iter_im = 370; //patient2liver2
   listimg.resize(0);
   listimgseg.resize(0);
   listimgklt.resize(0);
   listdepth.resize(0);
   listrtt.resize(0);

   listrttstress.resize(0);
   listrttstressplast.resize(0);
//   listpcd.resize(0);
//   listvisible.resize(0);
   pcl = false;
   disp = false;
   npasses = 1;
}

template <class DataTypes>
void DataIO<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event) {
    if (dynamic_cast<simulation::AnimateBeginEvent*>(event)) {
        sofa::helper::AdvancedTimer::stepBegin("ImageReading") ;
        readImages();
        sofa::helper::AdvancedTimer::stepEnd("ImageReading") ;

        sofa::helper::AdvancedTimer::stepBegin("ImageWriting") ;
        if ((int)this->getContext()->getTime() == nimages.getValue()*niterations.getValue()) {
            writeImages();
        }
        sofa::helper::AdvancedTimer::stepEnd("ImageWriting") ;
    }
}

template <class DataTypes>
void DataIO<DataTypes>::writeImages() {

    for (int frame_count = 0 ;frame_count < listimgseg.size()-5; frame_count++) {
        std::cout << "(DataIO) write image : " << frame_count << std::endl;

        RGBDFileSystemIO mydumper (outputPath.getValue(), frame_count) ;

        cv::Mat img = *listimg[frame_count];
        mydumper.write_image(img);

        cv::Mat imgseg = *listimgseg[frame_count],
                imgseg1 ;
        cvtColor(imgseg,imgseg1 ,CV_RGBA2RGB);
        mydumper.write_segimage(imgseg);

        cv::Mat deptht = *listdepth[frame_count],
                deptht1 ;
        deptht.convertTo (deptht1, CV_8UC1, 100);
        mydumper.write_depth_bis(deptht);
        mydumper.write_depth(deptht1);

        mydumper.write_depth_file(deptht) ;

//        if (useKLTPoints.getValue()){
//            cv::Mat imgklt = *listimgklt[frame_count];
//            mydumper.write_klt(imgklt);
//        }

        delete listimg[frame_count];
        delete listimgseg[frame_count];
        delete listdepth[frame_count];
    }


    for (int frame_count = 1 ;frame_count < listrtt.size(); frame_count++) {
        std::cout << "(DataIO) write rtt : " << frame_count << std::endl;

        RGBDFileSystemIO mydumper (outputPath.getValue(), frame_count-1) ;

        cv::Mat rtt,rtt1, rttstress, rttstressplast;
        rtt = *listrtt[frame_count];
        cvtColor(rtt,rtt1 ,CV_RGB2BGR);
        mydumper.write_rtt(rtt1);

        if (npasses == 2) {
            rtt = *listrttstress[frame_count];
            cvtColor(rtt,rttstress ,CV_RGB2BGR);
            mydumper.write_rtt_stress(rttstress);
        }
        // the fuck is that ??
        if (npasses == 2){
            rtt = *listrttstressplast[frame_count];
            cvtColor(rtt,rttstressplast ,CV_RGB2BGR);
            mydumper.write_rtt_stress_plast(rttstressplast);
        }
        delete listrtt[frame_count];
    }
}

} // rgbdtracking

} // namespace sofa

