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

#define SOFA_RGBDTRACKING_MESHPROCESSING_INL

#include <limits>
#include <iterator>
#include <sofa/helper/gl/Color.h>
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseVisual/InteractiveCamera.h>
#include <sofa/core/behavior/ForceField.inl>
#include  <sofa/simulation/Simulation.h>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/Mapping.inl>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/gui/BaseGUI.h>
#include <sofa/gui/BaseViewer.h>
#include <sofa/gui/GUIManager.h>
#include <sofa/simulation/InitVisitor.h>

#ifdef USING_OMP_PRAGMAS
#include <omp.h>
#endif

#include <SofaLoader/MeshObjLoader.h>
#include <limits>
#include <iterator>
#include <sofa/helper/gl/Color.h>


#ifdef Success
  #undef Success
#endif

#include <pcl/common/common_headers.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#include <algorithm>    // std::max
#include "MeshProcessing.h"


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
MeshProcessing<DataTypes>::MeshProcessing( )
 : Inherit()
 	, cameraIntrinsicParameters(initData(&cameraIntrinsicParameters,Vector4(),"cameraIntrinsicParameters","camera parameters"))
	, sourceSurfacePositions(initData(&sourceSurfacePositions,"sourceSurface","Points of the surface of the source mesh."))
	, sourcePositions(initData(&sourcePositions,"sourcePositions","Points of the mesh."))
        , sourceVisiblePositions(initData(&sourceVisiblePositions,"sourceVisiblePositions","Visible points of the surface of the mesh."))
        , sourceVisible(initData(&sourceVisible,"sourceVisible","Visibility of the points of the surface of the mesh."))
        , sourceBorder(initData(&sourceBorder,"sourceBorder","Points of the border of the mesh."))
        , indicesVisible(initData(&indicesVisible,"indicesVisible","Indices of the visible points of the mesh."))
        , sourceContourPositions(initData(&sourceContourPositions,"sourceContourPositions","Contour points of the surface of the mesh."))
        , sourceTriangles(initData(&sourceTriangles,"sourceTriangles","Triangles of the source mesh."))
        , sourceNormals(initData(&sourceNormals,"sourceNormals","Normals of the source mesh."))
	, sourceSurfaceNormals(initData(&sourceSurfaceNormals,"sourceSurfaceNormals","Normals of the surface of the source mesh."))
	, useContour(initData(&useContour,false,"useContour","Emphasize forces close to the target contours"))
	, useVisible(initData(&useVisible,true,"useVisible","Use the vertices of the viisible surface of the source mesh"))
	, visibilityThreshold(initData(&visibilityThreshold,(Real)0.001,"visibilityThreshold","Threshold to determine visible vertices"))
	, niterations(initData(&niterations,3,"niterations","Number of iterations in the tracking process"))
	, borderThdSource(initData(&borderThdSource,7,"borderThdSource","border threshold on the source silhouette"))
{
	
	this->f_listening.setValue(true); 
	iter_im = 0;

        hght = 0;
        wdth = 0;
	
	rectRtt.x = 0;
	rectRtt.y = 0;
	rectRtt.height = hght;
	rectRtt.width = wdth;

        timeMeshProcessing = 0;

}

template <class DataTypes>
MeshProcessing<DataTypes>::~MeshProcessing()
{
}

template <class DataTypes>
void MeshProcessing<DataTypes>::init()
{

    this->Inherit::init();
    core::objectmodel::BaseContext* context = this->getContext();
    mstate = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());
    sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
    root->get(renderingmanager);

    Vector4 camParam = cameraIntrinsicParameters.getValue();

    rgbIntrinsicMatrix(0,0) = camParam[0];
    rgbIntrinsicMatrix(1,1) = camParam[1];
    rgbIntrinsicMatrix(0,2) = camParam[2];
    rgbIntrinsicMatrix(1,2) = camParam[3];

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

}

template<class DataTypes>
void MeshProcessing<DataTypes>::getSourceVisible(double znear, double zfar)
{
	
    int t = (int)this->getContext()->getTime();
    //if (t%2 == 0)
    {
        cv::Mat _rtd0, depthr;
        renderingmanager->getDepths(depthr);
        depthrend = depthr.clone();

        wdth = depthr.cols;
        hght = depthr.rows;

        _rtd0.create(hght,wdth, CV_8UC1);
        double depthsN[hght * wdth ];

	std::vector<cv::Point> ptfgd;
	ptfgd.resize(0);
        cv::Point pt;
        double clip_z;
        //double timef0 = (double)getTickCount();

        for (int j = 0; j < wdth; j++)
            for (int i = 0; i< hght; i++)
            {
                if ((double)depthr.at<float>(hght-i-1,j) < 1  && (double)depthr.at<float>(hght-i-1,j)	> 0.001)
                {
                    //if (j >= rectRtt.x && j < rectRtt.x + rectRtt.width && i >= rectRtt.y && i < rectRtt.y + rectRtt.height) {
                    //if ((double)(float)depths1[j-rectRtt.x+(i-rectRtt.y)*(rectRtt.width)]	< 1){
                    _rtd0.at<uchar>(hght-i-1,j) = 255;

                    clip_z = (depthr.at<float>(hght-i-1,j) - 0.5) * 2.0;
                    //double clip_z = (depths1[j-rectRtt.x+(i-rectRtt.y)*(rectRtt.width)] - 0.5) * 2.0;
                    depthsN[j+i*wdth] = 2*znear*zfar/(clip_z*(zfar-znear)-(zfar+znear));

                    pt.x = j;
                    pt.y = i;
                    ptfgd.push_back(pt);
                }
                else
                {
                    _rtd0.at<uchar>(hght-i-1,j) = 0;
                    depthsN[j+i*wdth] = 0;
                }

            }
		
        if(t > 0)
        {
            rectRtt = cv::boundingRect(ptfgd);
            //std::cout << " ptfgd " << frame_count << " rect1 " << rectangle.x << " " << rectangle.y << " rect2 " << rectangle.width << " " << rectangle.height << std::endl;
            rectRtt.x -= 10;
            rectRtt.y -= 10;
            rectRtt.height += 20;
            rectRtt.width += 20;
        }
		
        depthMap = _rtd0.clone();
        //depthr.convertTo(depthu, CV_8UC1, 100);
        //cv::imwrite("depth001.png",depthu);
        //cv::imwrite("depth01.png", depthMap);
        const VecCoord& x = mstate->read(core::ConstVecCoordId::position())->getValue();

        helper::vector<bool> sourcevisible;
        sourcevisible.resize(x.size());
        VecCoord sourceVis;
        Vector3 pos;

        helper::vector< int > indicesvisible;
        indicesvisible.resize(0);


        for (unsigned int k = 0; k < x.size(); k++)
        {
            int x_u = (int)(x[k][0]*rgbIntrinsicMatrix(0,0)/x[k][2] + rgbIntrinsicMatrix(0,2));
            int x_v = (int)(x[k][1]*rgbIntrinsicMatrix(1,1)/x[k][2] + rgbIntrinsicMatrix(1,2));
           // std::cout << " depths00 " << (float)x_u << " " << (double)x_v << std::endl;
            //std::cout << " depths01 " << (float)depthsN[x_u+(hght-x_v-1)*wdth] << " " <<(float)x[k][2] << " " << rgbIntrinsicMatrix(1,1) << " " <<  rgbIntrinsicMatrix(1,2) << std::endl;
	
            if (x_u>=0 && x_u<wdth && x_v<hght && x_v >= 0){
                if((float)abs(depthsN[x_u+(hght-x_v-1)*wdth]+(float)x[k][2]) < visibilityThreshold.getValue() || (float)depthsN[x_u+(hght-x_v-1)*wdth] == 0)
                {
                    sourcevisible[k] = true;
                    pos = x[k];
                    sourceVis.push_back(pos);
                    indicesvisible.push_back(k);
                }
                else
                {
                    sourcevisible[k] = false;
                }
            }
            else {sourcevisible[k] = false;}
	
        }

        std::cout << " nvisible " << sourceVis.size() << " xsize " << sourcevisible.size() <<  std::endl;
        sourceVisiblePositions.setValue(sourceVis);
        sourceVisible.setValue(sourcevisible);
        indicesVisible.setValue(indicesvisible);

    }
}

template<class DataTypes>
void MeshProcessing<DataTypes>::updateSourceVisible()
{
    const VecCoord&  x = mstate->read(core::ConstVecCoordId::position())->getValue();
    VecCoord sourceVis;
    helper::vector<bool> sourcevisible = sourceVisible.getValue();
    Vector3 pos;			
        for (unsigned int i=0; i< x.size(); i++)
	{
            if (sourcevisible[i])
            {
                pos = x[i];
                sourceVis.push_back(pos);
            }
	}
    sourceVisiblePositions.setValue(sourceVis);
	
}

template<class DataTypes>
void MeshProcessing<DataTypes>::extractSourceContour()
{
		
    double cannyTh1 = 350;
    double cannyTh2 = 10;
    cv::Mat contour,dist,dist0;
    //cv::imwrite("depthmap.png", depthMap);
    cv::Canny( depthMap, contour, cannyTh1, cannyTh2, 3);
    contour = cv::Scalar::all(255) - contour;

    cv::distanceTransform(contour, dist, CV_DIST_L2, 3);
    dist.convertTo(dist0, CV_8U, 1, 0);
    pcl::PointCloud<pcl::PointXYZRGB> sourceContour;

    const VecCoord& x = mstate->read(core::ConstVecCoordId::position())->getValue();
    unsigned int nbs=x.size();
    cv::Mat contourpoints = cv::Mat::zeros(hght,wdth, CV_8UC3);
    contourpoints= cv::Mat(hght,wdth,CV_8UC3,cv::Scalar(255,255,255));
        for (int j = 0; j < wdth; j++)
            for (int i = 0; i< hght; i++)
            {
		contourpoints.at<Vec3b>(i,j)[0]= contour.at<uchar>(i,j);
            }
	
    pcl::PointXYZRGB newPoint;
    helper::vector< bool > sourceborder;
    sourceborder.resize(nbs);
    //cv::imwrite("dist0.png", dist0);
	
    int nsourcecontour = 0;
    sourceWeights.resize(0);
	
    double totalweights = 0;

        for (unsigned int i=0; i<nbs; i++)
	{
            int x_u = (int)(x[i][0]*rgbIntrinsicMatrix(0,0)/x[i][2] + rgbIntrinsicMatrix(0,2));
            int x_v = (int)(x[i][1]*rgbIntrinsicMatrix(1,1)/x[i][2] + rgbIntrinsicMatrix(1,2));
            int thickness = 1;
            int lineType = 2;

            if (dist0.at<uchar>(x_v,x_u) < borderThdSource.getValue()/*7*/)
            {
                newPoint.z = x[i][2];
                newPoint.x = x[i][0];
                newPoint.y = x[i][1];
                newPoint.r = 0;
                newPoint.g = 0;
                newPoint.b = 0;
                sourceContour.push_back(newPoint);
                sourceborder[i] = true;
                nsourcecontour++;
                circle( contourpoints,cv::Point(x_u,x_v),wdth/128.0,Scalar( 0, 0, 255 ),thickness,lineType );
            }
            else sourceborder[i] = false;

            //sourceWeights.push_back((double)1./(0.12*(1.0+sqrt(dist0.at<uchar>(x_v,x_u)))));
            sourceWeights.push_back((double)exp(-dist0.at<uchar>(x_v,x_u)/sigmaWeight.getValue()));
            totalweights += sourceWeights[i];
	}
	
        for (unsigned int i=0; i < sourceWeights.size();i++)
	{
            sourceWeights[i]*=((double)sourceWeights.size()/totalweights);
	}
	
    //cv::imwrite("contourpoints.png", contourpoints);
    //cv::imwrite("dist.png", dist);
	
    //std::cout << " n source " << nbs << " n source contour " << nsourcecontour << std::endl;
	
    VecCoord sourcecontourpos;
    sourcecontourpos.resize(sourceContour.size());
    Vector3 pos;

        for (unsigned int i=0; i<sourcecontourpos.size(); i++)
        {
            pos[0] = sourceContour[i].x;
            pos[1] = sourceContour[i].y;
            pos[2] = sourceContour[i].z;
            sourcecontourpos[i]=pos;
        }
    const VecCoord&  p = sourcecontourpos;
    sourceContourPositions.setValue(p);
    sourceBorder.setValue(sourceborder);
	
}

template<class DataTypes>
void MeshProcessing<DataTypes>::extractSourceVisibleContour()
{
	
    double cannyTh1 = 350;
    double cannyTh2 = 10;
    cv::Mat contour,dist,dist0;
    //cv::imwrite("depthmap.png", depthMap);
	
    cv::Canny( depthMap, contour, cannyTh1, cannyTh2, 3);
    contour = cv::Scalar::all(255) - contour;
    cv::distanceTransform(contour, dist, CV_DIST_L2, 3);
    dist.convertTo(dist0, CV_8U, 1, 0);
    pcl::PointCloud<pcl::PointXYZRGB> sourceContour;
	
    const VecCoord& x =  mstate->read(core::ConstVecCoordId::position())->getValue();
	
    unsigned int nbs=x.size();
    cv::Mat contourpoints = cv::Mat::zeros(hght,wdth, CV_8UC3);
    contourpoints= cv::Mat(hght,wdth,CV_8UC3,cv::Scalar(255,255,255));
	for (int j = 0; j < wdth; j++)
	  for (int i = 0; i< hght; i++)
          {
            contourpoints.at<Vec3b>(i,j)[0]= contour.at<uchar>(i,j);
          }
	
    pcl::PointXYZRGB newPoint;
    helper::vector< bool > sourceborder;
    sourceborder.resize(nbs);
    //cv::imwrite("dist0.png", dist0);
    int nsourcecontour = 0;
    sourceWeights.resize(0);
    double totalweights = 0;
        for (unsigned int i=0; i<nbs; i++)
	{
            int x_u = (int)(x[i][0]*rgbIntrinsicMatrix(0,0)/x[i][2] + rgbIntrinsicMatrix(0,2));
            int x_v = (int)(x[i][1]*rgbIntrinsicMatrix(1,1)/x[i][2] + rgbIntrinsicMatrix(1,2));
            int thickness = 1;
            int lineType = 2;

            if (dist0.at<uchar>(x_v,x_u) < borderThdSource.getValue() /*6*/ && (sourceVisible.getValue())[i])
            {
                newPoint.z = x[i][2];
                newPoint.x = x[i][0];
                newPoint.y = x[i][1];
                newPoint.r = 0;
                newPoint.g = 0;
                newPoint.b = 0;
                sourceContour.push_back(newPoint);
                sourceborder[i] = true;
                nsourcecontour++;
                circle( contourpoints,cv::Point(x_u,x_v),wdth/128.0,Scalar( 0, 0, 255 ),thickness,lineType );
            }
            else sourceborder[i] = false;
            //sourceWeights.push_back((double)1./(0.12*(1.0+sqrt(dist0.at<uchar>(x_v,x_u)))))
            sourceWeights.push_back((double)exp(-dist0.at<uchar>(x_v,x_u)/sigmaWeight.getValue()));
            totalweights += sourceWeights[i];
	}
	
        for (unsigned int i=0; i < sourceWeights.size();i++)
	{
            sourceWeights[i]*=((double)sourceWeights.size()/totalweights);
	}
	
    //cv::imwrite("contourpoints.png", contourpoints);
    //cv::imwrite("dist.png", dist);
    //std::cout << " n source " << nbs << " n source contour " << nsourcecontour << " " << borderThdSource.getValue() << std::endl;
	
    VecCoord sourcecontourpos;
    sourcecontourpos.resize(sourceContour.size());
    Vector3 pos;

        for (unsigned int i=0; i<sourcecontourpos.size(); i++)
        {
            pos[0] = sourceContour[i].x;
            pos[1] = sourceContour[i].y;
            pos[2] = sourceContour[i].z;
            sourcecontourpos[i]=pos;
        }
    const VecCoord&  p = sourcecontourpos;
    sourceContourPositions.setValue(p);
    sourceBorder.setValue(sourceborder);
	
}

template<class DataTypes>
void MeshProcessing<DataTypes>::updateSourceVisibleContour()
{
    // build k-d tree
    const VecCoord&  x = mstate->read(core::ConstVecCoordId::position())->getValue();
    VecCoord sourcecontourpos;
    sourcecontourpos.resize(sourceContourPositions.getValue().size());
    Vector3 pos;
    int k = 0;

    helper::vector< bool > sourceborder = sourceBorder.getValue();

        for (unsigned int i=0; i<x.size(); i++)
	{
            if (sourceborder[i])
            {
                pos = x[i];
                sourcecontourpos[k]=pos;
                k++;
	    }
	}
    const VecCoord&  p = sourcecontourpos;
    sourceContourPositions.setValue(p);
	
}

template <class DataTypes>
void MeshProcessing<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{
    if (dynamic_cast<simulation::AnimateBeginEvent*>(event))
    {
        int t = (int)this->getContext()->getTime();
        timeMeshProcessing = (double)getTickCount();
            if (!depthrend.empty() && t%niterations.getValue() == 0)
            {
		if (!useVisible.getValue())
		{
                    if(useContour.getValue())
                    extractSourceContour();
		}				
		else 
		{
                    //if (t%(npasses + niterations.getValue() - 1) ==0 )
                    {
                        double znear = renderingmanager->getZNear();
                        double zfar = renderingmanager->getZFar();
                        //std::cout << " znear01 " << znear << " zfar01 " << zfar << std::endl;
                        getSourceVisible(znear, zfar);
                    }

                    if(useContour.getValue())
                        extractSourceVisibleContour();
                }
	}
	
            if (!depthrend.empty() && useVisible.getValue() && t%niterations.getValue()!= 0)
            {
                if(useContour.getValue())
                    updateSourceVisibleContour();
                else updateSourceVisible();

            }
            timeMeshProcessing = ((double)getTickCount() - timeMeshProcessing)/getTickFrequency();
            cout << "TIME MESHPROCESSING " << timeMeshProcessing << endl;
		
    }
}

template <class DataTypes>
void MeshProcessing<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

}


}
}
} // namespace sofa

