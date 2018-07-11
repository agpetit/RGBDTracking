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

#define SOFA_RGBDTRACKING_VIRTUALCAMERA_INL

#include <limits>
#include <iterator>

#include <sofa/helper/gl/Color.h>
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseVisual/BaseCamera.h>
#include <SofaBaseVisual/InteractiveCamera.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/simulation/Simulation.h>

#include <algorithm>

#include "VirtualCamera.h"

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
VirtualCamera<DataTypes>::VirtualCamera( )
 : Inherit()
 	, cameraIntrinsicParameters(initData(&cameraIntrinsicParameters,Vector4(),"cameraIntrinsicParameters","camera parameters"))
        , cameraPosition(initData(&cameraPosition,"cameraPosition","Position of the camera w.r.t the point cloud"))
        , cameraOrientation(initData(&cameraOrientation,"cameraOrientation","Orientation of the camera w.r.t the point cloud"))
        , viewportWidth(initData(&viewportWidth,640,"viewportWidth","Width of the viewport"))
        , viewportHeight(initData(&viewportHeight,480,"viewportHeight","Height of the viewport"))
        , cameraChanged(initData(&cameraChanged,false,"cameraChanged","If the camera has changed or not"))
{
	this->f_listening.setValue(true); 
}

template <class DataTypes>
VirtualCamera<DataTypes>::~VirtualCamera()
{
}

template <class DataTypes>
void VirtualCamera<DataTypes>::init()
{		
}

template<class DataTypes>
void VirtualCamera<DataTypes>::setCamera()
{
    Vector4 camParam = cameraIntrinsicParameters.getValue();
    rgbIntrinsicMatrix(0,0) = camParam[0];
    rgbIntrinsicMatrix(1,1) = camParam[1];
    rgbIntrinsicMatrix(0,2) = camParam[2];
    rgbIntrinsicMatrix(1,2) = camParam[3];
    int hght = viewportHeight.getValue();
    int wdth = viewportWidth.getValue();
    sofa::gui::GUIManager::SetDimension(wdth,hght);
    sofa::component::visualmodel::BaseCamera::SPtr currentCamera;
    sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
    root->get(currentCamera);
    currentCamera->p_fieldOfView.setValue(atan((hght * 0.5 ) / rgbIntrinsicMatrix(1,1)) * 360.0 / M_PI);
    sofa::gui::BaseGUI *gui = sofa::gui::GUIManager::getGUI();
    sofa::gui::BaseViewer * viewer = gui->getViewer();
    viewer->setView(cameraPosition.getValue(),cameraOrientation.getValue());
}

template <class DataTypes>
void VirtualCamera<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{
        if (dynamic_cast<simulation::AnimateBeginEvent*>(event))
	{	
       double timeT = (double)getTickCount();
       if (cameraChanged.getValue())
           setCamera();
       std::cout << "TIME VIRTUAL CAMERA " << ((double)getTickCount() - timeT)/getTickFrequency() << std::endl;
        }
}

template <class DataTypes>
void VirtualCamera<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

}



}
}
} // namespace sofa

