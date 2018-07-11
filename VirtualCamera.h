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

#ifndef SOFA_RGBDTRACKING_VIRTUALCAMERA_H
#define SOFA_RGBDTRACKING_VIRTUALCAMERA_H

#include <RGBDTracking/config.h>
#include <boost/thread.hpp>

#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <SofaBaseTopology/TopologyData.h>
#include <sofa/gui/BaseGUI.h>
#include <sofa/gui/BaseViewer.h>
#include <sofa/gui/GUIManager.h>

#include <set>

#include <pthread.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>

#include <opencv/cv.h>
#include <opencv2/core.hpp>

#include <sys/times.h>

using namespace std;
using namespace cv;


namespace sofa
{

namespace core
{

namespace objectmodel
{

using helper::vector;
using namespace sofa::defaulttype;
using namespace sofa::component::topology;


template<class DataTypes>
class VirtualCameraInternalData
{
public:
};

template<class DataTypes>
class VirtualCamera : public sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(VirtualCamera,DataTypes),sofa::core::objectmodel::BaseObject);
	
    typedef sofa::core::objectmodel::BaseObject Inherit;

    typedef sofa::defaulttype::Vector4 Vector4;
    typedef sofa::defaulttype::Vector3 Vec3;

    Data<Vector4> cameraIntrinsicParameters;
    Eigen::Matrix3f rgbIntrinsicMatrix;

    Data<Vec3> cameraPosition;
    Data<Quat> cameraOrientation;

    Data<int> viewportHeight;
    Data<int> viewportWidth;

    Data<bool> cameraChanged;
    double timePCD;
	
    VirtualCamera();
    virtual ~VirtualCamera();

    void init();
    void handleEvent(sofa::core::objectmodel::Event *event);
    void setCamera();
    void draw(const core::visual::VisualParams* vparams) ;




};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(VirtualCamera_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_RGBDTRACKING_API VirtualCamera<defaulttype::Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_RGBDTRACKING_API VirtualCamera<defaulttype::Vec3fTypes>;
#endif
#endif


} //

} //

} // namespace sofa

#endif
