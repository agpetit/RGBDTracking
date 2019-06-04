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
*                               SOFA :: Plugins                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define SOFA_RGBDTRACKING_VIRTUALCAMERA_CPP


#include "VirtualCamera.inl"

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/Mapping.inl>
#include <sofa/simulation/Simulation.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/gui/BaseGUI.h>
#include <sofa/gui/BaseViewer.h>
#include <sofa/gui/GUIManager.h>
#include <iostream>
#include <map>

#ifdef USING_OMP_PRAGMAS
    #include <omp.h>
#endif


using std::cerr;
using std::endl;

namespace sofa {

namespace rgbdtracking {


    using namespace sofa::defaulttype;

      SOFA_DECL_CLASS(VirtualCamera)

      // Register in the Factory
      int VirtualCameraClass = core::RegisterObject("Defines a virtual camera component")
    #ifndef SOFA_FLOAT
        .add< VirtualCamera<Vec3dTypes> >()
    #endif
    #ifndef SOFA_DOUBLE
        .add< VirtualCamera<Vec3fTypes> >()
    #endif
    ;

    #ifndef SOFA_FLOAT
      template class SOFA_RGBDTRACKING_API VirtualCamera<Vec3dTypes>;
    #endif
    #ifndef SOFA_DOUBLE
      template class SOFA_RGBDTRACKING_API VirtualCamera<Vec3fTypes>;

    #endif

using namespace helper;

}
} // namespace sofa

