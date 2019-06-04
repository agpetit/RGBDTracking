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

#define SOFA_RGBDTRACKING_IMAGECONVERTER_CPP

#include "ImageConverter.inl"

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/objectmodel/BaseContext.h>

using std::cerr;
using std::endl;

namespace sofa
{

namespace rgbdtracking {

    using namespace sofa::defaulttype;

      SOFA_DECL_CLASS(ImageConverter)

      // Register in the Factory
      int ImageConverterClass = core::RegisterObject("converts input data each niter steps") //rescaling mainly ?
    #ifndef SOFA_FLOAT
        .add< ImageConverter<Vec3dTypes,ImageUC> >()
        .add< ImageConverter<Vec3dTypes,ImageUS> >()
        .add< ImageConverter<Vec3dTypes,ImageF> >()
    #endif
    #ifndef SOFA_DOUBLE
        .add< ImageConverter<Vec3fTypes,ImageUC> >()
        .add< ImageConverter<Vec3fTypes,ImageUS> >()
        .add< ImageConverter<Vec3fTypes,ImageF> >()
    #endif
    ;

    #ifndef SOFA_FLOAT
      template class SOFA_RGBDTRACKING_API ImageConverter<Vec3dTypes,ImageUC>;
      template class SOFA_RGBDTRACKING_API ImageConverter<Vec3dTypes,ImageUS>;
      template class SOFA_RGBDTRACKING_API ImageConverter<Vec3dTypes,ImageF>;
    #endif
    #ifndef SOFA_DOUBLE
      template class SOFA_RGBDTRACKING_API ImageConverter<Vec3fTypes,ImageUC>;
      template class SOFA_RGBDTRACKING_API ImageConverter<Vec3fTypes,ImageUS>;
      template class SOFA_RGBDTRACKING_API ImageConverter<Vec3fTypes,ImageF>;
    #endif


}

} // namespace sofa



