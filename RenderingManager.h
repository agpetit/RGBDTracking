/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2017 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#ifndef SOFA_RGBDTRACKING_RENDERINGMANAGER_H_
#define SOFA_RGBDTRACKING_RENDERINGMANAGER_H_
#include <RGBDTracking/config.h>
#include "DataIO.h"


#include <sofa/core/visual/VisualManager.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/gl/FrameBufferObject.h>
#include <SofaOpenglVisual/OglShader.h>
#include <sofa/core/objectmodel/DataFileName.h>

#include <opencv/cv.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>

#include <iostream>


using namespace std;
using namespace cv;

namespace sofa
{

namespace component
{

namespace visualmodel
{

class RenderingManager : public core::visual::VisualManager
{
public:
    SOFA_CLASS(RenderingManager,core::visual::VisualModel);

private:
    typedef sofa::defaulttype::Vector4 Vector4;

    static const std::string DEPTH_OF_FIELD_VERTEX_SHADER;
    static const std::string DEPTH_OF_FIELD_FRAGMENT_SHADER;
    Data<double> zNear, zFar;
    Data<bool> useBBox;
    Data<bool> useRenderAR;
    Data<Vector4> BBox;
    bool postProcessEnabled;
    float *depths;
    cv::Mat depthmat,texturemat;
    Data<int> niterations;

public:
    ///Files where vertex shader is defined
    sofa::core::objectmodel::DataFileName vertFilename;
    ///Files where fragment shader is defined
    sofa::core::objectmodel::DataFileName fragFilename;

protected:
    RenderingManager();
    virtual ~RenderingManager();
public:
    void init() override;
    void reinit() override { };
    void initVisual() override;

    void preDrawScene(core::visual::VisualParams* vp) override;
    bool drawScene(core::visual::VisualParams* vp) override;
    void postDrawScene(core::visual::VisualParams* vp) override;

    void handleEvent(sofa::core::objectmodel::Event* event) override;
    void getDepths(cv::Mat &depths_){depths_ = depthmat;}
    void getTexture(cv::Mat &texture_){texture_ = texturemat;}
    double getZNear(){return zNear.getValue();}
    double getZFar(){return zFar.getValue();}
};


} //visualmodel

} //component

} //sofa

#endif
