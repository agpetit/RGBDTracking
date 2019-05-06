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

#define SOFA_RGBDTRACKING_CLOSESTPOINTFORCEFIELD_CPP

#include <sofa/core/ObjectFactory.h>
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

#include <SofaLoader/MeshObjLoader.h>
#include <limits>
#include <iterator>
#include <sofa/helper/gl/Color.h>

#ifdef Success
  #undef Success
#endif

#include "ClosestPointForceField.h"
#include "ImageConverter.h"


using std::cerr;
using std::endl;

namespace sofa
{

namespace component
{

namespace forcefield
{

    using namespace sofa::defaulttype;

      SOFA_DECL_CLASS(ClosestPointForceField)

      // Register in the Factory
      int ClosestPointForceFieldClass = core::RegisterObject("Compute forces based on closest points from/to a target surface/point set")
    #ifndef SOFA_FLOAT
        .add< ClosestPointForceField<Vec3dTypes> >()
    #endif
    #ifndef SOFA_DOUBLE
        .add< ClosestPointForceField<Vec3fTypes> >()
    #endif
    ;

    #ifndef SOFA_FLOAT
      template class SOFA_RGBDTRACKING_API ClosestPointForceField<Vec3dTypes>;
    #endif
    #ifndef SOFA_DOUBLE
      template class SOFA_RGBDTRACKING_API ClosestPointForceField<Vec3fTypes>;
    #endif

using namespace helper;


template <class DataTypes>
ClosestPointForceField<DataTypes>::ClosestPointForceField(core::behavior::MechanicalState<DataTypes> *mm )
    : Inherit(mm)
    , ks(initData(&ks,(Real)0.0,"stiffness","uniform stiffness for the all springs."))
    , kd(initData(&kd,(Real)0.0,"damping","uniform damping for the all springs."))
    , blendingFactor(initData(&blendingFactor,(Real)1,"blendingFactor","blending between projection (=0) and attraction (=1) forces."))
    , projectToPlane(initData(&projectToPlane,false,"projectToPlane","project closest points in the plane defined by the normal."))
    , springs(initData(&springs,"spring","index, stiffness, damping"))
    , cameraIntrinsicParameters(initData(&cameraIntrinsicParameters,Vector4(),"cameraIntrinsicParameters","camera parameters"))
    , sourceSurfacePositions(initData(&sourceSurfacePositions,"sourceSurface","Points of the surface of the source mesh."))
    , sourcePositions(initData(&sourcePositions,"sourcePositions","Points of the mesh."))
    , sourceContourPositions(initData(&sourceContourPositions,"sourceContourPositions","Contour points of the surface of the mesh."))
    , sourceContourNormals(initData(&sourceContourNormals,"sourceContourNormals","Normals to the contour points of the visible surface of the mesh."))
    , sourceVisible(initData(&sourceVisible,"sourceVisible","Visibility of the points of the surface of the mesh."))
    , indicesVisible(initData(&indicesVisible,"indicesVisible","Indices of the visible points of the mesh."))
    , sourceVisiblePositions(initData(&sourceVisiblePositions,"sourceVisiblePositions","Visible points of the surface of the mesh."))
    , sourceBorder(initData(&sourceBorder,"sourceBorder","Points of the border of the mesh."))
    , sourceWeights(initData(&sourceWeights,"sourceWeights","Weights for the nodes in the mesh."))
    , sourceTriangles(initData(&sourceTriangles,"sourceTriangles","Triangles of the source mesh."))
    , sourceNormals(initData(&sourceNormals,"sourceNormals","Normals of the source mesh."))
    , sourceSurfaceNormals(initData(&sourceSurfaceNormals,"sourceSurfaceNormals","Normals of the surface of the source mesh."))
    , targetPositions(initData(&targetPositions,"targetPositions","Points of the target point cloud."))
    , targetContourPositions(initData(&targetContourPositions,"targetContourPositions","Contour points of the target point cloud."))
    , targetBorder(initData(&targetBorder,"targetBorder","Contour flag of the target point cloud."))
    , targetWeights(initData(&targetWeights,"targetWeights","Weights for the points of the target point cloud."))
    , curvatures(initData(&curvatures,"curvatures","curvatures."))
    , outlierThreshold(initData(&outlierThreshold,(Real)7,"outlierThreshold","suppress outliers when distance > (meandistance + threshold*stddev)."))
    , rejectBorders(initData(&rejectBorders,false,"rejectBorders","ignore border vertices."))
    , useDistContourNormal(initData(&useDistContourNormal,false,"useVisible","Use the vertices of the visible surface of the source mesh"))
    , drawContour(initData(&drawContour,false,"drawContour"," "))
    , showArrowSize(initData(&showArrowSize,0.01f,"showArrowSize","size of the axis."))
    , drawMode(initData(&drawMode,0,"drawMode","The way springs will be drawn:\n- 0: Line\n- 1:Cylinder\n- 2: Arrow."))
    , drawColorMap(initData(&drawColorMap,true,"drawColorMap","Hue mapping of distances to closest point"))
    , theCloserTheStiffer(initData(&theCloserTheStiffer,false,"theCloserTheStiffer","Modify stiffness according to distance"))
    , useContour(initData(&useContour,false,"useContour","Emphasize forces close to the target contours"))
    , useContourWeight(initData(&useContourWeight,false,"useContourWeight","Emphasize forces close to the contours"))
    , useVisible(initData(&useVisible,true,"useVisible","Use the vertices of the viisible surface of the source mesh"))
    , useRealData(initData(&useRealData,true,"useRealData","Use real data"))
    , startimage(initData(&startimage,1,"startimage","Frame index to start registration"))
    , startimageklt(initData(&startimageklt,50,"startimageklt","Frame index to start registration with KLT features"))
    , niterations(initData(&niterations,3,"niterations","Number of iterations in the tracking process"))
    , dataPath(initData(&dataPath,"dataPath","Path for data writings",false))
    , windowKLT(initData(&windowKLT,5,"windowKLT","window for the KLT tracker"))
    , useKLTPoints(initData(&useKLTPoints, false,"useKLTPoints","Use KLT Points"))
{
    iter_im = 0;
}

template <class DataTypes>
ClosestPointForceField<DataTypes>::~ClosestPointForceField()
{
}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::reinit()
{

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue(); 	//RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
    this->clearSprings(x.size());

    for(unsigned int i=0;i<x.size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());

}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::init()
{

    this->Inherit::init();
    core::objectmodel::BaseContext* context = this->getContext();

        if(!(this->mstate)) this->mstate = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());

    // add a spring for every input point
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();//RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
    this->clearSprings(x.size());

    npoints = x.size();

        for(unsigned int i=0;i<x.size();i++)
            this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());

    glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);

    closestpoint = new sofa::core::objectmodel::ClosestPoint<DataTypes>;

    closestpoint->init();
    closestpoint->blendingFactor.setValue(blendingFactor.getValue());
    closestpoint->outlierThreshold.setValue(outlierThreshold.getValue());
    closestpoint->rejectBorders.setValue((rejectBorders.getValue()));
    closestpoint->useContour.setValue(useContour.getValue());
    closestpoint->useVisible.setValue(useVisible.getValue());
    closestpoint->useDistContourNormal.setValue(useDistContourNormal.getValue());
    Vector4 camParam = cameraIntrinsicParameters.getValue();

    rgbIntrinsicMatrix(0,0) = camParam[0];
    rgbIntrinsicMatrix(1,1) = camParam[1];
    rgbIntrinsicMatrix(0,2) = camParam[2];
    rgbIntrinsicMatrix(1,2) = camParam[3];
    closestpoint->rgbIntrinsicMatrix = rgbIntrinsicMatrix;

    // Tracker parameters
    tracker.setTrackerId(1);
    //tracker.setOnMeasureFeature(&modifyFeature);
    tracker.setMaxFeatures(1000);
    tracker.setWindowSize(10);
    tracker.setQuality(0.02);
    tracker.setMinDistance(10);
    tracker.setHarrisFreeParameter(0.04);
    tracker.setBlockSize(9);
    tracker.setUseHarris(1);
    tracker.setPyramidLevels(3);

    tracker1.setTrackerId(1);
    //tracker.setOnMeasureFeature(&modifyFeature);
    tracker1.setMaxFeatures(1000);
    tracker1.setWindowSize(10);
    tracker1.setQuality(0.02);
    tracker1.setMinDistance(10);
    tracker1.setHarrisFreeParameter(0.04);
    tracker1.setBlockSize(9);
    tracker1.setUseHarris(1);
    tracker1.setPyramidLevels(3);

}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::resetSprings()
{
    this->clearSprings(sourceVisiblePositions.getValue().size());
        for(unsigned int i=0;i<sourceVisiblePositions.getValue().size();i++)
            this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());
}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::mapKLTPointsTriangles ( sofa::helper::vector< tri > &triangles)
{
    int outside = 0;
        const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    //const sofa::core::topology::BaseMeshTopology::SeqTriangles& triangles = this->fromTopology->getTriangles();
    sofa::helper::vector<Mat3x3d> bases;
    sofa::helper::vector<Vector3> centers;

    helper::vector< bool > sourcevisible = sourceVisible.getValue();

    {
        {
            int c0 = triangles.size();
            bases.resize ( triangles.size());
            centers.resize ( triangles.size());
            for ( unsigned int t = 0; t < triangles.size(); t++ )
            {
                Mat3x3d m,mt;
                                double xt, yt;
                             Vector3 xim0,xim1,xim2;
                                 bool visible = true;
                        if (sourcevisible[triangles[t][2]] && sourcevisible[triangles[t][1]] && sourcevisible[triangles[t][0]])
                                 {
                        int x_u_2 = (int)(x[triangles[t][2]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[t][2]][2] + rgbIntrinsicMatrix(0,2));
                        int x_v_2 = (int)(x[triangles[t][2]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[t][2]][2] + rgbIntrinsicMatrix(1,2));
                        int x_u_1 = (int)(x[triangles[t][1]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[t][1]][2] + rgbIntrinsicMatrix(0,2));
                        int x_v_1 = (int)(x[triangles[t][1]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[t][1]][2] + rgbIntrinsicMatrix(1,2));

                        int x_u_0 = (int)(x[triangles[t][0]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[t][0]][2] + rgbIntrinsicMatrix(0,2));
                        int x_v_0 = (int)(x[triangles[t][0]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[t][0]][2] + rgbIntrinsicMatrix(1,2));

                        //std::cout << " x_u_2 " << x_u_2 << " " << x_u_1 << " " << x_u_0 << std::endl;

                        xim0[0] = x_u_0;
                        xim0[1] = x_v_0;
                        xim0[2] = 0;
                        xim1[0] = x_u_1;
                        xim1[1] = x_v_1;
                        xim1[2] = 0;
                        xim2[0] = x_u_2;
                        xim2[1] = x_v_2;
                        xim2[2] = 0;

                                m[0] = xim1-xim0;
                m[1] = xim2-xim0;
                m[2] = cross ( m[0],m[1] );
                mt.transpose ( m );
                bases[t].invert ( mt );
                centers[t] = ( xim0+xim1+xim2 ) /3;

                        int index = -1;
                        double distance = 1e10;

                        Real tt,uu,vv;
            /*for ( unsigned int t1 = 0; t1 < triangles.size(); t1++ )
            {
                Mat3x3d m,mt;
                                double xt, yt;
                             Vector3 xim0,xim1,xim2;

                                 bool intersect = true;

                                 {

    tt = 0; uu = 0; vv = 0;
    Vector3 edge1 = x[triangles[t1][1]] - x[triangles[t1][0]];
    Vector3 edge2 = x[triangles[t1][2]] - x[triangles[t1][0]];
    Vector3 tvec, pvec, qvec;
    Real det, inv_det;
                        //std::cout << " x_u_2 " << x_u_2 << " " << x_u_1 << " " << x_u_0 << std::endl;
    pvec = centers[t].cross(edge2);
    det = dot(edge1, pvec);
                                //std::cout << " dot " << det << std::endl;
    if(det<=1.0e-20 && det >=-1.0e-20)
    {
        intersect = false;
    }
    inv_det = 1.0 / det;
    tvec = - x[triangles[t1][0]];
    uu = dot(tvec, pvec) * inv_det;
    if (uu < -0.0000001 || uu > 1.0000001)
        intersect = false;
    qvec = tvec.cross(edge1);
    vv = dot(centers[t], qvec) * inv_det;
    if (vv < -0.0000001 || (uu + vv) > 1.0000001)
        intersect = false;
    tt = dot(edge2, qvec) * inv_det;

        //std::cout << " dot " << tt << std::endl;
    if (tt < 0.0000001 || tt!=tt || vv!=vv || uu!=uu)                   std::cout << " index 0 " <<  x[triangles[index][0]][0] << " x_v " << x[triangles[index][0]][1] << " x_v " << x[triangles[index][0]][2] << std::endl;

        intersect = false;
                        }

                        if (intersect)
                        {
                                if (centers[t][2] < x[triangles[t1][0]][2] && centers[t][2] < x[triangles[t1][1]][2] && centers[t][2] < x[triangles[t1][2]][2])
                                        visible = true;
                                        else visible = false;
                        }
                        if (visible){
                                sourceVisible[triangles[t][2]] = true;
                                sourceVisible[triangles[t][1]] = true;
                                sourceVisible[triangles[t][0]] = true;

                        }
                        else{
                                sourceVisible[triangles[t][2]] = false;
                                sourceVisible[triangles[t][1]] = false;
                                sourceVisible[triangles[t][0]] = false;
                        }

                }*/


                                 }
                        }


float xp, yp;
long id;

         for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
                Mat3x3d m,mt;
                                double xt, yt;
                             Vector3 xim0,xim1,xim2;
                                 const int kk = k;
                                 tracker.getFeature(kk, id, xp, yp);
                                int n0 = (int)yp;
                                int m0 = (int)xp;
                Vector3 pos;
                                pos[0] = xp;
                                pos[1] = yp;
                                pos[2] = 0;
                Vector3 coefs;

                int index = -1;
                double distance = 1e10;
                                //std::cout << "  xp yp " << xp << " " << yp << std::endl;
            for ( unsigned int t = 0; t < triangles.size(); t++ )
            {
                Mat3x3d m,mt;
                                double xt, yt;
                             Vector3 xim0,xim1,xim2;
                        if (sourcevisible[triangles[t][2]] && sourcevisible[triangles[t][1]] && sourcevisible[triangles[t][0]])
                                 {
                           // std::cout << " ok visible " << std::endl;
                        int x_u_2 = (int)(x[triangles[t][2]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[t][2]][2] + rgbIntrinsicMatrix(0,2));
                        int x_v_2 = (int)(x[triangles[t][2]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[t][2]][2] + rgbIntrinsicMatrix(1,2));
                        int x_u_1 = (int)(x[triangles[t][1]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[t][1]][2] + rgbIntrinsicMatrix(0,2));
                        int x_v_1 = (int)(x[triangles[t][1]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[t][1]][2] + rgbIntrinsicMatrix(1,2));

                        int x_u_0 = (int)(x[triangles[t][0]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[t][0]][2] + rgbIntrinsicMatrix(0,2));
                        int x_v_0 = (int)(x[triangles[t][0]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[t][0]][2] + rgbIntrinsicMatrix(1,2));

                        xim0[0] = x_u_0;
                        xim0[1] = x_v_0;
                        xim0[2] = 0;
                        xim1[0] = x_u_1;
                        xim1[1] = x_v_1;
                        xim1[2] = 0;
                        xim2[0] = x_u_2;
                        xim2[1] = x_v_2;
                        xim2[2] = 0;
                    Vec3d v = bases[t] * ( pos - xim0 );
                                        v[0] = ( pos - xim0 ).norm2()/(( pos - xim0 ).norm2() + ( pos - xim1 ).norm2() + ( pos - xim2 ).norm2());
                                        v[1] = ( pos - xim1 ).norm2()/(( pos - xim0 ).norm2() + ( pos - xim1 ).norm2() + ( pos - xim2 ).norm2());
                                        v[2] = ( pos - xim2 ).norm2()/(( pos - xim0 ).norm2() + ( pos - xim1 ).norm2() + ( pos - xim2 ).norm2());
                    double d = std::max ( std::max ( -v[0],-v[1] ),std::max ( ( v[2]<0?-v[2]:v[2] )-0.01,v[0]+v[1]-1 ) );
                    /*if ( d>0 )*/ d = ( pos-centers[t] ).norm2();
                    if ( d<distance ) { coefs = v; distance = d; index = t; }
                        }
                }
                if ( distance>0 )
                {
                    ++outside;
                }
                //if ( index < c0 )
                                {
                                        /*mapping map;
                                        map.coef = coefs;
                                        map.triangle = index;*/
                                        //std::cout << " mapping " << kk << " id " << id << " index " << index << " coefs " << coefs[0] << " " << coefs[1] << " " << coefs[2] << std::endl;
                                        mappingkltcoef[id] = coefs;
                                        mappingkltind[id] = index;
                        int x_u_0 = (int)(x[triangles[index][0]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[index][0]][2] + rgbIntrinsicMatrix(0,2));
                        int x_v_0 = (int)(x[triangles[index][0]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[index][0]][2] + rgbIntrinsicMatrix(1,2));
                        //std::cout << " x_u_0 " << x_u_0 << " x_v_0 " << x_v_0 <<  "  xp yp " << xp << " " << yp << " coefs " << coefs[0] << " " << coefs[1] << " " << coefs[2] << std::endl;

                                }

            }
        }
        }
}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::KLTPointsTo3D()
{

    int outside = 0;

        sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
        sofa::component::visualmodel::BaseCamera::SPtr currentCamera;// = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();
        root->get(currentCamera);

        float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
        float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy

        const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();

        VecCoord targetpos;
        targetpos.resize(tracker.getMaxFeatures());

        double znear = currentCamera->getZNear();
        double zfar = currentCamera->getZFar();

        Vector3 pos;
        Vector3 col;
        float xp, yp;
        long id;

         for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
                Mat3x3d m,mt;
                double xt, yt;
                Vector3 xim0,xim1,xim2;
                const int kk = k;
                tracker.getFeature(kk, id, xp, yp);
        //std::cout << " id " << id << " xp " << xp << " yp " << yp << std::endl;
                                int n0 = (int)yp;
                                 int m0 = (int)xp;
                                 float depthValue;
                                if (!useRealData.getValue())
                                depthValue = (float)depth.at<float>(2*yp,2*xp);
                                else depthValue = (float)depth.at<float>(yp,xp);
                        int avalue = (int)color.at<Vec4b>(yp,xp)[3];
                        if ( depthValue>0 && depthValue < 1)                // if depthValue is not NaN
                        {
                                double clip_z = (depthValue - 0.5) * 2.0;
                                if (!useRealData.getValue()) pos[2] = -2*znear*zfar/(clip_z*(zfar-znear)-(zfar+znear));
                                else pos[2] = depthValue;
                                pos[0] = (xp - rgbIntrinsicMatrix(0,2)) * pos[2] * rgbFocalInvertedX;
                                pos[1] = (yp - rgbIntrinsicMatrix(1,2)) * pos[2] * rgbFocalInvertedY;
                                targetpos[id]=pos;

                        }

            }
    const VecCoord&  p = targetpos;
        targetKLTPositions.setValue(p);

}




template <class DataTypes>
void ClosestPointForceField<DataTypes>::addForce(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{

    double timeaddforce = (double)getTickCount();
    addForceMesh(mparams, _f, _x, _v);
    std::cout << "TIME ADDFORCE " <<  (getTickCount() - timeaddforce)/getTickFrequency() << std::endl;

}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::addForceMesh(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
    int t = (int)this->getContext()->getTime();

    sofa::helper::vector< tri > triangles;
    triangles = sourceTriangles.getValue();

    bool reinitv = false;

    helper::vector< bool > sourcevisible = sourceVisible.getValue();
    helper::vector< int > indicesvisible = indicesVisible.getValue();
    helper::vector< bool > sourceborder = sourceBorder.getValue();

        if (t%niterations.getValue() == 0)
        {

            if (npoints != (this->mstate->read(core::ConstVecCoordId::position())->getValue()).size())
            {
                reinit();
                reinitv = true;
            }
            npoints = (this->mstate->read(core::ConstVecCoordId::position())->getValue()).size();

        }

    double time = (double)getTickCount();
    double timef0 = (double)getTickCount();

        if(ks.getValue()==0) return;

    VecDeriv& f = *_f.beginEdit();       //WDataRefVecDeriv f(_f);
    const VecCoord& x = _x.getValue();			//RDataRefVecCoord x(_x);
    const VecDeriv& v = _v.getValue();			//RDataRefVecDeriv v(_v);
    ReadAccessor< Data< VecCoord > > tn(targetNormals);
    ReadAccessor< Data< VecCoord > > tp(targetPositions);
    ReadAccessor< Data< VecCoord > > tcp(targetContourPositions);

        if (t%niterations.getValue() == 0)
        {
            f_.resize(f.size());
            x_.resize(x.size());
            v_.resize(v.size());
        }

    const vector<Spring>& s = this->springs.getValue();
    this->dfdx.resize(s.size());
    this->closestPos.resize(s.size());

    dfdx1.resize(s.size());
    m_potentialEnergy = 0;

    // get attraction/ projection factors
    Real attrF=(Real) blendingFactor.getValue();
    if(attrF<(Real)0.) attrF=(Real)0.;
    if(attrF>(Real)1.) attrF=(Real)1.;
    Real projF=((Real)1.-attrF);

    //closestpoint->updateClosestPointsGt();

    sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
    typename sofa::core::objectmodel::RGBDDataProcessing<DataTypes>::SPtr rgbddataprocessing;
    root->get(rgbddataprocessing);

    cv::Mat distimg = rgbddataprocessing->seg.distImage;
    cv::Mat foreg = rgbddataprocessing->foreground;


    if (useKLTPoints.getValue())
    {

        sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
        typename sofa::core::objectmodel::ImageConverter<DataTypes,DepthTypes>::SPtr imconv;
        root->get(imconv);

        typename sofa::core::objectmodel::DataIO<DataTypes>::SPtr dataio;
        root->get(dataio);

        bool okimages =false;

        //if (useRealData.getValue())
        {
        /*if (useSensor.getValue()){
                //color_1 = color.clone();
                //depth_1 = depth.clone();
            if(!((imconv->depth).empty()) && !((imconv->color).empty()))
            {
                if (scaleImages.getValue() > 1)
                {
                cv::resize(imconv->depth, depth, cv::Size(imconv->depth.cols/scaleImages.getValue(), imconv->depth.rows/scaleImages.getValue()), 0, 0);
                cv::resize(imconv->color, color, cv::Size(imconv->color.cols/scaleImages.getValue(), imconv->color.rows/scaleImages.getValue()), 0, 0);
                }
                else
                {
                color = imconv->color;
                depth = imconv->depth;
                }
                okimages = true;
            }
                //cv::imwrite("depth22.png", depth);
        }
        else*/
        {
                color = dataio->color;
                depth = dataio->depth;
                color_1 = dataio->color_1;
                okimages = true;
         }

        cvtColor(foreg,gray,CV_BGR2GRAY);
        cv::Mat gray1;
        cv::flip(gray,gray1,0);
        vpImageConvert::convert(gray,vpI);

        if (t == 0)
        {
            display.init(vpI, 100, 100,"Display...") ;
            // Display the image
            vpDisplay::display(vpI) ;
            vpDisplay::flush(vpI) ;

        // Point detection using Harris. In input we have an OpenCV image
        tracker.initTracking(gray);
        tracker1.initTracking(gray);
        tracker.display(vpI, vpColor::red);
        }
        else if (t >= 3){
            if (t == 3 || t%(windowKLT.getValue()*1) ==0 )
            {
    tracker.initTracking(gray);
    tracker1.initTracking(gray);
    mappingkltcoef.resize(tracker.getMaxFeatures());
    mappingkltind.resize(tracker.getMaxFeatures());
    mapKLTPointsTriangles(triangles);
            }

      cvtColor(foreg,gray,CV_BGR2GRAY);

    cv::flip(gray,gray1,0);
    vpImageConvert::convert(gray,vpI);
// Display the image
vpDisplay::display(vpI) ;
// Tracking of the detected points
tracker.track(gray);
tracker1.track(gray);

// Display the tracked points
tracker.display(vpI, vpColor::red);
vpDisplay::flush(vpI) ;
vpImage<vpRGBa> Icol;
cv::Mat icol;
display.getImage(Icol);
vpImageConvert::convert(Icol,icol);
/*imgklt = new cv::Mat;
//*imgl = color;
*imgklt = icol;

if (t%npasses == 0)
dataio->listimgklt.push_back(imgklt);*/

KLTPointsTo3D();

        }

    }
    }

    closestpoint->sourcePositions.setValue(this->mstate->read(core::ConstVecCoordId::position())->getValue());

    if (useVisible.getValue())
    closestpoint->sourceVisiblePositions.setValue(sourceVisiblePositions.getValue());

    closestpoint->targetPositions.setValue(targetPositions.getValue());
    closestpoint->targetBorder = targetBorder.getValue();
    closestpoint->sourceSurfacePositions.setValue(sourceSurfacePositions.getValue());
    closestpoint->sourceBorder = sourceBorder.getValue();
    helper::vector<double> sourceweights = sourceWeights.getValue();

    time = (double)getTickCount();

    double error = 0;
    int nerror = 0;

        if(tp.size()==0)
            for (unsigned int i=0; i<s.size(); i++) closestPos[i]=x[i];
        else
        {
            if (!useContour.getValue())
                closestpoint->updateClosestPoints();
            else
            {
                if ((targetContourPositions.getValue()).size() > 0 && (sourceContourPositions.getValue()).size()>0 )
                {
                closestpoint->targetContourPositions.setValue(targetContourPositions.getValue());
                closestpoint->sourceContourPositions.setValue(sourceContourPositions.getValue());
                closestpoint->sourceContourNormals.setValue(sourceContourNormals.getValue());
                closestpoint->updateClosestPointsContours();
                }
                else closestpoint->updateClosestPoints();
            }

            double timeClosestPoint = ((double)getTickCount() - timef0)/getTickFrequency();

            std::cout << "TIME CLOSESTPOINT " << timeClosestPoint << std::endl;
            indices = closestpoint->getIndices();

            // count number of attractors
            cnt.resize(s.size()); cnt.fill(0);
                if(attrF>0)
                    if (!useVisible.getValue())
                    {
                        for (unsigned int i=0; i<tp.size(); i++)
                            if(!closestpoint->targetIgnored[i])
                                cnt[closestpoint->closestTarget[i].begin()->second]++;

                    }
                    else
                    {
                            if ((sourceVisiblePositions.getValue()).size()>0)
                            {
                                for (unsigned int i=0; i<tp.size(); i++)
                                {
                                    if(!closestpoint->targetIgnored[i])
                                        cnt[indicesvisible[closestpoint->closestTarget[i].begin()->second]]++;
                                }
                            }
                            else for (unsigned int i=0; i<tp.size(); i++) cnt[closestpoint->closestTarget[i].begin()->second]++;

                    }

                    if(theCloserTheStiffer.getValue())
                    {
                        // find the min and the max distance value from source point to target point
                        min=0;
                        max=0;
                            for (unsigned int i=0; i<x.size(); i++)
                            {
                                if(min==0 || min>closestpoint->closestSource[i].begin()->first) min=closestpoint->closestSource[i].begin()->first;
                                if(max==0 || max<closestpoint->closestSource[i].begin()->first) max=closestpoint->closestSource[i].begin()->first;
                            }
                    }

            // compute targetpos = projF*closestto + attrF* sum closestfrom / count
            int ivis=0;
            int kk = 0;
            unsigned int id;
            {
            if(projF>0)
            {
                if (!useVisible.getValue())
                {
                    if (useContour.getValue())
                    {
                        if (targetContourPositions.getValue().size() > 0)
                        for (unsigned int i=0; i<s.size(); i++)
                        {
                            unsigned int id=closestpoint->closestSource[i].begin()->second;
                                if(!closestpoint->sourceIgnored[i])
                                {
                                    if(!sourceborder[i])
                                    {
                                        id=closestpoint->closestSource[i].begin()->second;
                                            if(projectToPlane.getValue() && tn.size()!=0)
                                                closestPos[i]=/*(1-(Real)sourceweights[i])**/(x[i]+tn[id]*dot(tp[id]-x[i],tn[id]))*projF;
                                            else
                                                closestPos[i]=/*(1-(Real)sourceweights[i])**/tp[id]*projF;

                                            if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                    }
                                    else
                                    {
                                        id=indices[kk];
                                        closestPos[i]=tcp[id]*projF;

                                            if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                        kk++;
                                    }

                                }
                                else
                                {
                                    closestPos[i]=x[i]*projF;
                                        if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                }
                            }
                        }
                        else
                        {
                            for (unsigned int i=0; i<s.size(); i++)
                            {
                                unsigned int id=closestpoint->closestSource[i].begin()->second;
                                    if(!closestpoint->sourceIgnored[i])
                                    {
                                        if(projectToPlane.getValue() && tn.size()!=0)	closestPos[i]=(x[i]+tn[id]*dot(tp[id]-x[i],tn[id]))*projF;
                                        else closestPos[i]=tp[id]*projF;
                                            /*if (sourceSurface[i])
                                            {
                                            closestPos[i]=(x[i]+ssn[i]*dot(tp[id]-x[i],ssn[i]))*projF;
                                            }
                                            else closestPos[i]=tp[id]*projF;*/
                                            if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                            //closestPos[i]+=x[i]*attrF;
                                    }
                                    else
                                    {
                                        closestPos[i]=x[i]*projF;
                                        if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                    }
                            }
                        }
                        }
                        else
                        {
                            if (useContour.getValue())
                            {
                                if (targetContourPositions.getValue().size() > 0)
                                for (unsigned int i=0; i<s.size(); i++)
                                {
                                    std::cout << "ii " << i<< " " << sourcevisible.size() << std::endl;


                                    //if(/*!closestpoint->sourceIgnored[i] &&*/ sourcevisible[i])
                                    {
                                        if (sourcevisible[i])
                                        {
                                            if(!sourceborder[i])
                                            {
                                            id=closestpoint->closestSource[ivis].begin()->second;

                                                if(projectToPlane.getValue() && tn.size()!=0)	closestPos[i]=/*(1-(Real)sourceweights[i])**/(x[i]+tn[id]*dot(tp[id]-x[i],tn[id]))*projF;
                                                else closestPos[i]=/*(1-(Real)sourceweights[i])**/tp[id]*projF;
                                                    if(!cnt[i]) closestPos[i]+=x[i]*attrF;

                                            }
                                            else
                                            {
                                                unsigned int id=indices[kk];
                                                closestPos[i]=tcp[id]*projF;
                                                    if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                                kk++;

                                            }
                                            ivis++;
                                        }
                                        else
                                        {
                                            closestPos[i]=x[i]*projF;
                                                if(!cnt[i]) {closestPos[i]+=x[i]*attrF;}
                                        }

                                    }
                                }
                            }
                            else
                            {
                                for (unsigned int i=0; i<s.size(); i++)
                                {
                                    //if(/*!closestpoint->sourceIgnored[i] &&*/ sourcevisible[i])
                                    {
                                        if (sourcevisible[i])
                                        {
                                            unsigned int id=closestpoint->closestSource[ivis].begin()->second;
                                                if(!closestpoint->sourceIgnored[ivis])
                                                {
                                                    closestPos[i]=tp[id]*projF;
                                                    error += this->computeError(x[i],tp[id]);
                                                    nerror++;
                                                }
                                                else closestPos[i]=x[i]*projF;
                                                ivis++;
                                        }
                                        else
                                        {
                                            closestPos[i]=x[i]*projF;
                                            //closestPos[i]=x[i];
                                            /*closestPos[i]=x[i]*projF;
                                            if(!cnt[i]) {closestPos[i]+=x[i]*attrF;}*/

                                        }
                                    if(!cnt[i]) {closestPos[i]+=x[i]*attrF;}

                                    }
                                }
                            }
                        }
                    }
                    else for (unsigned int i=0; i<s.size(); i++) { if(!cnt[i]) closestPos[i]=x[i]; else closestPos[i].fill(0); }

                // attraction

                if(attrF>0)
                        if(!useVisible.getValue())
                        {
                            for (unsigned int i=0; i<tp.size(); i++)
                            {
                                unsigned int id=closestpoint->closestTarget[i].begin()->second;
                                    if( !useContour.getValue() && !closestpoint->targetIgnored[i] && t > niterations.getValue()) //&& !targetBackground[i])
                                    {
                                        /*if (sourceSurface[id])
                                        {
                                                //std::cout << " ssn " << i << " " << ssn[id][1] << std::endl;

                                        closestPos[id]+=(tp[i]+ssn[id]*dot(x[id]-tp[i],ssn[id]))*attrF/(Real)cnt[id];
                                        }
                                        else*/
                                        {
                                                //closestPos[i]=(tp[i]+sn[id]*dot(x[id]-tp[i],sn[id]))*attrF/(Real)cnt[id];
                                                closestPos[id]+=tp[i]*attrF/(Real)cnt[id];
                                        }
                                    }
                                    else if (targetContourPositions.getValue().size() > 0)
                                    {
                                        if (useContour.getValue() && t > niterations.getValue() )//&& t%niterations.getValue() > 0)
                                            closestPos[id]+=tp[i]*attrF/(Real)cnt[id];
                                    }
                                }

                        }
                        else
                        {

                            {
                                if( !useContour.getValue())
                                {
                                    int kkt = 0;
                                    for (unsigned int i=0; i<tp.size(); i++)
                                    {
                                        unsigned int id=closestpoint->closestTarget[i].begin()->second;
                                        if(!closestpoint->targetIgnored[i]) //&& !targetBackground[i])
                                        {
                                            /*if (sourceSurface[id])
                                            {
                                                //std::cout << " ssn " << i << " " << ssn[id][1] << std::endl;
                                                closestPos[id]+=(tp[i]+ssn[id]*dot(x[id]-tp[i],ssn[id]))*attrF/(Real)cnt[id];
                                            }
                                            else*/
                                            {
                                                //closestPos[i]=(tp[i]+sn[id]*dot(x[id]-tp[i],sn[id]))*attrF/(Real)cnt[id];
                                                //closestPos[id]+=tp[i]*attrF/(Real)cnt[id];
                                                unsigned int id1 = indicesvisible[id];
                                                closestPos[id1]+=tp[i]*attrF/(Real)cnt[id1];
                                            }
                                        }
                                    }
                                }
                                else if (targetContourPositions.getValue().size() > 0)
                                {
                                    for (unsigned int i=0; i<tp.size(); i++)
                                    {
                                        unsigned int id=closestpoint->closestTarget[i].begin()->second;
                                        if(!closestpoint->targetIgnored[i])
                                        {
                                            unsigned int id1;
                                            //if (!rgbddataprocessing->targetBorder[i])
                                            id1 = indicesvisible[id];
                                            /*else
                                            {
                                                id1 = indicesTarget[kkt];
                                                kkt++;
                                            }*/
                                        //if (sourceborder[id1] && rgbddataprocessing->targetBorder[i])
                                        closestPos[id1]+=tp[i]*attrF/(Real)cnt[id1];
                                        }
                                        //if(rgbddataprocessing->targetBorder[i])
                                    }
                                }
                            }
                        }
                    }
                }


        ind = 0;
        int ivis =0;
        sourcew.resize(s.size());
        int npointspen = 0;



        /*VecCoord tpos = targetPositions.getValue();

        if (t%(niterations.getValue()) == 0){

            std::string pathvisible = "out/imagesPatient3Liver1_1/visiblemesh%06d.txt";
            std::string pathvisibleindices = "out/imagesPatient3Liver1_1/visiblemeshindices%06d.txt";
            std::string pathpcdmatch = "out/imagesPatient3Liver1_1/pcdmatch%06d.txt";
            std::string pathpcd = "out/imagesPatient3Liver1_1/pcd%06d.txt";


            char buf5[FILENAME_MAX];
            sprintf(buf5, pathvisible.c_str(),iter_im);
            std::string filename5(buf5);
            //std::cout << "Write: " << filename4 << std::endl;

            char buf6[FILENAME_MAX];
            sprintf(buf6, pathvisibleindices.c_str(), iter_im);
            std::string filename6(buf6);

            char buf7[FILENAME_MAX];
            sprintf(buf7, pathpcdmatch.c_str(), iter_im);
            std::string filename7(buf7);

            char buf8[FILENAME_MAX];
            sprintf(buf8, pathpcd.c_str(), iter_im);
            std::string filename8(buf8);

            iter_im++;


            filePCD.open(filename8.c_str(), std::ofstream::out);
            fileVisible.open(filename5.c_str(), std::ofstream::out);
            fileindicesVisible.open(filename6.c_str(), std::ofstream::out);
            filePCDMatch.open(filename7.c_str(), std::ofstream::out);
        for (unsigned int i=0; i<tpos.size(); i++)
        {
            filePCD << tpos[i][0];
            filePCD << " ";
            filePCD << tpos[i][1];
            filePCD << " ";
            filePCD << tpos[i][2];
            filePCD << "\n";

        }

        for (unsigned int i=0; i<s.size(); i++)
        {
            //serr<<"addForce() between "<<springs[i].m1<<" and "<<closestPos[springs[i].m1]<<sendl;
                if(sourcevisible[i])
                {
                fileVisible << x[s[i].m1][0];
                fileVisible << " ";
                fileVisible << x[s[i].m1][1];
                fileVisible << " ";
                fileVisible << x[s[i].m1][2];
                fileVisible << "\n";

                fileindicesVisible << s[i].m1;
                fileindicesVisible << "\n";

                filePCDMatch << this->closestPos[i][0];
                filePCDMatch << " ";
                filePCDMatch << this->closestPos[i][1];
                filePCDMatch << " ";
                filePCDMatch << this->closestPos[i][2];
                filePCDMatch << "\n";
                }

        }
        fileVisible.close();
        fileindicesVisible.close();
        filePCDMatch.close();
        filePCD.close();

        }*/


            for (unsigned int i=0; i<s.size(); i++)
            {
                //serr<<"addForce() between "<<springs[i].m1<<" and "<<closestPos[springs[i].m1]<<sendl;
                if (t%(niterations.getValue()) == 0 && t >= startimage.getValue() )
                    {
                    if (!useContourWeight.getValue())// || targetContourPositions.getValue().size()==0 )
                    this->addSpringForce(m_potentialEnergy,f,x,v, i, s[i]);
                    else this->addSpringForceWeight(m_potentialEnergy,f,x,v, i,ivis, s[i]);//this->addSpringForceWeight(m_potentialEnergy,f,x,v, i, s[i]);
                    if (sourcevisible[i]) ivis++;
                }
            }

            VecCoord  targetKLTPos;
            //if (useKLTPoints.getValue() && t%niterations.getValue() == 0 )
            if (useKLTPoints.getValue())
            targetKLTPos = targetKLTPositions.getValue();

            Vector3 coefs;
            int index;
            long id;
            float xp0, yp0;
            int x0, y0;

            /*sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
            typename sofa::core::objectmodel::RGBDDataProcessing<DataTypes>::SPtr rgbddataprocessing;
            root->get(rgbddataprocessing);

            cv::Mat distimg = rgbddataprocessing->seg.distImage;
            cv::Mat foreg = rgbddataprocessing->foreground;*/

            //cv::imwrite("foreg.png", distimg);

            //if (useKLTPoints.getValue() && t >= startimageklt.getValue() && t%(niterations.getValue()) == 0){
            if (useKLTPoints.getValue() && t >= startimageklt.getValue()){
            for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
                    //std::cout << " k " << k << std::endl;
                                          // std::cout << " index " << index << std::endl;
                                          // std::cout << " triangleindex " << triangles[index][0] << std::endl;
                                           //std::cout << " targetKLTPos " << triangles[index][0] << std::endl;

                           int x_u = (int)(x[triangles[index][0]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[index][0]][2] + rgbIntrinsicMatrix(0,2));
                           int x_v = (int)(x[triangles[index][0]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[index][0]][2] + rgbIntrinsicMatrix(1,2));
                       tracker.getFeature(k, id, xp0, yp0);
                       x0 = (int)xp0;
                       y0 = (int)yp0;
                           coefs = mappingkltcoef[id];
                           index = mappingkltind[id];
                           float depthValue;
                           int avalue;
                           if (!useRealData.getValue()) {depthValue = (float)depth.at<float>(2*yp0,2*xp0);
                                   avalue = (int)distimg.at<uchar>(yp0,xp0);
                           }
                           else {depthValue = (float)depth.at<float>(yp0,xp0);
                                   avalue = (int)distimg.at<uchar>(yp0,xp0);
                           }

               if( depthValue < 1 && depthValue > 0 && avalue > 10 && foreg.at<cv::Vec4b>(y0,x0)[3] > 0){
                           //std::cout << " x_u " <<targetKLTPos[id][2] << " kltfeat " << x[triangles[index][0]][2] << std::endl;
                   //std::cout << " x_u " << x_u << " x_v " << x_v << " x0 " << x0 << " y0 " << y0 << std::endl;
                   //std::cout << " x_u " << targetKLTPos[id][0] << " x_v " << targetKLTPos[id][1] << " x_v " << targetKLTPos[id][2] << std::endl;
                   //std::cout << " index 0 " <<  x[triangles[index][0]][0] << " x_v " << x[triangles[index][0]][1] << " x_v " << x[triangles[index][0]][2] << std::endl;
                   //std::cout << " index 1 " <<  x[triangles[index][1]][0] << " x_v " << x[triangles[index][1]][1] << " x_v " << x[triangles[index][1]][2] << std::endl;
                   //std::cout << " index 2 " <<  x[triangles[index][2]][0] << " x_v " << x[triangles[index][2]][1] << " x_v " << x[triangles[index][2]][2] << std::endl;

                           /*if (!meshprocessing->sourceBorder[triangles[index][0]])*/ this->addSpringForceKLTA(m_potentialEnergy,f,x,v, targetKLTPos[id], triangles[index][0], s[triangles[index][0]], 0.3*coefs[0]);
                           /*if (!meshprocessing->sourceBorder[triangles[index][1]])*/ this->addSpringForceKLTA(m_potentialEnergy,f,x,v, targetKLTPos[id], triangles[index][1], s[triangles[index][1]], 0.3*coefs[1]);
                           /*if (!meshprocessing->sourceBorder[triangles[index][2]])*/ this->addSpringForceKLTA(m_potentialEnergy,f,x,v, targetKLTPos[id], triangles[index][2], s[triangles[index][2]], 0.3*coefs[2]);
                           }

               }
                           //getchar();
                            }


    _f.endEdit();

}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
{
    int a = spring.m1;
    Coord u = this->closestPos[i]-p[a];
    Real d = u.norm();
    if( d>1.0e-4 )
    {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        potentialEnergy += elongation * elongation * spring.ks / 2;
        /*serr<<"addSpringForce, p = "<<p<<sendl;
        serr<<"addSpringForce, new potential energy = "<<potentialEnergy<<sendl;*/
        Deriv relativeVelocity = -v[a];
        Real elongationVelocity = dot(u,relativeVelocity);
        Real forceIntensity;
        if(theCloserTheStiffer.getValue())
        {
            Real ks_max=spring.ks;
            Real ks_min=spring.ks/10;
            Real ks_mod = ks_min*(max-elongation)/(max-min)+ks_max*(elongation-min)/(max-min);
            forceIntensity = (Real)(ks_mod*elongation+spring.kd*elongationVelocity);
        }
        else {
                if (elongation < 0.02)
                forceIntensity = (Real)(spring.ks*elongation+spring.kd*elongationVelocity);
                else forceIntensity = (Real)(spring.ks*elongation+spring.kd*elongationVelocity);
                }
        Deriv force = u*forceIntensity;
        f[a]+=force;
        Mat& m = this->dfdx[i];
        Real tgt = forceIntensity * inverseLength;
        for( int j=0; j<N; ++j )
        {
            // anisotropic
            //for( int k=0; k<N; ++k ) m[j][k] = tgt * u[j] * u[k];

            // isotropic
            for( int k=0; k<N; ++k ) m[j][k] = ((Real)spring.ks-tgt) * u[j] * u[k];
            m[j][j] += tgt;
        }
        //dfdx1[i] = m;

    }
    else // null length, no force and no stiffness
    {
        Mat& m = this->dfdx[i];
        for( int j=0; j<N; ++j )
        {
            for( int k=0; k<N; ++k )
            {
                m[j][k] = 0;
            }
        }
    }
}


template <class DataTypes>
void ClosestPointForceField<DataTypes>::addStoredSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
{
    int a = spring.m1;

        f[a]=f_[a];
        Mat& m = this->dfdx[i];
                m = dfdx1[i];
}

template <class DataTypes>
double ClosestPointForceField<DataTypes>::computeError(Vector3 sourcePoint, Vector3 targetPoint)
{
    //int a = spring.m1;
        Real elongation;
    Coord u = sourcePoint-targetPoint;
    Real d = u.norm();
    //if( d>1.0e-8 )
    {
        /*Real inverseLength = 1.0f/d;
        u *= inverseLength;*/
        elongation = (Real)d;
        }
        return elongation;
}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::addSpringForceWeight(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, int ivis, const Spring& spring)
{
    int a = spring.m1;
    Coord u = this->closestPos[i]-p[a];
    Real d = u.norm();
    if( d>1.0e-4 )
    {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        double stiffweight = 1;
        helper::vector< bool > sourceborder = sourceBorder.getValue();
        helper::vector< bool > sourcevisible = sourceVisible.getValue();
        helper::vector< int > indicesvisible = indicesVisible.getValue();
        helper::vector< double > sourceweights = sourceWeights.getValue();


        /*for (int k = 0; k < targetPositions.getValue().size(); k++){
                if(closestpoint->closestTarget[k].begin()->second == i || closestpoint->closestSource[i].begin()->second == k)
                {
                        if(sourceborder[i])
                stiffweight = (double)sourceWeights[i];
                    else stiffweight = (double)sourceWeights[i];
                //stiffweight = (double)targetWeights[k];
                }
                //else if (closestpoint->closestSource[i].begin()->second == k){
                //stiffweight = (double)combinedWeights[ind];
                //stiffweight = (double)sourceWeights[i]*targetWeights[k];
                //stiffweight = (double)targetWeights[k];
                //}
                ind++;
                }*/

        if(sourcevisible[i]){

        int k = (int)closestpoint->closestSource[ivis].begin()->second;

        if(!closestpoint->targetIgnored[k]) stiffweight = (double)targetWeights.getValue()[k];//*exp(-curvatures.getValue()[k]);
        else stiffweight = 1;

        }
        else stiffweight = 1;

        //std::cout << " ok ok " << stiffweight << std::endl;

        sourcew[i] = stiffweight;

                //if (sourceborder[i]) stiffweight*=1;
                //double stiffweight = (double)1/targetWeights[(int)closestpoint->closestSource[i].begin()->second];

        potentialEnergy += stiffweight*elongation * elongation * spring.ks / 2;
        /*serr<<"addSpringForce, p = "<<p<<sendl;
        serr<<"addSpringForce, new potential energy = "<<potentialEnergy<<sendl;*/
        Deriv relativeVelocity = -v[a];
        Real elongationVelocity = dot(u,relativeVelocity);
        Real forceIntensity;

        if(theCloserTheStiffer.getValue())
        {
            Real ks_max=stiffweight*spring.ks;
            Real ks_min=stiffweight*spring.ks/10;
            Real ks_mod = ks_min*(max-elongation)/(max-min)+ks_max*(elongation-min)/(max-min);
            forceIntensity = (Real)(ks_mod*elongation+spring.kd*elongationVelocity);
        }
        else {
                        if (elongation < 0.02)
                forceIntensity = (Real)(stiffweight*spring.ks*elongation+spring.kd*elongationVelocity);
                else forceIntensity = (Real)(stiffweight*spring.ks*elongation+spring.kd*elongationVelocity);
                }
        Deriv force = u*forceIntensity;
        f[a]+=force;
        Mat& m = this->dfdx[i];
        Real tgt = forceIntensity * inverseLength;
        for( int j=0; j<N; ++j )
        {
            // anisotropic
            //for( int k=0; k<N; ++k ) m[j][k] = tgt * u[j] * u[k];

            // isotropic
            for( int k=0; k<N; ++k ) m[j][k] = ((Real)stiffweight*spring.ks-tgt) * u[j] * u[k];
            m[j][j] += tgt;
        }
    }
    else // null length, no force and no stiffness
    {
        Mat& m = this->dfdx[i];
        for( int j=0; j<N; ++j )
        {
            for( int k=0; k<N; ++k )
            {
                m[j][k] = 0;
            }
        }
    }
}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::addSpringForceKLTA(double& potentialEnergy, VecDeriv& f, const  VecCoord& p,const VecDeriv& v, Coord& KLTtarget, int i, const Spring& spring, double coef)
{
    int a = spring.m1;
    Coord u = KLTtarget-p[a];
    Real d = u.norm();
    if( d>1.0e-4 )
    {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        potentialEnergy += coef*elongation * elongation *spring.ks / 2;
                //std::cout << " u " << potentialEnergy<< std::endl;
        /*          serr<<"addSpringForce, p = "<<p<<sendl;
        serr<<"addSpringForce, new potential energy = "<<potentialEnergy<<sendl;*/
        Deriv relativeVelocity = -v[a];
        Real elongationVelocity = dot(u,relativeVelocity);
        Real forceIntensity;
        if(theCloserTheStiffer.getValue())
        {
            Real ks_max=coef*spring.ks;
            Real ks_min=coef*spring.ks/10;
            Real ks_mod = ks_min*(max-elongation)/(max-min)+ks_max*(elongation-min)/(max-min);
            forceIntensity = (Real)(ks_mod*elongation+spring.kd*elongationVelocity);
        }
        else {
                        if (elongation < 0.02)
                forceIntensity = (Real)(coef*spring.ks*elongation+coef*spring.kd*elongationVelocity);
                else forceIntensity = (Real)(coef*spring.ks*elongation+coef*spring.kd*elongationVelocity);
                }
        Deriv force = u*forceIntensity;
        f[a]+=force;
        Mat& m = this->dfdx[i];
        Real tgt = forceIntensity * inverseLength;
                //std::cout << " u " << -v[a] << std::endl;
        for( int j=0; j<N; ++j )
        {
            // anisotropic
            //for( int k=0; k<N; ++k ) m[j][k] = tgt * u[j] * u[k];

            // isotropic
            for( int k=0; k<N; ++k ) m[j][k] = ((Real)coef*spring.ks-tgt) * u[j] * u[k];
            m[j][j] += tgt;
        }
                        //dfdx1[i] = m;

    }
    else // null length, no force and no stiffness
    {
        Mat& m = this->dfdx[i];
        for( int j=0; j<N; ++j )
        {
            /*for( int k=0; k<N; ++k )
            {
                m[j][k] = 0;
            }*/
        }
    }
}

template<class DataTypes>
void ClosestPointForceField<DataTypes>::addSpringDForce(VecDeriv& df,const  VecDeriv& dx, int i, const Spring& spring, double kFactor, double /*bFactor*/)
{
    const int a = spring.m1;
    const Coord d = -dx[a];
    Deriv dforce = this->dfdx[i]*d;
    dforce *= kFactor;
    df[a]+=dforce;
    //serr<<"addSpringDForce, a="<<a<<", b="<<b<<", dforce ="<<dforce<<sendl;
}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams,DataVecDeriv& _df , const DataVecDeriv&  _dx )
{

    VecDeriv& df = *_df.beginEdit();		//WDataRefVecDeriv df(_df);
    const VecDeriv&  dx = _dx.getValue();	// RDataRefVecDeriv dx(_dx);

    double kFactor 		 =  mparams->kFactor();
    double bFactor       =  mparams->bFactor();

    if(ks.getValue()==0) return;

    const vector<Spring>& s = this->springs.getValue();

    //serr<<"addDForce, dx = "<<dx<<sendl;
    //serr<<"addDForce, df before = "<<f<<sendl;
    for (unsigned int i=0; i<s.size(); i++)
    {
        this->addSpringDForce(df,dx, i, s[i], kFactor, bFactor);
    }
    //serr<<"addDForce, df = "<<f<<sendl;
    _df.endEdit();

}

template<class DataTypes>
void ClosestPointForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *m, SReal kFactor, unsigned int &offset)
{
            if(ks.getValue()==0) return;

    double kFact = kFactor;
    const vector<Spring >& ss = this->springs.getValue();
    const unsigned int n = ss.size() < this->dfdx.size() ? ss.size() : this->dfdx.size();
    for (unsigned int e=0; e<n; e++)
    {
        const Spring& s = ss[e];
        unsigned p1 = offset+Deriv::total_size*s.m1;
        const Mat& mt = this->dfdx[e];
        for(int i=0; i<N; i++)
            for (int j=0; j<N; j++)
            {
                Real k = (Real)(mt[i][j]*kFact);
                m->add(p1+i,p1+j, -k);
            }
    }
}

template<class DataTypes>
void ClosestPointForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    const Vec<4,float> c(1,0,0,1);
    std::vector< Vector3 > points;
    //const vector<Spring>& springs = this->springs.getValue();
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    points.resize(0);

    sofa::helper::vector< tri > triangles;
    triangles = sourceTriangles.getValue();

    std::cout << " XSIZE " << x.size() << std::endl;

    int t = (int)this->getContext()->getTime();

    if (targetPositions.getValue().size()>0 && sourceVisiblePositions.getValue().size()>0)
        for (unsigned int i=0; i<x.size(); i++)
        {
            //if(closestpoint->sourceIgnored[ivis] )
            {

                points.resize(0);
                Vector3 point = DataTypes::getCPos(x[i]);
                points.push_back(point);
               // std::cout << curvatures.getValue()[i] << std::endl;
                //if (targetWeights.getValue().size()>0) vparams->drawTool()->drawPoints(points, 10, sofa::defaulttype::Vec<4,float>(0.5*sourcew[i],0,0,1));
            }
        }

    points.resize(0);

    VecCoord  targetKLTPos;
    if (useKLTPoints.getValue() && t%niterations.getValue() == 0 )
    targetKLTPos = targetKLTPositions.getValue();

    Vector3 coefs;
    int index;
    long id;
    float xp0, yp0;
    int x0, y0;

    sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
    typename sofa::core::objectmodel::RGBDDataProcessing<DataTypes>::SPtr rgbddataprocessing;
    root->get(rgbddataprocessing);

    cv::Mat distimg = rgbddataprocessing->seg.distImage;
    cv::Mat foreg = rgbddataprocessing->foreground;

    if (useKLTPoints.getValue() && t >= startimageklt.getValue() && t%(niterations.getValue()) == 0){
    for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
            //std::cout << " k " << k << std::endl;

               tracker.getFeature(k, id, xp0, yp0);
               x0 = (int)xp0;
               y0 = (int)yp0;
                   coefs = mappingkltcoef[id];
                   index = mappingkltind[id];
                   int x_u = (int)(x[triangles[index][0]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[index][0]][2] + rgbIntrinsicMatrix(0,2));
                   int x_v = (int)(x[triangles[index][0]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[index][0]][2] + rgbIntrinsicMatrix(1,2));
                   float depthValue;
                   int avalue;
                   if (!useRealData.getValue()) {depthValue = (float)depth.at<float>(2*yp0,2*xp0);
                           avalue = (int)distimg.at<uchar>(yp0,xp0);
                   }
                   else {depthValue = (float)depth.at<float>(yp0,xp0);
                           avalue = (int)distimg.at<uchar>(yp0,xp0);
                   }


       if( depthValue < 1 && depthValue > 0 && avalue > 2 && foreg.at<cv::Vec4b>(y0,x0)[3] > 0){

           //if (useContour.getValue() && drawContour.getValue())
           Vector3 point1 = DataTypes::getCPos(targetKLTPos[id]);
           Vector3 point2 = (coefs[0]*DataTypes::getCPos(x[triangles[index][0]]) + coefs[1]*DataTypes::getCPos(x[triangles[index][1]]) + coefs[2]*DataTypes::getCPos(x[triangles[index][2]]));
           //std::cout << " pt " << point2[0] << " " << point2[1] << " " << point2[2] << std::endl;
           {
               points.push_back(point1);
               points.push_back(point2);
           }

           }

       }
    //for (unsigned int i=0;i<points.size()/2;++i) vparams->drawTool()->drawArrow(points[2*i+1], points[2*i], 0.0005, c);

                   //getchar();
                    }

   /* if (targetPositions.getValue().size()>0 && sourceVisiblePositions.getValue().size()>0)
    for (unsigned int i=0; i<x.size(); i++)
    {
        //if(closestpoint->sourceIgnored[ivis] )
        {

            bool dispP = true;
            //if (useContour.getValue() && drawContour.getValue())
            Vector3 point1 = DataTypes::getCPos(x[i]);
            Vector3 point2 = DataTypes::getCPos(this->closestPos[i]);
            //std::cout << " pt " << point2[0] << " " << point2[1] << " " << point2[2] << std::endl;

            if (dispP)
            {
                points.push_back(point1);
                points.push_back(point2);
            }
        }
    }
    if (drawMode.getValue() == 1)
                    for (unsigned int i=0;i<points.size()/2;++i) vparams->drawTool()->drawCylinder(points[2*i+1], points[2*i], showArrowSize.getValue(), c);
    else if (drawMode.getValue() == 2)
                    for (unsigned int i=0;i<points.size()/2;++i) vparams->drawTool()->drawArrow(points[2*i+1], points[2*i], 0.0005, c);

    std::cout << " XSIZE " << x.size() << std::endl;*/

}

}
}
} // namespace sofa

//#endif  /* SOFA_COMPONENT_CLOSESTPOINTFORCEFIELD_INL */


