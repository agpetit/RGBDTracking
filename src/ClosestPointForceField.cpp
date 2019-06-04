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

#include "ClosestPointForceField.h"
#include "ImageConverter.h"


using std::cerr;
using std::endl;


namespace sofa {

namespace rgbdtracking {

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(ClosestPointForceField)

// Register in the Factory
int ClosestPointForceFieldClass = core::RegisterObject("Compute forces based on closest points from/to a target surface/point set")
//#ifndef SOFA_FLOAT
    .add< ClosestPointForceField<sofa::defaulttype::Vec3dTypes> >()
//#endif
//#ifndef SOFA_DOUBLE
//    .add< ClosestPointForceField<Vec3fTypes> >()
//#endif
;

//#ifndef SOFA_FLOAT
    template class SOFA_RGBDTRACKING_API ClosestPointForceField<Vec3dTypes>;
//#endif
//#ifndef SOFA_DOUBLE
//    template class SOFA_RGBDTRACKING_API ClosestPointForceField<Vec3fTypes>;
//#endif

using namespace helper;


template <class DataTypes>
ClosestPointForceField<DataTypes>::ClosestPointForceField(core::behavior::MechanicalState<DataTypes> *mm )
    : Inherit(mm)
    , l_dataio(initLink("dataio", "Link to dataio component"))
    , l_rgbddataprocess(initLink("target", "Link to RGBDDataProcessing component"))
    , l_meshprocessing(initLink("source", "Link to MeshProcessing component"))

    , ks(initData(&ks,(Real)0.0,"stiffness","uniform stiffness for the all springs."))
    , kd(initData(&kd,(Real)0.0,"damping","uniform damping for the all springs."))
    , projectToPlane(initData(&projectToPlane,false,"projectToPlane","project closest points in the plane defined by the normal."))
    , springs(initData(&springs,"spring","index, stiffness, damping"))

    //closest point parameters
    , blendingFactor(initData(&blendingFactor,(Real)1,"blendingFactor","blending between projection (=0) and attraction (=1) forces."))
    , outlierThreshold(initData(&outlierThreshold,(Real)7,"outlierThreshold","suppress outliers when distance > (meandistance + threshold*stddev)."))
    , rejectBorders(initData(&rejectBorders,false,"rejectBorders","ignore border vertices."))
    , useContour(initData(&useContour,false,"useContour","Emphasize forces close to the target contours"))
    , useVisible(initData(&useVisible,true,"useVisible","Use the vertices of the viisible surface of the source mesh"))
    , useDistContourNormal(initData(&useDistContourNormal,false,"useVisible","Use the vertices of the visible surface of the source mesh"))
    , cameraIntrinsicParameters(initData(&cameraIntrinsicParameters,Vector4(),"cameraIntrinsicParameters","camera parameters"))


    , theCloserTheStiffer(initData(&theCloserTheStiffer,false,"theCloserTheStiffer","Modify stiffness according to distance"))
    , useContourWeight(initData(&useContourWeight,false,"useContourWeight","Emphasize forces close to the contours"))
    , startimage(initData(&startimage,1,"startimage","Frame index to start registration"))
    , niterations(initData(&niterations,3,"niterations","Number of iterations in the tracking process"))

    // KLT
    , startimageklt(initData(&startimageklt,50,"startimageklt","Frame index to start registration with KLT features"))
    , windowKLT(initData(&windowKLT,5,"windowKLT","window for the KLT tracker"))
    , useKLTPoints(initData(&useKLTPoints, false,"useKLTPoints","Use KLT Points"))
    , targetKLTPositions(initData(&targetKLTPositions, "targetkltpositions", "target KLT positions"))
{}

template <class DataTypes>
ClosestPointForceField<DataTypes>::~ClosestPointForceField()
{
}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::reinit()
{
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    //RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));

    this->clearSprings(x.size());

    for(unsigned int i=0;i<x.size();i++) {
        this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());
    }

}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::init()
{

    this->Inherit::init();
    core::objectmodel::BaseContext* context = this->getContext();

    if(!(this->mstate)) {
        this->mstate = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());
    }
    // add a spring for every input point
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    //RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
    this->clearSprings(x.size());

    npoints = x.size();

    for(unsigned int i=0;i<x.size();i++) {
        this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());
    }

    glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);

    closestpoint = new rgbdtracking::ClosestPoint<DataTypes>;

    closestpoint->init();
    closestpoint->blendingFactor.setValue(blendingFactor.getValue());
    closestpoint->outlierThreshold.setValue(outlierThreshold.getValue());
    closestpoint->rejectBorders.setValue((rejectBorders.getValue()));
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
void ClosestPointForceField<DataTypes>::mapKLTPointsTriangles ( sofa::helper::vector< tri > &triangles)
{
    int outside = 0;
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    //const sofa::core::topology::BaseMeshTopology::SeqTriangles& triangles = this->fromTopology->getTriangles();
    sofa::helper::vector<Mat3x3d> bases;
    sofa::helper::vector<Vector3> centers;

    helper::vector< bool > sourcevisible = l_meshprocessing->sourceVisible.getValue();



    int c0 = triangles.size();
    bases.resize ( triangles.size());
    centers.resize ( triangles.size());
    for ( unsigned int t = 0; t < triangles.size(); t++ ) {
        Mat3x3d m,mt;
        Vector3 xim0,xim1,xim2;
        if (sourcevisible[triangles[t][2]] &&
            sourcevisible[triangles[t][1]] &&
            sourcevisible[triangles[t][0]]
        ) {
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
        }
    }


    float xp, yp;
    long id;

    for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
        const int kk = k;
        tracker.getFeature(kk, id, xp, yp);
        Vector3 pos = Vector3 (xp, yp, 0);

        Vector3 coefs;

        int index = -1;
        double distance = 1e10;
        //std::cout << "  xp yp " << xp << " " << yp << std::endl;
        for ( unsigned int t = 0; t < triangles.size(); t++ ) {
            if (sourcevisible[triangles[t][2]] &&
                sourcevisible[triangles[t][1]] &&
                sourcevisible[triangles[t][0]]
            ) {
                Vector3 xim[3];
               // std::cout << " ok visible " << std::endl;
                for (size_t i = 0 ; i < 3 ; i++) {
                    int x_u = (int)(x[triangles[t][i]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[t][i]][2] + rgbIntrinsicMatrix(0,2));
                    int x_v = (int)(x[triangles[t][i]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[t][i]][2] + rgbIntrinsicMatrix(1,2));
                    xim[i] = Vector3(x_u, x_v, 0) ;
                }

                Vec3d v = bases[t] * ( pos - xim[0] );
                for (size_t i = 0 ; i < 3 ; i++) {
                    v[i] = ( pos - xim[i] ).norm2()/(( pos - xim[0] ).norm2() + ( pos - xim[1] ).norm2() + ( pos - xim[2] ).norm2());
                }

                double d = std::max ( std::max ( -v[0],-v[1] ),std::max ( ( v[2]<0?-v[2]:v[2] )-0.01,v[0]+v[1]-1 ) );
                /*if ( d>0 )*/
                d = ( pos-centers[t] ).norm2();
                if ( d<distance ) {
                    coefs = v;
                    distance = d;
                    index = t;
                }
            }
        }
        if ( distance>0 ) {
            ++outside;
        }

        mappingkltcoef[id] = coefs;
        mappingkltind[id] = index;
    }
}


template <class DataTypes>
void ClosestPointForceField<DataTypes>::KLTPointsTo3D()
{

    float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
    float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy

    VecCoord targetpos;
    targetpos.resize(tracker.getMaxFeatures());


    for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
        long id;
        float xp, yp;
        Vector3 pos;
        const int kk = k;
        tracker.getFeature(kk, id, xp, yp);

        cv::Mat depth = l_dataio->depth ;
        float depthValue = (float)depth.at<float>(yp,xp);
        if ( depthValue>0 && depthValue < 1) {
            // if depthValue is not NaN
            pos[2] = depthValue;
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

    helper::AdvancedTimer::stepBegin("AddForceMesh") ;
    addForceMesh(mparams, _f, _x, _v);
    helper::AdvancedTimer::stepEnd("AddForceMesh") ;

}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::addForceMesh(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{

    std::cout << "IN" << std::endl ;
    int t = (int)this->getContext()->getTime();

    sofa::helper::vector< tri > triangles;
    triangles = l_meshprocessing->sourceTriangles.getValue();

    helper::vector< bool > sourcevisible = l_meshprocessing->sourceVisible.getValue();
    helper::vector< int > indicesvisible = l_meshprocessing->indicesVisible.getValue();
    helper::vector< bool > sourceborder = l_meshprocessing->sourceBorder.getValue();

    if (t%niterations.getValue() == 0) {
        if (npoints != (this->mstate->read(core::ConstVecCoordId::position())->getValue()).size()) {
            reinit();
        }
        npoints = (this->mstate->read(core::ConstVecCoordId::position())->getValue()).size();
    }

    VecDeriv& f = *_f.beginEdit();       //WDataRefVecDeriv f(_f);
    const VecCoord& x = _x.getValue();			//RDataRefVecCoord x(_x);
    const VecDeriv& v = _v.getValue();			//RDataRefVecDeriv v(_v);
    ReadAccessor< Data< VecCoord > > tn(l_rgbddataprocess->targetNormals);
    ReadAccessor< Data< VecCoord > > tp(l_rgbddataprocess->targetPositions);
    ReadAccessor< Data< VecCoord > > tcp(l_rgbddataprocess->targetContourPositions);

    if (t%niterations.getValue() == 0) {
        f_.resize(f.size());
        x_.resize(x.size());
        v_.resize(v.size());
    }
    std::cout << "IN" << std::endl ;

    const vector<Spring>& s = this->springs.getValue();
    this->dfdx.resize(s.size());
    this->closestPos.resize(s.size());

    dfdx1.resize(s.size());
    double potentialEnergy = 0;

    // get attraction/ projection factors
    Real attrF=(Real) blendingFactor.getValue();
    if(attrF<(Real)0.) {
        attrF=(Real)0.;
    } else if(attrF>(Real)1.) {
        attrF=(Real)1.;
    }
    Real projF=((Real)1.-attrF);

    //closestpoint->updateClosestPointsGt();

    cv::Mat distimg = l_rgbddataprocess->seg.distImage;
    cv::Mat foreg = l_rgbddataprocess->foreground;

    cv::Mat depth = l_dataio->depth;
    if (useKLTPoints.getValue()) {
        cv::Mat gray;

        cvtColor(foreg,gray,CV_BGR2GRAY);
        cv::Mat gray1;
        cv::flip(gray,gray1,0);
        vpImage<unsigned char> vpI ;
        vpImageConvert::convert(gray,vpI);

        vpDisplayX display;
        if (t == 0) {
            display.init(vpI, 100, 100,"Display...") ;
            // Display the image
            vpDisplay::display(vpI) ;
            vpDisplay::flush(vpI) ;

            // Point detection using Harris. In input we have an OpenCV image
            tracker.initTracking(gray);
            tracker1.initTracking(gray);
            tracker.display(vpI, vpColor::red);
        } else if (t >= 3){
            if (t == 3 || t%(windowKLT.getValue()*1) ==0 ) {
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

            KLTPointsTo3D();

        }
    }

    closestpoint->sourcePositions.setValue(this->mstate->read(core::ConstVecCoordId::position())->getValue());

    if (useVisible.getValue()) {
        closestpoint->sourceVisiblePositions.setValue(l_meshprocessing->sourceVisiblePositions.getValue());
    }

    closestpoint->targetPositions.setValue(l_rgbddataprocess->targetPositions.getValue());
    closestpoint->targetBorder = l_rgbddataprocess->targetBorder.getValue();
    closestpoint->sourceSurfacePositions.setValue(l_meshprocessing->sourceSurfacePositions.getValue());
    closestpoint->sourceBorder = l_meshprocessing->sourceBorder.getValue();

//    time = (double)getTickCount();

    double error = 0;
    int nerror = 0;

    if(tp.size()==0) {
        for (unsigned int i=0; i<s.size(); i++) {
            closestPos[i]=x[i];
        }
    } else {
        //closestpoint process
        helper::AdvancedTimer::stepBegin("ClosestPoint") ;
        if (!useContour.getValue()) {
            closestpoint->updateClosestPoints();
        } else if ((l_rgbddataprocess->targetContourPositions.getValue()).size() > 0 &&
                   (l_meshprocessing->sourceContourPositions.getValue()).size()>0 ) {
            closestpoint->targetContourPositions.setValue(l_rgbddataprocess->targetContourPositions.getValue());
            closestpoint->sourceContourPositions.setValue(l_meshprocessing->sourceContourPositions.getValue());
            closestpoint->sourceContourNormals.setValue(l_meshprocessing->sourceContourNormals.getValue());
            closestpoint->updateClosestPointsContours();
        } else {
            closestpoint->updateClosestPoints();
        }
        helper::AdvancedTimer::stepEnd("ClosestPoint") ;

        std::vector<int> indices = closestpoint->getIndices();

        // count number of attractors
        vector<unsigned int>  cnt;
        cnt.resize(s.size()); cnt.fill(0);
        if(attrF>0) {
            for (unsigned int i=0; i<tp.size(); i++) {
                if (useVisible.getValue() &&
                    l_meshprocessing->sourceVisiblePositions.getValue().size()>0 &&
                    !closestpoint->targetIgnored[i]
                ) {
                    cnt[indicesvisible[closestpoint->closestTarget[i].begin()->second]]++;
                } else {
                    cnt[closestpoint->closestTarget[i].begin()->second]++;
                }
            }
        }

        if(theCloserTheStiffer.getValue()) {
            // find the min and the max distance value from source point to target point
            min=0;
            max=0;
            for (unsigned int i=0; i<x.size(); i++) {
                if(min==0 || min>closestpoint->closestSource[i].begin()->first) {
                    min=closestpoint->closestSource[i].begin()->first;
                }
                if(max==0 || max<closestpoint->closestSource[i].begin()->first) {
                    max=closestpoint->closestSource[i].begin()->first;
                }
            }
        }

        // compute targetpos = projF*closestto + attrF* sum closestfrom / count
        int ivis=0;

        if(projF>0) {
            int kk = 0;
            unsigned int id;
            if (!useVisible.getValue() &&
                useContour.getValue() &&
                l_rgbddataprocess->targetContourPositions.getValue().size() > 0
            ) {
                for (unsigned int i=0; i<s.size(); i++) {
                    unsigned int id=closestpoint->closestSource[i].begin()->second;
                    if(!closestpoint->sourceIgnored[i]) {
                        if(!sourceborder[i]) {
                            id=closestpoint->closestSource[i].begin()->second;
                            if(projectToPlane.getValue() && tn.size()!=0) {
                                closestPos[i]=(x[i]+tn[id]*dot(tp[id]-x[i],tn[id]));
                            } else {
                                closestPos[i]=tp[id];
                            }
                        } else {
                            id=indices[kk++];
                            closestPos[i]=tcp[id];
                        }
                    } else {
                        closestPos[i]=x[i] ;
                    }
                    closestPos[i] *= projF ;
                    if(!cnt[i]) {
                        closestPos[i]+=x[i]*attrF;
                    }
                }
            } else if (!useVisible.getValue()) {
                for (unsigned int i=0; i<s.size(); i++) {
                    unsigned int id=closestpoint->closestSource[i].begin()->second;
                    if(!closestpoint->sourceIgnored[i]) {
                        if(projectToPlane.getValue() && tn.size()!=0) {
                            closestPos[i]=(x[i]+tn[id]*dot(tp[id]-x[i],tn[id]))*projF;
                        } else {
                            closestPos[i]=tp[id]*projF;
                        }
                    } else {
                        closestPos[i]=x[i]*projF;
                    }
                    if(!cnt[i]) {
                        closestPos[i]+=x[i]*attrF;
                    }
                }
            } else if (useContour.getValue() && l_rgbddataprocess->targetContourPositions.getValue().size() > 0) {
                for (unsigned int i=0; i<s.size(); i++) {
                    std::cout << "ii " << i<< " " << sourcevisible.size() << std::endl;

                    if (sourcevisible[i] && !sourceborder[i]) {
                        id=closestpoint->closestSource[ivis++].begin()->second;

                        if(projectToPlane.getValue() && tn.size()!=0) {
                            closestPos[i]=(x[i]+tn[id]*dot(tp[id]-x[i],tn[id]))*projF;
                        }else {
                            closestPos[i]=tp[id]*projF;
                        }
                    } else if (sourcevisible[i]) {
                        id=indices[kk++];
                        closestPos[i]=tcp[id]*projF;
                    } else {
                        closestPos[i]=x[i]*projF;
                    }
                    if(!cnt[i]) {
                         closestPos[i]+=x[i]*attrF;
                     }
                }
            } else {
                for (unsigned int i=0; i<s.size(); i++) {
                    if (sourcevisible[i] && !closestpoint->sourceIgnored[ivis]) {
                        unsigned int id=closestpoint->closestSource[ivis++].begin()->second;
                        closestPos[i]=tp[id]*projF;
                        error += this->computeError(x[i],tp[id]);
                        nerror++;
                    } else {
                        closestPos[i]=x[i]*projF;
                    }
                    if(!cnt[i]) {
                    // is this necessary ?
                        closestPos[i]+=x[i]*attrF;
                    }
                }
            }
        } else {
            for (unsigned int i=0; i<s.size(); i++) {
                if(!cnt[i]) {
                    closestPos[i]=x[i];
                }else {
                    closestPos[i].fill(0);
                }
            }
        }

        // attraction
        if(attrF>0 && t > niterations.getValue() &&
          (!useVisible.getValue() || l_rgbddataprocess->targetContourPositions.getValue().size() > 0)
        ) {
            for (unsigned int i=0; i<tp.size(); i++) {
                unsigned int id=closestpoint->closestTarget[i].begin()->second;
                if( !closestpoint->targetIgnored[i]) {
                    id = (!useContour.getValue()) ? indicesvisible[id] : id ;
                    closestPos[id]+=tp[i]*attrF/(Real)cnt[id];
                }
            }
        }

    }

    int ivis =0;

    for (unsigned int i=0; i<s.size(); i++) {
        //serr<<"addForce() between "<<springs[i].m1<<" and "<<closestPos[springs[i].m1]<<sendl;
        if (t%(niterations.getValue()) == 0 &&
            t >= startimage.getValue()
        ) {
            if (!useContourWeight.getValue()) {
            // || targetContourPositions.getValue().size()==0 )
                this->addSpringForce(potentialEnergy,f,x,v, i, s[i]);
            } else {
                this->addSpringForceWeight(potentialEnergy,f,x,v, i,ivis, s[i]);
            }
            if (sourcevisible[i]) {
                ivis++;
            }
        }
    }


    if (useKLTPoints.getValue() &&
        t >= startimageklt.getValue()
    ){
        VecCoord  targetKLTPos = targetKLTPositions.getValue();
        for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
            long id;
            float xp0, yp0;
            tracker.getFeature(k, id, xp0, yp0);
            int x0 = (int)xp0,
                y0 = (int)yp0;
            Vector3 coefs = mappingkltcoef[id];
            int index = mappingkltind[id];

            float depthValue = (float)depth.at<float>(yp0,xp0);
            int avalue = (int)distimg.at<uchar>(yp0,xp0);

           if( depthValue < 1 &&
               depthValue > 0 &&
               avalue > 10 &&
               foreg.at<cv::Vec4b>(y0,x0)[3] > 0
            ){
               this->addSpringForceKLTA(potentialEnergy,f,x,v, targetKLTPos[id], triangles[index][0], s[triangles[index][0]], 0.3*coefs[0]);
               this->addSpringForceKLTA(potentialEnergy,f,x,v, targetKLTPos[id], triangles[index][1], s[triangles[index][1]], 0.3*coefs[1]);
               this->addSpringForceKLTA(potentialEnergy,f,x,v, targetKLTPos[id], triangles[index][2], s[triangles[index][2]], 0.3*coefs[2]);
            }

        }
       //getchar();
    }
    m_potentialEnergy = potentialEnergy ;
    _f.endEdit();
}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
{
    int a = spring.m1;
    Coord u = this->closestPos[i]-p[a];
    Real d = u.norm();
    if( d>1.0e-4 ) {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        potentialEnergy += elongation * elongation * spring.ks / 2;
        /*serr<<"addSpringForce, p = "<<p<<sendl;
        serr<<"addSpringForce, new potential energy = "<<potentialEnergy<<sendl;*/
        Deriv relativeVelocity = -v[a];
        Real elongationVelocity = dot(u,relativeVelocity);
        Real forceIntensity;
        if(theCloserTheStiffer.getValue()){
            Real ks_max=spring.ks;
            Real ks_min=spring.ks/10;
            Real ks_mod = ks_min*(max-elongation)/(max-min)+ks_max*(elongation-min)/(max-min);
            forceIntensity = (Real)(ks_mod*elongation+spring.kd*elongationVelocity);
        } else if (elongation < 0.02) {
            forceIntensity = (Real)(spring.ks*elongation+spring.kd*elongationVelocity);
        } else {
            forceIntensity = (Real)(spring.ks*elongation+spring.kd*elongationVelocity);
        }

        Deriv force = u*forceIntensity;
        f[a]+=force;
        Mat& m = this->dfdx[i];
        Real tgt = forceIntensity * inverseLength;
        for( int j=0; j<N; ++j ) {
            // anisotropic
            //for( int k=0; k<N; ++k ) m[j][k] = tgt * u[j] * u[k];

            // isotropic
            for( int k=0; k<N; ++k ) {
                m[j][k] = ((Real)spring.ks-tgt) * u[j] * u[k];
            }
            m[j][j] += tgt;
        }
        //dfdx1[i] = m;

    } else {
    // null length, no force and no stiffness
        Mat& m = this->dfdx[i];
        for( int j=0; j<N; ++j ) {
            for( int k=0; k<N; ++k ) {
                m[j][k] = 0;
            }
        }
    }
}

template <class DataTypes>
double ClosestPointForceField<DataTypes>::computeError(Vector3 sourcePoint, Vector3 targetPoint)
{
    //int a = spring.m1;
    Real elongation;
    Coord u = sourcePoint-targetPoint;
    Real d = u.norm();
    elongation = (Real)d;
    return elongation;
}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::addSpringForceWeight(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, int ivis, const Spring& spring)
{
    int a = spring.m1;
    Coord u = this->closestPos[i]-p[a];
    Real d = u.norm();
    if( d>1.0e-4 ) {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        double stiffweight = 1;
        helper::vector< bool > sourcevisible = l_meshprocessing->sourceVisible.getValue();

        if(sourcevisible[i]){
            int k = (int)closestpoint->closestSource[ivis].begin()->second;

            if(!closestpoint->targetIgnored[k]) {
                stiffweight = (double)l_rgbddataprocess->targetWeights.getValue()[k];
            } else {
                stiffweight = 1;
            }
        } else {
            stiffweight = 1;
        }

        //std::cout << " ok ok " << stiffweight << std::endl;

        //if (sourceborder[i]) stiffweight*=1;
        //double stiffweight = (double)1/targetWeights[(int)closestpoint->closestSource[i].begin()->second];

        potentialEnergy += stiffweight*elongation * elongation * spring.ks / 2;
        /*serr<<"addSpringForce, p = "<<p<<sendl;
        serr<<"addSpringForce, new potential energy = "<<potentialEnergy<<sendl;*/
        Deriv relativeVelocity = -v[a];
        Real elongationVelocity = dot(u,relativeVelocity);
        Real forceIntensity;

        if(theCloserTheStiffer.getValue()) {
            Real ks_max=stiffweight*spring.ks;
            Real ks_min=stiffweight*spring.ks/10;
            Real ks_mod = ks_min*(max-elongation)/(max-min)+ks_max*(elongation-min)/(max-min);
            forceIntensity = (Real)(ks_mod*elongation+spring.kd*elongationVelocity);
        } else {
            if (elongation < 0.02) {
                forceIntensity = (Real)(stiffweight*spring.ks*elongation+spring.kd*elongationVelocity);
            } else {
                forceIntensity = (Real)(stiffweight*spring.ks*elongation+spring.kd*elongationVelocity);
            }
        }
        Deriv force = u*forceIntensity;
        f[a]+=force;
        Mat& m = this->dfdx[i];
        Real tgt = forceIntensity * inverseLength;
        for( int j=0; j<N; ++j ) {
            // anisotropic
            //for( int k=0; k<N; ++k ) m[j][k] = tgt * u[j] * u[k];

            // isotropic
            for( int k=0; k<N; ++k ) {
                m[j][k] = ((Real)stiffweight*spring.ks-tgt) * u[j] * u[k];
            }
            m[j][j] += tgt;
        }
    } else {
    // null length, no force and no stiffness
        Mat& m = this->dfdx[i];
        for( int j=0; j<N; ++j ) {
            for( int k=0; k<N; ++k ) {
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
    if( d>1.0e-4 ) {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        potentialEnergy += coef*elongation * elongation *spring.ks / 2;

        Deriv relativeVelocity = -v[a];
        Real elongationVelocity = dot(u,relativeVelocity);
        Real forceIntensity;
        if(theCloserTheStiffer.getValue()) {
            Real ks_max=coef*spring.ks;
            Real ks_min=coef*spring.ks/10;
            Real ks_mod = ks_min*(max-elongation)/(max-min)+ks_max*(elongation-min)/(max-min);
            forceIntensity = (Real)(ks_mod*elongation+spring.kd*elongationVelocity);
        } else if (elongation < 0.02){
            forceIntensity = (Real)(coef*spring.ks*elongation+coef*spring.kd*elongationVelocity);
        } else {
            forceIntensity = (Real)(coef*spring.ks*elongation+coef*spring.kd*elongationVelocity);
        }

        Deriv force = u*forceIntensity;
        f[a]+=force;
        Mat& m = this->dfdx[i];
        Real tgt = forceIntensity * inverseLength;
                //std::cout << " u " << -v[a] << std::endl;
        for( int j=0; j<N; ++j ) {
            // anisotropic
            //for( int k=0; k<N; ++k ) m[j][k] = tgt * u[j] * u[k];

            // isotropic
            for( int k=0; k<N; ++k ) {
                m[j][k] = ((Real)coef*spring.ks-tgt) * u[j] * u[k];
            }
            m[j][j] += tgt;
        }
                        //dfdx1[i] = m;

    } else {
    // null length, no force and no stiffness
        Mat& m = this->dfdx[i];
        for( int j=0; j<N; ++j ) {
            /*for( int k=0; k<N; ++k )
            {
                m[j][k] = 0;
            }*/
        }
    }
}

template <class DataTypes>
void ClosestPointForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams,DataVecDeriv& _df , const DataVecDeriv&  _dx )
{

    VecDeriv& df = *_df.beginEdit();		//WDataRefVecDeriv df(_df);
    const VecDeriv&  dx = _dx.getValue();	// RDataRefVecDeriv dx(_dx);

    double kFactor 		 =  mparams->kFactor();
//    double bFactor       =  mparams->bFactor();

    if(ks.getValue()==0) {
        return;
    }

    const vector<Spring>& s = this->springs.getValue();

    for (unsigned int i=0; i<s.size(); i++) {
        const int a = s[i].m1;
        const Coord d = -dx[a];
        Deriv dforce = this->dfdx[i]*d;
        dforce *= kFactor;
        df[a]+=dforce;
    }
    //serr<<"addDForce, df = "<<f<<sendl;
    _df.endEdit();

}

template<class DataTypes>
void ClosestPointForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *m, SReal kFactor, unsigned int &offset)
{
    if(ks.getValue()==0) {
        return;
    }

    double kFact = kFactor;
    const vector<Spring >& ss = this->springs.getValue();
    const unsigned int n = (ss.size() < this->dfdx.size()) ? ss.size() : this->dfdx.size();
    for (unsigned int e=0; e<n; e++) {
        const Spring& s = ss[e];
        unsigned p1 = offset+Deriv::total_size*s.m1;
        const Mat& mt = this->dfdx[e];
        for(int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                Real k = (Real)(mt[i][j]*kFact);
                m->add(p1+i,p1+j, -k);
            }
        }
    }
}

// this function might be useless
template<class DataTypes>
void ClosestPointForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (! l_rgbddataprocess || ! l_meshprocessing || ! l_dataio) {
        std::cerr << "ClosestPtForceField : link for meshprocessing and/or rbgddataprocessing and/or dataio is not available" << std::endl ;
        return ;
    }
    std::vector< Vector3 > points;
    //const vector<Spring>& springs = this->springs.getValue();
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    points.resize(0);

    sofa::helper::vector< tri > triangles;
    triangles = l_meshprocessing->sourceTriangles.getValue();

//    std::cout << " XSIZE " << x.size() << std::endl;

    int t = (int)this->getContext()->getTime();

    if (l_rgbddataprocess->targetPositions.getValue().size()>0 && l_meshprocessing->sourceVisiblePositions.getValue().size()>0) {
        for (unsigned int i=0; i<x.size(); i++) {
            points.resize(0);
            Vector3 point = DataTypes::getCPos(x[i]);
            points.push_back(point);
        }
    }

    points.resize(0);

    Vector3 coefs;
    int index;
    long id;
    float xp0, yp0;
    int x0, y0;

//    sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
//    typename rgbdtracking::RGBDDataProcessing<DataTypes>::SPtr rgbddataprocessing;
//    root->get(rgbddataprocessing);

    cv::Mat distimg = l_rgbddataprocess->seg.distImage;
    cv::Mat foreg = l_rgbddataprocess->foreground;

    cv::Mat depth = l_dataio->depth ;
    if (useKLTPoints.getValue() && t >= startimageklt.getValue() && t%(niterations.getValue()) == 0){
        VecCoord targetKLTPos = targetKLTPositions.getValue();
        for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
           tracker.getFeature(k, id, xp0, yp0);
           x0 = (int)xp0;
           y0 = (int)yp0;
           coefs = mappingkltcoef[id];
           index = mappingkltind[id];
           float depthValue;
           int avalue;
           depthValue = (float)depth.at<float>(yp0,xp0);
           avalue = (int)distimg.at<uchar>(yp0,xp0);

           if( depthValue < 1 && depthValue > 0 && avalue > 2 && foreg.at<cv::Vec4b>(y0,x0)[3] > 0){
               Vector3 point1 = DataTypes::getCPos(targetKLTPos[id]);
               Vector3 point2 = (coefs[0]*DataTypes::getCPos(x[triangles[index][0]]) + coefs[1]*DataTypes::getCPos(x[triangles[index][1]]) + coefs[2]*DataTypes::getCPos(x[triangles[index][2]]));
               points.push_back(point1);
               points.push_back(point2);
           }

        }
    }
}

}

} // namespace sofa

//#endif  /* SOFA_COMPONENT_CLOSESTPOINTFORCEFIELD_INL */


