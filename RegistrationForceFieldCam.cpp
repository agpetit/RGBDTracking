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

#define SOFA_RGBDTRACKING_REGISTRATIONFORCEFIELDCAM_CPP

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

#include <pcl/common/common_headers.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/registration/icp.h>
#include <pcl/common/transforms.h>

#include "RegistrationForceFieldCam.h"
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

      SOFA_DECL_CLASS(RegistrationForceFieldCam)

      // Register in the Factory
      int RegistrationForceFieldCamClass = core::RegisterObject("Compute forces based on closest points from/to a target surface/point set")
    #ifndef SOFA_FLOAT
        .add< RegistrationForceFieldCam<Vec3dTypes> >()
    #endif
    #ifndef SOFA_DOUBLE
        .add< RegistrationForceFieldCam<Vec3fTypes> >()
    #endif
    ;

    #ifndef SOFA_FLOAT
      template class SOFA_RGBDTRACKING_API RegistrationForceFieldCam<Vec3dTypes>;
    #endif
    #ifndef SOFA_DOUBLE
      template class SOFA_RGBDTRACKING_API RegistrationForceFieldCam<Vec3fTypes>;
    #endif

using namespace helper;


template <class DataTypes>
RegistrationForceFieldCam<DataTypes>::RegistrationForceFieldCam(core::behavior::MechanicalState<DataTypes> *mm )
    : Inherit(mm)
    , ks(initData(&ks,(Real)0.0,"stiffness","uniform stiffness for the all springs."))
    , kd(initData(&kd,(Real)0.0,"damping","uniform damping for the all springs."))
    , blendingFactor(initData(&blendingFactor,(Real)1,"blendingFactor","blending between projection (=0) and attraction (=1) forces."))
    , projectToPlane(initData(&projectToPlane,false,"projectToPlane","project closest points in the plane defined by the normal."))
    , springs(initData(&springs,"spring","index, stiffness, damping"))
    , sourceSurfacePositions(initData(&sourceSurfacePositions,"sourceSurface","Points of the surface of the source mesh."))
    , sourcePositions(initData(&sourcePositions,"sourcePositions","Points of the mesh."))
    , sourceContourPositions(initData(&sourceContourPositions,"sourceContourPositions","Contour points of the surface of the mesh."))
    , sourceVisible(initData(&sourceVisible,"sourceVisible","Visibility of the points of the surface of the mesh."))
    , indicesVisible(initData(&indicesVisible,"indicesVisible","Indices of the visible points of the mesh."))
    , sourceVisiblePositions(initData(&sourceVisiblePositions,"sourceVisiblePositions","Visible points of the surface of the mesh."))
    , sourceBorder(initData(&sourceBorder,"sourceBorder","Points of the border of the mesh."))
    , sourceTriangles(initData(&sourceTriangles,"sourceTriangles","Triangles of the source mesh."))
    , sourceNormals(initData(&sourceNormals,"sourceNormals","Normals of the source mesh."))
    , sourceSurfaceNormals(initData(&sourceSurfaceNormals,"sourceSurfaceNormals","Normals of the surface of the source mesh."))
    , targetPositions(initData(&targetPositions,"targetPositions","Points of the target point cloud."))
    , targetContourPositions(initData(&targetContourPositions,"targetContourPositions","Contour points of the target point cloud."))
    , targetBorder(initData(&targetBorder,"targetBorder","Contour flag of the target point cloud."))
    , targetWeights(initData(&targetWeights,"targetWeights","Weights for the points of the target point cloud."))
    , drawSource(initData(&drawSource,false,"drawSource"," "))
    , drawTarget(initData(&drawTarget,false,"drawTarget"," "))
    , drawContour(initData(&drawContour,false,"drawContour"," "))
    , showArrowSize(initData(&showArrowSize,0.01f,"showArrowSize","size of the axis."))
    , drawMode(initData(&drawMode,0,"drawMode","The way springs will be drawn:\n- 0: Line\n- 1:Cylinder\n- 2: Arrow."))
    , drawColorMap(initData(&drawColorMap,true,"drawColorMap","Hue mapping of distances to closest point"))
    , theCloserTheStiffer(initData(&theCloserTheStiffer,false,"theCloserTheStiffer","Modify stiffness according to distance"))
    , useContour(initData(&useContour,false,"useContour","Emphasize forces close to the target contours"))
    , useVisible(initData(&useVisible,true,"useVisible","Use the vertices of the viisible surface of the source mesh"))
    , useRealData(initData(&useRealData,true,"useRealData","Use real data"))
    , useSensor(initData(&useSensor,false,"useSensor","Use the sensor"))
    , niterations(initData(&niterations,3,"niterations","Number of iterations in the tracking process"))
    , errorfunction(initData(&errorfunction,"errorfunction", "error"))
{
        iter_im = 0;

}

template <class DataTypes>
RegistrationForceFieldCam<DataTypes>::~RegistrationForceFieldCam()
{
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::reinit()
{

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue(); 	//RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
    this->clearSprings(x.size());

    for(unsigned int i=0;i<x.size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());

}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::init()
{

    this->Inherit::init();
    core::objectmodel::BaseContext* context = this->getContext();

        if(!(this->mstate)) this->mstate = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());

        // Get source normals
        if(!sourceNormals.getValue().size()) serr<<"normals of the source model not found"<<sendl;

    // add a spring for every input point
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue(); 			//RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
    this->clearSprings(x.size());

    npoints = x.size();

        for(unsigned int i=0;i<x.size();i++)
            this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());

    glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);

    sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());

    root->get(closestpoint);
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::resetSprings()
{

    this->clearSprings(sourceVisiblePositions.getValue().size());
        for(unsigned int i=0;i<sourceVisiblePositions.getValue().size();i++)
            this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());

}


template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::addForce(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{

    double timeaddforce = (double)getTickCount();
    addForceMesh(mparams, _f, _x, _v);
    std::cout << " TIME ADDFORCE " <<  (getTickCount() - timeaddforce)/getTickFrequency() << std::endl;
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::addForceMesh(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
    int t = (int)this->getContext()->getTime();


    sofa::helper::vector< tri > triangles;
    triangles = sourceTriangles.getValue();

    bool reinitv = false;

    helper::vector< bool > sourcevisible = sourceVisible.getValue();
    helper::vector< int > indicesvisible = indicesVisible.getValue();
    helper::vector< bool > sourceborder = sourceBorder.getValue();

                 if (t > 0 && t%niterations.getValue() == 0){

                if (npoints != (this->mstate->read(core::ConstVecCoordId::position())->getValue()).size())
                {
                        reinit();
                        reinitv = true;
                }

                npoints = (this->mstate->read(core::ConstVecCoordId::position())->getValue()).size();

        }

    double time = (double)getTickCount();

    if(ks.getValue()==0) return;

    double timef0 = (double)getTickCount();


    VecDeriv&        f = *_f.beginEdit();       //WDataRefVecDeriv f(_f);
    const VecCoord&  x = _x.getValue();			//RDataRefVecCoord x(_x);
    const VecDeriv&  v = _v.getValue();			//RDataRefVecDeriv v(_v);
    ReadAccessor< Data< VecCoord > > tn(targetNormals);
    ReadAccessor< Data< VecCoord > > tp(targetPositions);
    ReadAccessor< Data< VecCoord > > tcp(targetContourPositions);

                        if (t%niterations.getValue() == 0) {
                                f_.resize(f.size());
                                x_.resize(x.size());
                                v_.resize(v.size());
                        }

    const vector<Spring>& s = this->springs.getValue();
    this->dfdx.resize(s.size());
    this->closestPos.resize(s.size());

    dfdx1.resize(s.size());

    //closestpoint->updateClosestPointsGt();

        if (useVisible.getValue())
        closestpoint->sourceVisiblePositions.setValue(sourceVisiblePositions.getValue());

        closestpoint->timer = t;
        closestpoint->targetPositions.setValue(targetPositions.getValue());
        closestpoint->sourceSurfacePositions.setValue(sourceSurfacePositions.getValue());
        closestpoint->sourceBorder = sourceBorder.getValue();
        closestpoint->targetBorder = targetBorder.getValue();

    if (!useContour.getValue())
                closestpoint->updateClosestPoints();
        else
        {
                if (t<=2)
                closestpoint->updateClosestPoints();
                else
                {
                closestpoint->targetContourPositions.setValue(targetContourPositions.getValue());
                closestpoint->sourceContourPositions.setValue(sourceContourPositions.getValue());
                closestpoint->normalsContour = normalsContour;
                closestpoint->updateClosestPointsContours();
                }
        }

    double timeClosestPoint = ((double)getTickCount() - timef0)/getTickFrequency();

    std::cout << " TIME CLOSESTPOINT " << timeClosestPoint << std::endl;
        indices = closestpoint->getIndices();

    m_potentialEnergy = 0;

    // get attraction/ projection factors
    Real attrF=(Real) blendingFactor.getValue();
    if(attrF<(Real)0.) attrF=(Real)0.;
    if(attrF>(Real)1.) attrF=(Real)1.;
    Real projF=((Real)1.-attrF);

        //std::cout << " tp size " << tp.size() << std::endl;

        time = (double)getTickCount();


    if(tp.size()==0)
        for (unsigned int i=0; i<s.size(); i++)
            closestPos[i]=x[i];
    else {

        // count number of attractors
        cnt.resize(s.size()); cnt.fill(0);
                if(attrF>0)
                        if (!useVisible.getValue())
                        {
                                if (useContour.getValue() && t >= niterations.getValue())
                                {
                                        for (unsigned int i=0; i<tp.size(); i++)
                                        if(!closestpoint->targetIgnored[i])// && !rgbddataprocessing->targetBorder[i])
                                        cnt[closestpoint->closestTarget[i].begin()->second]++;
                                }
                                else
                                {
                                for (unsigned int i=0; i<tp.size(); i++)
                                if(!closestpoint->targetIgnored[i])
                                        cnt[closestpoint->closestTarget[i].begin()->second]++;
                                }
                        }
                        else
                        {
                        if (useContour.getValue()){
                                        if (t >2 )
                                        {
                                        for (unsigned int i=0; i<tp.size(); i++)
                                        {
                                                if(!closestpoint->targetIgnored[i])
                                                cnt[indicesvisible[closestpoint->closestTarget[i].begin()->second]]++;
                                        }
                                        }
                                        else
                                        for (unsigned int i=0; i<tp.size(); i++) cnt[closestpoint->closestTarget[i].begin()->second]++;

                                }
                                else
                                {
                                if (t > 2 )
                                        {
                                                for (unsigned int i=0; i<tp.size(); i++)
                                                {

                                                //std::cout << " ind " << indicesvisible[closestpoint->closestTarget[i].begin()->second] << " " << closestpoint->closestTarget[i].begin()->second << std::endl;
                                                        if(!closestpoint->targetIgnored[i])// && !targetBackground[i])
                                                        cnt[indicesvisible[closestpoint->closestTarget[i].begin()->second]]++;
                                                }
                                        }
                                        else
                                        {
                                        for (unsigned int i=0; i<tp.size(); i++) cnt[closestpoint->closestTarget[i].begin()->second]++;
                                        }
                                }
                        }

                std::cout << " tp size0 " << tp.size() << std::endl;

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

                double error = 0;
                int nerror = 0;

                int ivis=0;
                int kk = 0;
                unsigned int id;
                        {
        if(projF>0) {
                if (!useVisible.getValue())
                {


                        if (useContour.getValue() && t > niterations.getValue() )//&& t%niterations.getValue() == 0)
                        {
                                for (unsigned int i=0; i<s.size(); i++)
                                {
                                unsigned int id=closestpoint->closestSource[i].begin()->second;
                                        if(!closestpoint->sourceIgnored[i])
                                        {
                                                if(!sourceborder[i])
                                                {
                                                id=closestpoint->closestSource[i].begin()->second;
                                                if(projectToPlane.getValue() && tn.size()!=0)	closestPos[i]=/*(1-(Real)sourceWeights[i])**/(x[i]+tn[id]*dot(tp[id]-x[i],tn[id]))*projF;
                                                else closestPos[i]=/*(1-(Real)sourceWeights[i])**/tp[id]*projF;
                                        /*id=indices[kk];
                                        closestPos[i]+=(Real)sourceWeights[i]*tcp[id]*projF;
                                        kk++;*/

                    if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                        //closestPos[i] = x[i];//*projF;

                                        }
                                        else {
                                                id=indices[kk];
                                                closestPos[i]=tcp[id]*projF;

                                                if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                                kk++;
                                                }

                                                }
                                        else {
                    closestPos[i]=x[i]*projF;
                    if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                        }
                                                }
                                //std::cout << " tp size2 " << tp.size() << std::endl;

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

                                    if (t >2 ){
                                if (useContour.getValue())
                                                {
                                for (unsigned int i=0; i<s.size(); i++)
                                {
               // if(/*!closestpoint->sourceIgnored[i] &&*/ sourcevisible[i])
                                        {
                                        if (sourcevisible[i]){
                                        if(sourceborder[i])
                                                {
                                        id=closestpoint->closestSource[i].begin()->second;
                    if(projectToPlane.getValue() && tn.size()!=0)	closestPos[i]=/*(1-(Real)sourceWeights[i])**/(x[i]+tn[id]*dot(tp[id]-x[i],tn[id]))*projF;
                    else closestPos[i]=/*(1-(Real)sourceWeights[i])**/tp[id]*projF;

                                        /*id=indices[kk];
                                        closestPos[i]+=(Real)sourceWeights[i]*tcp[id]*projF;
                                        kk++;*/

                    if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                        //closestPos[i] = x[i];//*projF;
                                        //if(!cnt[i]) closestPos[i]+=x[i]*attrF;

                                        }
                                        else {

                                                unsigned int id=indices[kk];
                        //std::cout << " tp size2 " << tcp.size() << " " << id << " " << kk << " " << indices.size()<< std::endl;

                                                closestPos[i]=tcp[id]*projF;

                                                if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                        kk++;
                        //std::cout << " tp size20 " << id << std::endl;

                                                }
                                                }
                                                else
                                                {
                                        //closestPos[i]=x[i];
                                        closestPos[i]=x[i]*projF;
                                        if(!cnt[i]) {closestPos[i]+=x[i]*attrF;}
                                                }

                                        }
                                        }

                                                }
                                                else{
                                        //std::cout << " tp size10 " << tp.size() << std::endl;
                                                for (unsigned int i=0; i<s.size(); i++)
                                                {
               // if(/*!closestpoint->sourceIgnored[i] &&*/ sourcevisible[i])
                                                    //std::cout << " source visible " << (int)sourcevisible[i] << std::endl;

                                        {

                                        if (sourcevisible[i]){
                                                unsigned int id=closestpoint->closestSource[ivis].begin()->second;
                                                if(!closestpoint->sourceIgnored[ivis])
                                                        {
                                                closestPos[i]=tp[id]*projF;

                                   error += this->computeError(x[i],tp[id]);
                                                nerror++;}
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

                                        //std::cout << " tp size11 " << tp.size() << std::endl;

                                        if(!cnt[i]) {closestPos[i]+=x[i]*attrF;}


                                        }
                                        }
                                                        }

                                        }
                                        else
                                        {
                                        std::cout << " tp size12 " << s.size() << std::endl;

                                for (unsigned int i=0; i<s.size(); i++)
                                {
               // if(/*!closestpoint->sourceIgnored[i] &&*/ sourcevisible[i])
                                        {

                                        //unsigned int id=closestpoint->closestSource[i].begin()->second;

                                        //std::cout << " tp size12 " << id << std::endl;

                    if(projectToPlane.getValue() && tn.size()!=0)	closestPos[i]=(x[i]+tn[id]*dot(tp[id]-x[i],tn[id]))*projF;
                    else closestPos[i]=tp[id]*projF;

                    if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                                        }
                                        }
                                        //std::cout << " tp size13 " << tp.size() << std::endl;

                                        }
                        }
        }
        else for (unsigned int i=0; i<s.size(); i++) { if(!cnt[i]) closestPos[i]=x[i]; else closestPos[i].fill(0); }

                // attraction

        if(attrF>0)
                        if(!useVisible.getValue())
                        {
            for (unsigned int i=0; i<tp.size(); i++){
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
                                else
                                {
                                        if (useContour.getValue() && t > niterations.getValue() )//&& t%niterations.getValue() > 0)
                                        {
                                        //if (sourceborder[id] && rgbddataprocessing->targetBorder[i])
                                        closestPos[id]+=tp[i]*attrF/(Real)cnt[id];
                                        }


                                }
                }

                        }
                        else
                        {
                                        if (t >2){
                                if( !useContour.getValue()) //&& !targetBackground[i])
                                {
                                int kkt = 0;
                                for (unsigned int i=0; i<tp.size(); i++)
                                        {
                                unsigned int id=closestpoint->closestTarget[i].begin()->second;

                if(!closestpoint->targetIgnored[i]) //&& !targetBackground[i])
                                {
                                        //std::cout << " id target " << i << " id source " << id << std::endl;

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
                                else
                                {

                                for (unsigned int i=0; i<tp.size(); i++)
                                        {
                                unsigned int id=closestpoint->closestTarget[i].begin()->second;
                if(!closestpoint->targetIgnored[i]) //&& !targetBackground[i])
                                {
                                        unsigned int id1;
                    //if (!rgbddataprocessing->targetBorder[i])
                                        id1 = indicesvisible[id];
                                        /*else {
                                        id1 = indicesTarget[kkt];
                                        kkt++;
                                        }*/
                                        //if (sourceborder[id1] && rgbddataprocessing->targetBorder[i])
                                        closestPos[id1]+=tp[i]*attrF/(Real)cnt[id1];
                                }
                                                                //if(rgbddataprocessing->targetBorder[i])
                                                                        {
                                        }

                                        }
                                        }
                                        }
                                        else
                                        {
                                int kkt = 0;
                                for (unsigned int i=0; i<tp.size(); i++)
                                        {

                                unsigned int id=closestpoint->closestTarget[i].begin()->second;

                if(!closestpoint->targetIgnored[i]) //&& !targetBackground[i])
                                        closestPos[id]+=tp[i]*attrF/(Real)cnt[id];
                                        }

                                        }
                        }



        }
    }


        ind = 0;
        int kc = 0;

    for (unsigned int i=0; i<s.size(); i++)
    {
        //serr<<"addForce() between "<<springs[i].m1<<" and "<<closestPos[springs[i].m1]<<sendl;
                                if (t > 2*niterations.getValue() && t%(niterations.getValue()) == 0)
                                        {

                                        if( !useContour.getValue())
                                        this->addSpringForce(m_potentialEnergy,f,x,v, i, s[i]);
                                        else this->addSpringForceWeight(m_potentialEnergy,f,x,v, i, s[i]);


                                        }
    }

    _f.endEdit();


}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
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
        /*          serr<<"addSpringForce, p = "<<p<<sendl;
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
void RegistrationForceFieldCam<DataTypes>::addStoredSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
{
    int a = spring.m1;

        f[a]=f_[a];
        Mat& m = this->dfdx[i];
                m = dfdx1[i];
}

template <class DataTypes>
double RegistrationForceFieldCam<DataTypes>::computeError(Vector3 sourcePoint, Vector3 targetPoint)
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
void RegistrationForceFieldCam<DataTypes>::addSpringForceWeight(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
{
    int a = spring.m1;
    Coord u = this->closestPos[i]-p[a];
    Real d = u.norm();
    if( d>1.0e-4 )
    {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        double stiffweight;
        helper::vector< bool > sourceborder = sourceBorder.getValue();


                for (int k = 0; k < targetPositions.getValue().size(); k++){
                        if(closestpoint->closestTarget[k].begin()->second == i || closestpoint->closestSource[i].begin()->second == k)
                        {
                                if(sourceborder[i])
                        stiffweight = (double)sourceWeights[i];
                            else stiffweight = (double)sourceWeights[i];
                        //stiffweight = (double)targetWeights[k];
                        }
                        /*else if (closestpoint->closestSource[i].begin()->second == k){
                        //stiffweight = (double)combinedWeights[ind];
                        //stiffweight = (double)sourceWeights[i]*targetWeights[k];
                        stiffweight = (double)targetWeights[k];
                        }*/
                        ind++;
                        }

                if (sourceborder[i]) stiffweight*=1;
                        //double stiffweight = (double)1/targetWeights[(int)closestpoint->closestSource[i].begin()->second];

        potentialEnergy += stiffweight*elongation * elongation * spring.ks / 2;
        /*          serr<<"addSpringForce, p = "<<p<<sendl;
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

template<class DataTypes>
void RegistrationForceFieldCam<DataTypes>::addSpringDForce(VecDeriv& df,const  VecDeriv& dx, int i, const Spring& spring, double kFactor, double /*bFactor*/)
{
    const int a = spring.m1;
    const Coord d = -dx[a];
    Deriv dforce = this->dfdx[i]*d;
    dforce *= kFactor;
    df[a]+=dforce;
    //serr<<"addSpringDForce, a="<<a<<", b="<<b<<", dforce ="<<dforce<<sendl;
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::addDForce(const core::MechanicalParams* mparams,DataVecDeriv& _df , const DataVecDeriv&  _dx )
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
void RegistrationForceFieldCam<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *m, SReal kFactor, unsigned int &offset)
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
void RegistrationForceFieldCam<DataTypes>::draw(const core::visual::VisualParams* vparams)
{



}

}
}
} // namespace sofa

//#endif  /* SOFA_COMPONENT_INTERACTIONFORCEFIELD_RegistrationForceFieldCam_INL */


