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

#define SOFA_RGBDTRACKING_REGISTRATIONRIGID_CPP

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/Mapping.inl>
#include <sofa/simulation/Simulation.h>
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
#include <pcl/io/pcd_io.h>
#include <pcl/registration/icp.h>
#include <pcl/common/transforms.h>
#include <pcl/registration/icp.h>

#include "RegistrationRigid.h"
#include "ImageConverter.h"

using std::cerr;
using std::endl;

namespace sofa
{

namespace core
{


    using namespace sofa::defaulttype;

      SOFA_DECL_CLASS(RegistrationRigid)

      // Register in the Factory
      int RegistrationRigidClass = core::RegisterObject("Compute forces based on closest points from/to a target surface/point set")
    #ifndef SOFA_FLOAT
        .add< RegistrationRigid<Vec3dTypes> >()
    #endif
    #ifndef SOFA_DOUBLE
        .add< RegistrationRigid<Vec3fTypes> >()
    #endif
    ;

    #ifndef SOFA_FLOAT
      template class SOFA_RGBDTRACKING_API RegistrationRigid<Vec3dTypes>;
    #endif
    #ifndef SOFA_DOUBLE
      template class SOFA_RGBDTRACKING_API RegistrationRigid<Vec3fTypes>;
    #endif

using namespace helper;


template <class DataTypes>
RegistrationRigid<DataTypes>::RegistrationRigid()
    : Inherit()
	, sourceSurfacePositions(initData(&sourceSurfacePositions,"sourceSurface","Points of the surface of the source mesh."))
	, sourcePositions(initData(&sourcePositions,"sourcePositions","Points of the mesh."))
        , sourceVisiblePositions(initData(&sourceVisiblePositions,"sourceVisiblePositions","Visible points of the surface of the mesh."))
        , targetPositions(initData(&targetPositions,"targetPositions","Points of the surface of the source mesh."))
        , sourceTriangles(initData(&sourceTriangles,"sourceTriangles","Triangles of the source mesh."))
        , sourceNormals(initData(&sourceNormals,"sourceNormals","Normals of the source mesh."))
	, sourceSurfaceNormals(initData(&sourceSurfaceNormals,"sourceSurfaceNormals","Normals of the surface of the source mesh."))
	, useVisible(initData(&useVisible,true,"useVisible","Use the vertices of the viisible surface of the source mesh"))
        , forceRegistration(initData(&forceRegistration,true,"forceRegistration","soft registration through ICP based forces"))
	//, useKalman(initData(&useKalman,false,"useKalman","Use the Kalman filter"))
        ,niterations(initData(&niterations,3,"niterations","Number of iterations in the tracking process"))
        ,startimage(initData(&startimage,1,"startimage","Frame index to start rigid registration"))
	,translation(initData(&translation,"translation", "translation parameters"))
	,rotation(initData(&rotation,"rotation", "rotation parameters"))
        ,rigidForces(initData(&rigidForces,"rigidforces", "rigid forces"))
        ,rigidState(initData(&rigidState,"rigidstate", "rigid state"))
        ,stopAfter(initData(&stopAfter,300000,"stopafter", "rigid state"))
        ,MeshToPointCloud(initData(&MeshToPointCloud,true,"meshToPointCloud", "rigid state"))
{
	this->f_listening.setValue(true); 

	iter_im = 0;
}

template <class DataTypes>
RegistrationRigid<DataTypes>::~RegistrationRigid()
{
}

template <class DataTypes>
void RegistrationRigid<DataTypes>::reinit()
{

}

template <class DataTypes>
void RegistrationRigid<DataTypes>::init()
{
    this->Inherit::init();
    core::objectmodel::BaseContext* context = this->getContext();

    mstate = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());
    // Get source normals
    if(!sourceNormals.getValue().size()) serr<<"normals of the source model not found"<<sendl;
	
    // add a spring for every input point
    const VecCoord& x = mstate->read(core::ConstVecCoordId::position())->getValue(); 			//RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));

    rigidForces.setValue(x);
	
    sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
    root->get(rgbddataprocessing);
}

template <class DataTypes>
void RegistrationRigid<DataTypes>::determineRigidTransformation ()
{
	
    const VecCoord& x = mstate->read(core::ConstVecCoordId::position())->getValue();
    const VecCoord&  tp = targetPositions.getValue();
	
    unsigned int nbs=x.size(),nbt=tp.size();
	
    source.reset(new pcl::PointCloud<pcl::PointXYZ>);
    source_registered.reset(new pcl::PointCloud<pcl::PointXYZ>);
	
    pcl::PointXYZ newPoint;
        for (unsigned int i=0; i < nbs; i++)
        {
            newPoint.z = x[i][2];
            newPoint.x = x[i][0];
            newPoint.y = x[i][1];
            /*newPoint.r = 0;
            newPoint.g = 0;
            newPoint.b = 0;*/
            source->points.push_back(newPoint);
            //std::cout << " x source  " << x[i][0] << " " << x[i][1] << " " << x[i][2] << std::endl;
        }

    target.reset(new pcl::PointCloud<pcl::PointXYZ>);

        for (unsigned int i=0; i < nbt; i++)
        {
            newPoint.z = tp[i][2];
            newPoint.x = tp[i][0];
            newPoint.y = tp[i][1];
            /*newPoint.r = 0;
            newPoint.g = 0;
            newPoint.b = 0;*/
            target->points.push_back(newPoint);
            //std::cout << " target  " << tp[i][0] << " " << tp[i][1] << " " << tp[i][2] << std::endl;
        }
	
    pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> registration;
    //registration.setInputCloud(rgbddataprocessing->target);
    registration.setInputSource(target);
    registration.setInputTarget(source);
    //registration->setInputCloud(source_segmented_);
    registration.setMaxCorrespondenceDistance(0.10);
    registration.setTransformationEpsilon (0.000001);
    registration.setMaximumIterations (1000);

  // Register
  registration.align (*source_registered);

  Eigen::Matrix4f transformation_matrix1 = registration.getFinalTransformation();
  transformation_matrix = transformation_matrix1.inverse();

  std::cout << " rigid registration pcd size " << transformation_matrix << std::endl;
  
Eigen::Matrix3f mat;
  for (int i = 0; i < 3; i++)
	  for (int j = 0; j< 3;j ++)
	  {
       mat(i,j) = transformation_matrix(i,j);
	  }
	  
vpHomogeneousMatrix cMoin;
for (int i = 0; i < 4; i++)
	  for (int j = 0; j< 4;j ++)
	  cMoin[i][j] = transformation_matrix(i,j);
	  
cMo = cMoin*(vpExponentialMap::direct(0*kalman.predictedVelPose).inverse());
	   	   
Eigen::Vector3f ea = mat.eulerAngles(0,1,2);

Eigen::Quaternionf q(mat);

VecReal trans;
VecReal rot;
trans.resize(3);
rot.resize(9);
trans[0] = (double)transformation_matrix(0,3);
trans[1] = (double)transformation_matrix(1,3);
trans[2] = (double)transformation_matrix(2,3);

rot[0] = transformation_matrix(0,0);
rot[1] = transformation_matrix(0,1);
rot[2] = transformation_matrix(0,2);

rot[3] = transformation_matrix(1,0);
rot[4] = transformation_matrix(1,1);
rot[5] = transformation_matrix(1,2);

rot[6] = transformation_matrix(2,0);
rot[7] = transformation_matrix(2,1);
rot[8] = transformation_matrix(2,2);

translation.setValue(trans);
rotation.setValue(rot);

pcl::transformPointCloud(*source, *source_registered, transformation_matrix);

helper::WriteAccessor< Data<VecCoord> > x1 = *mstate->write(core::VecCoordId::position());
VecCoord xrigid;
xrigid.resize(nbs);
double stiffness = 0.05;
double normerror = 0;

if (forceRegistration.getValue())
{
for (unsigned int i=0; i<nbs; i++)
{
    newPoint = source_registered->points[i];
    /*xrigid[i][0] = stiffness*(newPoint.x - x[i][0]);
    xrigid[i][1] = stiffness*(newPoint.y - x[i][1]);
    xrigid[i][2] = stiffness*(newPoint.z - x[i][2]);*/
    xrigid[i][0] = x1[i][0];
    xrigid[i][1] = x1[i][1];
    xrigid[i][2] = x1[i][2];
    normerror+=xrigid[i].norm();
}
//std::cout << " normerror " << normerror << std::endl;
}
else
{
for (unsigned int i=0; i<nbs; i++)
{
    newPoint = source_registered->points[i];
    x1[i][0] = newPoint.x;
    x1[i][1] = newPoint.y;
    x1[i][2] = newPoint.z;
}
}
rigidForces.setValue(xrigid);

}

template <class DataTypes>
void RegistrationRigid<DataTypes>::determineRigidTransformationVisible ()
{
    const VecCoord& x0 = mstate->read(core::ConstVecCoordId::position())->getValue();
    const VecCoord& x = sourceVisiblePositions.getValue();
    const VecCoord&  tp = targetPositions.getValue();
	
    unsigned int nbs=x.size(),nbt=tp.size(), nbs0 = x0.size();

    source.reset(new pcl::PointCloud<pcl::PointXYZ>);
    source_registered.reset(new pcl::PointCloud<pcl::PointXYZ>);
	
        pcl::PointXYZ newPoint;
	for (unsigned int i=0; i<nbs; i++)
	{
	newPoint.z = x[i][2];
	newPoint.x = x[i][0];
	newPoint.y = x[i][1];
        /*newPoint.r = 0;
	newPoint.g = 0;
        newPoint.b = 0;*/
	source->points.push_back(newPoint);
        //std::cout << " x source  " << x[i][0] << " " << x[i][1] << " " << x[i][2] << std::endl;

	}

   source0.reset(new pcl::PointCloud<pcl::PointXYZ>);
   source_registered0.reset(new pcl::PointCloud<pcl::PointXYZ>);
 	
	for (unsigned int i=0; i<nbs0; i++)
	{
	newPoint.z = x0[i][2];
	newPoint.x = x0[i][0];
	newPoint.y = x0[i][1];
        /*newPoint.r = 0;
	newPoint.g = 0;
        newPoint.b = 0;*/
	source0->points.push_back(newPoint);
        }

    target.reset(new pcl::PointCloud<pcl::PointXYZ>);

        for (unsigned int i=0; i < nbt; i++)
        {
        newPoint.z = tp[i][2];
        newPoint.x = tp[i][0];
        newPoint.y = tp[i][1];
        /*newPoint.r = 0;
        newPoint.g = 0;
        newPoint.b = 0;*/
        target->points.push_back(newPoint);
        //std::cout << " x source  " << x[i][0] << " " << x[i][1] << " " << x[i][2] << std::endl;
        }


  pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> registration;
  if (MeshToPointCloud.getValue()){
    //registration.setInputCloud(rgbddataprocessing->target);
  registration.setInputSource(target);
  //registration->setInputCloud(source_segmented_);
  registration.setInputTarget (source);}
  else
  {
      //registration.setInputCloud(rgbddataprocessing->target);
      registration.setInputSource(source);
      //registration->setInputCloud(source_segmented_);
      registration.setInputTarget (target);
  }
  registration.setMaxCorrespondenceDistance(0.05);
  //registration.setMaxCorrespondenceDistance(0.04);
  registration.setTransformationEpsilon (0.000001);
  registration.setMaximumIterations (100);

  // Register
  registration.align (*source_registered);

  Eigen::Matrix4f transformation_matrix1 = registration.getFinalTransformation();
  if (MeshToPointCloud.getValue())
  transformation_matrix = transformation_matrix1.inverse();
  else transformation_matrix = transformation_matrix1;
  
  Eigen::Matrix3f mat;
  for (int i = 0; i < 3; i++)
	  for (int j = 0; j< 3;j ++)
       mat(i,j) = transformation_matrix(i,j);
	   	   
Eigen::Vector3f ea = mat.eulerAngles(0,1,2);
Eigen::Quaternionf q(mat);

VecReal trans;
VecReal rot;
trans.resize(3);
rot.resize(3);

trans[0] = (double)transformation_matrix(0,3);
trans[1] = (double)transformation_matrix(1,3);
trans[2] = (double)transformation_matrix(2,3);

rot[0] = ea(0);
rot[1] = ea(1);
rot[2] = ea(2);

translation.setValue(trans);
rotation.setValue(rot);

pcl::transformPointCloud(*source0, *source_registered0, transformation_matrix);

helper::WriteAccessor< Data<VecCoord> > x1 = *mstate->write(core::VecCoordId::position());
helper::WriteAccessor< Data<VecCoord> > v1 = *mstate->write(core::VecDerivId::velocity());
helper::WriteAccessor< Data<VecCoord> > f1 = *mstate->write(core::VecDerivId::force());

//mstateRigid->setName("rigidstate");
//static_cast< core::objectmodel::BaseObject *>(mstateRigid) == static_cast< const core::objectmodel::BaseObject *>(this);
//*mstateRigid = *dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());
//helper::WriteAccessor< Data<VecCoord> > x1r = *mstateRigid->write(core::VecCoordId::position());

VecCoord xrigid;
xrigid.resize(nbs0);

if (forceRegistration.getValue())
{
for (unsigned int i=0; i<nbs0; i++)
{
    newPoint = source_registered0->points[i];

    xrigid[i][0] = newPoint.x;
    xrigid[i][1] = newPoint.y;
    xrigid[i][2] = newPoint.z;

    /*xrigid[i][0] = stiffness*(newPoint.x - x[i][0]);
    xrigid[i][1] = stiffness*(newPoint.y - x[i][1]);
    xrigid[i][2] = stiffness*(newPoint.z - x[i][2]);*/
    /*v1[i][0] = stiffness*(newPoint.x - x[i][0]);
    v1[i][1] = stiffness*(newPoint.y - x[i][1]);
    v1[i][2] = stiffness*(newPoint.z - x[i][2]);*/
    /*f1[i][0] += stiffness*(newPoint.x - x[i][0]);
    f1[i][1] += stiffness*(newPoint.y - x[i][1]);
    f1[i][2] += stiffness*(newPoint.z - x[i][2]);*/
    //normerror+=xrigid[i].norm();
}
}
else
{
for (unsigned int i=0; i<nbs0; i++)
{
    newPoint = source_registered0->points[i];
    x1[i][0] = newPoint.x;
    x1[i][1] = newPoint.y;
    x1[i][2] = newPoint.z;
}
}
rigidForces.setValue(xrigid);


//rigidState.setValue(mstateRigid->getName());

  	//for (unsigned int i=0; i<nbs; i++)
	{
		//mstate->applyTranslation(transformation_matrix(0,3),transformation_matrix(1,3),transformation_matrix(2,3));
		
		//mstate->applyTranslation(-0.01,0.0,0.0);
		//mstate->applyRotation(ea(0),ea(1),ea(2));

        /*(mstate->read(core::ConstVecCoordId::position())->getValue())[i][0] = 0;
        (mstate->read(core::ConstVecCoordId::position())->getValue())[i][1] = 0;
        (mstate->read(core::ConstVecCoordId::position())->getValue())[i][2] = 0;*/

	}
 // cout << "OK" <<  transformation_matrix(0,3) << " " << transformation_matrix(1,3) << " " << transformation_matrix(2,3) << endl;
}

template <class DataTypes>
double RegistrationRigid<DataTypes>::determineErrorICP ()
{
	
    const VecCoord& x = sourceSurfacePositions.getValue();
    const VecCoord&  tp = targetPositions.getValue();
	
    unsigned int nbs=x.size(),nbt=tp.size();
	
    sourceSurfacePointCloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
    sourceSurfacePointCloud_registered.reset(new pcl::PointCloud<pcl::PointXYZ>);
	
    pcl::PointXYZ newPoint;
	for (unsigned int i=0; i < nbs; i++)
	{
            newPoint.z = x[i][2];
            newPoint.x = x[i][0];
            newPoint.y = x[i][1];
            /*newPoint.r = 0;
            newPoint.g = 0;
            newPoint.b = 0;*/
            sourceSurfacePointCloud->points.push_back(newPoint);
            //std::cout << "  " << x[i][0] << " " << x[i][1] << " " << x[i][2] << std::endl;
	} 

        target.reset(new pcl::PointCloud<pcl::PointXYZ>);

            for (unsigned int i=0; i < nbt; i++)
            {
            newPoint.z = tp[i][2];
            newPoint.x = tp[i][0];
            newPoint.y = tp[i][1];
            /*newPoint.r = 0;
            newPoint.g = 0;
            newPoint.b = 0;*/
            target->points.push_back(newPoint);
            //std::cout << " x source  " << x[i][0] << " " << x[i][1] << " " << x[i][2] << std::endl;
            }
	
    cout << "final registration..." << std::flush;
    pcl::Registration<pcl::PointXYZ, pcl::PointXYZ>::Ptr registration1 (new pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ>);
    //registration1->setInputCloud(rgbddataprocessing->targetPointCloud);
    registration1->setInputSource(target);
    //registration->setInputCloud(source_segmented_);
    registration1->setInputTarget (sourceSurfacePointCloud);
    registration1->setMaxCorrespondenceDistance(0.10);
    registration1->setRANSACOutlierRejectionThreshold (0.1);
    registration1->setTransformationEpsilon (0.000001);
    registration1->setMaximumIterations (100);

    registration1->align(*sourceSurfacePointCloud_registered);

    double fitnessscore;
    fitnessscore = registration1->getFitnessScore(1000);
        return fitnessscore;

}

template <class DataTypes>
void RegistrationRigid<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{
    if (dynamic_cast<simulation::AnimateBeginEvent*>(event)) RegisterRigid();
}


template <class DataTypes>
void RegistrationRigid<DataTypes>::RegisterRigid()
{

	int t = (int)this->getContext()->getTime();
		
                 if (t >= startimage.getValue() && t%niterations.getValue() == 0){
			
		double time1 = (double)getTickCount();
		
		if (!useVisible.getValue()) determineRigidTransformation();
		else 
		{

		    if (t < stopAfter.getValue()){
                    if (t < 10)
                    determineRigidTransformation();
                    else determineRigidTransformationVisible();}
				
                }
				
			
			time1 = ((double)getTickCount() - time1)/getTickFrequency();
                        cout << "TIME RIGID ICP " << time1 << endl;
	
	}
	
		
}
            

}
} // namespace sofa



