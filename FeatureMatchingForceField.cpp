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

#define SOFA_RGBDTRACKING_FeatureMatchingForceField_CPP

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

#include <pcl/io/pcd_io.h>
#include <pcl/registration/transforms.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/keypoints/sift_keypoint.h>
#include <pcl/keypoints/harris_3d.h>

#include <pcl/features/fpfh_omp.h>
#include <pcl/features/pfh.h>
#include <pcl/features/pfhrgb.h>
#include <pcl/features/3dsc.h>
#include <pcl/features/shot_omp.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/kdtree/impl/kdtree_flann.hpp>

#include <pcl/registration/correspondence_estimation.h>
#include <pcl/registration/correspondence_rejection.h>
#include <pcl/registration/correspondence_rejection_sample_consensus.h>

#ifdef Success
  #undef Success
#endif

#include "FeatureMatchingForceField.h"


using std::cerr;
using std::endl;

namespace sofa
{

namespace component
{

namespace forcefield
{

    using namespace sofa::defaulttype;

      SOFA_DECL_CLASS(FeatureMatchingForceField)

      // Register in the Factory
      int FeatureMatchingForceFieldClass = core::RegisterObject("Compute forces based on closest points from/to a target surface/point set")
    #ifndef SOFA_FLOAT
        .add< FeatureMatchingForceField<Vec3dTypes> >()
    #endif
    #ifndef SOFA_DOUBLE
        .add< FeatureMatchingForceField<Vec3fTypes> >()
    #endif
    ;

    #ifndef SOFA_FLOAT
      template class SOFA_RGBDTRACKING_API FeatureMatchingForceField<Vec3dTypes>;
    #endif
    #ifndef SOFA_DOUBLE
      template class SOFA_RGBDTRACKING_API FeatureMatchingForceField<Vec3fTypes>;
    #endif

using namespace helper;


template <class DataTypes>
FeatureMatchingForceField<DataTypes>::FeatureMatchingForceField(core::behavior::MechanicalState<DataTypes> *mm )
    : Inherit(mm)
    , ks(initData(&ks,(Real)0.0,"stiffness","uniform stiffness for the all springs."))
    , kd(initData(&kd,(Real)0.0,"damping","uniform damping for the all springs."))
    , blendingFactor(initData(&blendingFactor,(Real)1,"blendingFactor","blending between projection (=0) and attraction (=1) forces."))
    , projectToPlane(initData(&projectToPlane,false,"projectToPlane","project closest points in the plane defined by the normal."))
    , springs(initData(&springs,"spring","index, stiffness, damping"))
    , cameraIntrinsicParameters(initData(&cameraIntrinsicParameters,Vector4(),"cameraIntrinsicParameters","camera parameters"))
    , sourceSurfacePositions(initData(&sourceSurfacePositions,"sourceSurface","Points of the surface of the source mesh."))
    , sourcePositions(initData(&sourcePositions,"sourcePositions","Points of the mesh."))
    , sourceVisible(initData(&sourceVisible,"sourceVisible","Visibility of the points of the surface of the mesh."))
    , indicesVisible(initData(&indicesVisible,"indicesVisible","Indices of the visible points of the mesh."))
    , sourceVisiblePositions(initData(&sourceVisiblePositions,"sourceVisiblePositions","Visible points of the surface of the mesh."))
    , sourceBorder(initData(&sourceBorder,"sourceBorder","Points of the border of the mesh."))
    , sourceTriangles(initData(&sourceTriangles,"sourceTriangles","Triangles of the source mesh."))
    , sourceNormals(initData(&sourceNormals,"sourceNormals","Normals of the source mesh."))
    , sourceSurfaceNormals(initData(&sourceSurfaceNormals,"sourceSurfaceNormals","Normals of the surface of the source mesh."))
    , targetPositions(initData(&targetPositions,"targetPositions","Points of the target point cloud."))
    , descriptor_type(initData(&descriptor_type,0,"descriptor","Descriptor type"))
    , keypoint_type(initData(&keypoint_type,0,"keypoint","Keypoint type"))
    , drawMode(initData(&drawMode,0,"drawMode","The way springs will be drawn:\n- 0: Line\n- 1:Cylinder\n- 2: Arrow."))
    , outlierThreshold(initData(&outlierThreshold,(Real)7,"outlierThreshold","suppress outliers when distance > (meandistance + threshold*stddev)."))
    , rejectBorders(initData(&rejectBorders,false,"rejectBorders","ignore border vertices."))
    , showArrowSize(initData(&showArrowSize,0.01f,"showArrowSize","size of the axis."))
    , drawColorMap(initData(&drawColorMap,true,"drawColorMap","Hue mapping of distances to closest point"))
    , theCloserTheStiffer(initData(&theCloserTheStiffer,false,"theCloserTheStiffer","Modify stiffness according to distance"))
{
    iter_im = 0;
}

template <class DataTypes>
FeatureMatchingForceField<DataTypes>::~FeatureMatchingForceField()
{
}

template <class DataTypes>
void FeatureMatchingForceField<DataTypes>::reinit()
{

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue(); 	//RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
    this->clearSprings(x.size());

    for(unsigned int i=0;i<x.size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());

}

template <class DataTypes>
void FeatureMatchingForceField<DataTypes>::init()
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

}

template <class DataTypes>
void FeatureMatchingForceField<DataTypes>::resetSprings()
{
    this->clearSprings(sourceVisiblePositions.getValue().size());
        for(unsigned int i=0;i<sourceVisiblePositions.getValue().size();i++)
            this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());
}


template <class DataTypes>
void FeatureMatchingForceField<DataTypes>::addForce(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{

    double timeaddforce = (double)getTickCount();
    int t = (int)this->getContext()->getTime();
    if (t > 5)
    addForceMesh(mparams, _f, _x, _v);
    std::cout << "TIME ADDFORCE " <<  (getTickCount() - timeaddforce)/getTickFrequency() << std::endl;

}

template <class DataTypes>
void FeatureMatchingForceField<DataTypes>::addForceMesh(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
    sofa::helper::vector< tri > triangles;
    triangles = sourceTriangles.getValue();

    bool reinitv = false;

    helper::vector< bool > sourcevisible = sourceVisible.getValue();
    helper::vector< int > indicesvisible = indicesVisible.getValue();
    helper::vector< bool > sourceborder = sourceBorder.getValue();


    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_1 (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_2 (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_temp(new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::PointXYZ>::Ptr filter_cloud_1 (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::PointXYZ>::Ptr filter_cloud_2 (new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::PointXYZ>::Ptr keypoints (new pcl::PointCloud<pcl::PointXYZ>);

            pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_1n (new pcl::PointCloud<pcl::PointXYZ>);
                    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_2n (new pcl::PointCloud<pcl::PointXYZ>);

    const VecCoord& xv = sourceVisiblePositions.getValue();
    const VecCoord&  tp = targetPositions.getValue();

    unsigned int nbs=xv.size(),nbt=tp.size();

        pcl::PointXYZ newPoint;
        for (unsigned int i=0; i<nbs; i++)
        {
        newPoint.z = xv[i][2];
        newPoint.x = xv[i][0];
        newPoint.y = xv[i][1];
        /*newPoint.r = 0;
        newPoint.g = 0;
        newPoint.b = 0;*/
        cloud_1->points.push_back(newPoint);

        }

        for (unsigned int i=0; i<nbt; i++)
        {
        newPoint.z = tp[i][2];
        newPoint.x = tp[i][0];
        newPoint.y = tp[i][1];
        /*newPoint.r = 0;
        newPoint.g = 0;
        newPoint.b = 0;*/
        cloud_2->points.push_back(newPoint);
        }


                std::cout << "\n remove NAN-Points" << std::endl;

                //remove NAN-Points
                std::vector<int> indices1, indices2;
                /*pcl::removeNaNFromPointCloud(*cloud_1, *cloud_1, indices1);
                pcl::removeNaNFromPointCloud(*cloud_2, *cloud_2, indices2);*/

                std::cout << "\n Filter PCD Files" << std::endl;

                //	Filter
                /*pcl::VoxelGrid<pcl::PointXYZRGB> vg;
                vg.setLeafSize(0.01f,0.01f,0.01f);
                vg.setInputCloud(cloud_1);
                vg.filter(*filter_cloud_1);
                vg.setInputCloud(cloud_2);
                vg.filter(*filter_cloud_2);*/

                std::cout << "\n Estimate Normals" << std::endl;

                // Normal-Estimation
                pcl::PointCloud<pcl::PointNormal>::Ptr norm_in(new pcl::PointCloud<pcl::PointNormal>);
                pcl::PointCloud<pcl::PointNormal>::Ptr norm_out(new pcl::PointCloud<pcl::PointNormal>);

                pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_in(new pcl::search::KdTree<pcl::PointXYZ>());
                pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_out(new pcl::search::KdTree<pcl::PointXYZ>());

                pcl::PointCloud<pcl::Normal>::Ptr norm_in1(new pcl::PointCloud<pcl::Normal>);
                pcl::PointCloud<pcl::Normal>::Ptr norm_out1(new pcl::PointCloud<pcl::Normal>);

                pcl::NormalEstimation<pcl::PointXYZ, pcl::PointNormal> ne;

                //Source-Cloud
                ne.setInputCloud(cloud_1);
                ne.setSearchSurface(cloud_1);
                ne.setSearchMethod(tree_in);
                ne.setRadiusSearch(0.04);
                ne.compute(*norm_in);


                //Target-Cloud
                ne.setInputCloud(cloud_2);
                ne.setSearchSurface(cloud_2);
                ne.setSearchMethod(tree_out);
                ne.setRadiusSearch(0.04);
                ne.compute(*norm_out);

            const float min_scale = 0.01;
            const int nr_octaves = 3;
            const int nr_scales_per_octave = 4;
            const float min_contrast = 0.001;
            const float radius = 0.05;

            // Copy the xyz info from cloud_xyz and add it to cloud_normals as the xyz field in PointNormals estimation is zero
            for(size_t i = 0; i<norm_in->points.size(); ++i)
            {
              norm_in->points[i].x = cloud_1->points[i].x;
              norm_in->points[i].y = cloud_1->points[i].y;
              norm_in->points[i].z = cloud_1->points[i].z;

              pcl::Normal norm;

              norm.normal_x = norm_in->points[i].normal_x;
              norm.normal_y = norm_in->points[i].normal_y;
              norm.normal_z = norm_in->points[i].normal_z;
              //if (!(norm.normal_x!=norm.normal_x || norm.normal_y!=norm.normal_y || norm.normal_z!=norm.normal_z))
              {
              norm_in1->points.push_back(norm);
              cloud_1n->points.push_back(cloud_1->points[i]);
              }
            }

            for(size_t i = 0; i<norm_out->points.size(); ++i)
            {
              norm_out->points[i].x = cloud_2->points[i].x;
              norm_out->points[i].y = cloud_2->points[i].y;
              norm_out->points[i].z = cloud_2->points[i].z;

              pcl::Normal norm;
              norm.normal_x = norm_out->points[i].normal_x;
              norm.normal_y = norm_out->points[i].normal_y;
              norm.normal_z = norm_out->points[i].normal_z;

              //std::cout << "norm " << norm.normal_x << " " << norm.normal_y << " " << norm.normal_z << std::endl;
              //if (!(norm.normal_x!=norm.normal_x || norm.normal_y!=norm.normal_y || norm.normal_z!=norm.normal_z))
              {
              norm_out1->points.push_back(norm);
              cloud_2n->points.push_back(cloud_2->points[i]);
              }
            }

            std::cout << "size " << cloud_2n->points.size() << " " << norm_out1->points.size() << std::endl;
            std::cout << "size " << cloud_1n->points.size() << " " << norm_in1->points.size() << std::endl;
            /*pcl::removeNaNFromPointCloud(*norm_in1, *norm_in1, indices1);
            pcl::removeNaNFromPointCloud(*norm_out1, *norm_out1, indices2);
            pcl::removeNaNFromPointCloud(*cloud_2, *cloud_2, indices1);
            pcl::removeNaNFromPointCloud(*cloud_1, *cloud_1, indices2);
            pcl::removeNaNFromPointCloud(*norm_in, *norm_in, indices1);
            pcl::removeNaNFromPointCloud(*norm_out, *norm_out, indices2);*/


            std::cout << "\n compute SIFT Keypoints" << std::endl;

        // Compute the SIFT keypoints
            pcl::SIFTKeypoint<pcl::PointNormal,pcl:: PointWithScale> sift_detector_1;
            pcl::SIFTKeypoint<pcl::PointNormal,pcl:: PointWithScale> sift_detector_2;
            pcl::search::KdTree<pcl::PointNormal>::Ptr tree (new pcl::search::KdTree<pcl::PointNormal>);
            pcl::PointCloud<pcl::PointWithScale> keypoints_temp_in;
            pcl::PointCloud<pcl::PointWithScale> keypoints_temp_out;

            sift_detector_1.setInputCloud(norm_in);
            sift_detector_2.setInputCloud(norm_out);
            sift_detector_1.setSearchMethod (tree);
            sift_detector_2.setSearchMethod (tree);
            sift_detector_1.setScales (min_scale, nr_octaves, nr_scales_per_octave);
            sift_detector_2.setScales (min_scale, nr_octaves, nr_scales_per_octave);
            sift_detector_1.setMinimumContrast (min_contrast);
            sift_detector_2.setMinimumContrast (min_contrast);

            /*sift_detector_1.setSearchSurface(norm_in);
            sift_detector_2.setSearchSurface(norm_out);*/
            sift_detector_1.setRadiusSearch (radius);
            sift_detector_2.setRadiusSearch (radius);

            sift_detector_1.compute(keypoints_temp_in);
            sift_detector_2.compute(keypoints_temp_out);

            std::cout << "No of SIFT points in the result are " << keypoints_temp_in.points.size () << std::endl;
            std::cout << "No of SIFT points in the result are " << keypoints_temp_out.points.size () << std::endl;


            std::vector<int> inliers_in;
            std::vector<int> inliers_out;


          pcl::PointCloud<pcl::PointXYZ>::Ptr keypoints_ptr_in(new pcl::PointCloud<pcl::PointXYZ>);
          pcl::PointCloud<pcl::PointXYZ>::Ptr keypoints_ptr_out(new pcl::PointCloud<pcl::PointXYZ>);
          pcl::PointCloud<pcl::PointXYZ>::Ptr keypoints_ptr_inliers_in(new   pcl::PointCloud<pcl::PointXYZ>);
          pcl::PointCloud<pcl::PointXYZ>::Ptr keypoints_ptr_inliers_out(new pcl::PointCloud<pcl::PointXYZ>);

          copyPointCloud (keypoints_temp_in , *keypoints_ptr_in);
          copyPointCloud (keypoints_temp_out , *keypoints_ptr_out);

          /*pcl::io::savePCDFileASCII("Keypoints_in.pcd", *keypoints_ptr_in);
          pcl::io::savePCDFileASCII("Keypoints_out.pcd",*keypoints_ptr_out);*/

          // Create the FPFH estimation class, and pass the input dataset+normals to it
          pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> pfh;
          pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> cloud_in_pfh;
          pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> cloud_out_pfh;
          pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);

          // Output datasets
          pcl::PointCloud<pcl::FPFHSignature33>::Ptr pfhs_in (new pcl::PointCloud<pcl::FPFHSignature33> ());
          pcl::PointCloud<pcl::FPFHSignature33>::Ptr pfhs_out (new pcl::PointCloud<pcl::FPFHSignature33> ());
          pcl::search::KdTree<pcl::PointXYZ>::Ptr tree_pfh_in(new pcl::search::KdTree<pcl::PointXYZ>());

          std::cout << "\n compute features" << std::endl;

          // Compute the features
          pfh.setInputCloud(keypoints_ptr_in);
          pfh.setSearchSurface(cloud_1n);
          pfh.setInputNormals(norm_in1);
          pfh.setRadiusSearch (0.05);
          pfh.setSearchMethod(tree_pfh_in);
          pfh.compute (*pfhs_in);

          pfh.setInputCloud (keypoints_ptr_out);
          pfh.setInputNormals(norm_out1);
          pfh.setSearchSurface(cloud_2n);
          pfh.setRadiusSearch (0.05);
          pfh.compute (*pfhs_out);

          pcl::registration::CorrespondenceEstimation<pcl::FPFHSignature33,pcl::FPFHSignature33> est;
          est.setInputCloud (pfhs_in);
          est.setInputTarget (pfhs_out);
          pcl::CorrespondencesPtr corr (new pcl::Correspondences);
          est.determineReciprocalCorrespondences (*corr);

          double inlierThreshold = 100;

          Eigen::Matrix4f transformation;
          boost::shared_ptr<pcl::Correspondences> corr_inliers(new pcl::Correspondences);
          pcl::registration::CorrespondenceRejectorSampleConsensus<pcl::PointXYZ> reg;
          reg.setInputSource(keypoints_ptr_in);
          reg.setInputTarget(keypoints_ptr_out);
          reg.setInlierThreshold(inlierThreshold);
          reg.setMaximumIterations(2000);
          reg.setInputCorrespondences(corr);
          reg.getCorrespondences(*corr_inliers);
          transformation = reg.getBestTransformation();

          std::cout << "\n compute Correspondences" << keypoints_ptr_in->size() << " " << keypoints_ptr_out->size() << std::endl;

          for (size_t i = 0; i < corr->size(); ++i){
                  std::cout << " \n Correspondences = " << corr_inliers->size() << (*corr)[i] << "\n" << std::endl;
                  std::cout << "\n compute Correspondences" << keypoints_ptr_in->points[(*corr)[i].index_query] << " " << keypoints_ptr_out->points[(*corr)[i].index_match] << std::endl;

          }



    /*detectKeypoints (source, source_keypoints_);
    detectKeypoints (target, target_keypoints_);

    extractDescriptors (source, source_keypoints_, source_features_);
    extractDescriptors (target, target_keypoints_, target_features_);

    findCorrespondences (sourcefeatures_, target_features_, source2target_);
    findCorrespondences (target_features_, source_features_, target2source_);

    filterCorrespondences ();*/


    /*    if (t%niterations.getValue() == 0)
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


    time = (double)getTickCount();

    double error = 0;
    int nerror = 0;

        if(tp.size()==0)
            for (unsigned int i=0; i<s.size(); i++) closestPos[i]=x[i];

        ind = 0;
        int ivis =0;
        sourcew.resize(s.size());
        int npointspen = 0;

            for (unsigned int i=0; i<s.size(); i++)
            {
                //serr<<"addForce() between "<<springs[i].m1<<" and "<<closestPos[springs[i].m1]<<sendl;
                if (t%(niterations.getValue()) == 0 && t > 1)
                {
                    this->addSpringForce(m_potentialEnergy,f,x,v, i, s[i]);
                    if (sourcevisible[i]) ivis++;
                }
            }

    _f.endEdit();*/

}

template <class DataTypes>
void FeatureMatchingForceField<DataTypes>::addSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
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
void FeatureMatchingForceField<DataTypes>::addStoredSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
{
    int a = spring.m1;

        f[a]=f_[a];
        Mat& m = this->dfdx[i];
                m = dfdx1[i];
}

template <class DataTypes>
double FeatureMatchingForceField<DataTypes>::computeError(Vector3 sourcePoint, Vector3 targetPoint)
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



template<class DataTypes>
void FeatureMatchingForceField<DataTypes>::addSpringDForce(VecDeriv& df,const  VecDeriv& dx, int i, const Spring& spring, double kFactor, double /*bFactor*/)
{
    const int a = spring.m1;
    const Coord d = -dx[a];
    Deriv dforce = this->dfdx[i]*d;
    dforce *= kFactor;
    df[a]+=dforce;
    //serr<<"addSpringDForce, a="<<a<<", b="<<b<<", dforce ="<<dforce<<sendl;
}

template <class DataTypes>
void FeatureMatchingForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams,DataVecDeriv& _df , const DataVecDeriv&  _dx )
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
void FeatureMatchingForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *m, SReal kFactor, unsigned int &offset)
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
void FeatureMatchingForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    const Vec<4,float> c(1,0,0,1);
    std::vector< Vector3 > points;
    //const vector<Spring>& springs = this->springs.getValue();
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    points.resize(0);

    std::cout << " XSIZE " << x.size() << std::endl;

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

}

}
}
} // namespace sofa

//#endif  /* SOFA_COMPONENT_FeatureMatchingForceField_INL */


