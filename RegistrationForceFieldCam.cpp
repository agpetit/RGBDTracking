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
    , cacheSize(initData(&cacheSize,(unsigned int)10,"cacheSize","number of closest points used in the cache to speed up closest point computation."))
	, cameraIntrinsicParameters(initData(&cameraIntrinsicParameters,Vector4(),"cameraIntrinsicParameters","camera parameters"))
	, blendingFactor(initData(&blendingFactor,(Real)1,"blendingFactor","blending between projection (=0) and attraction (=1) forces."))
    , outlierThreshold(initData(&outlierThreshold,(Real)7,"outlierThreshold","suppress outliers when distance > (meandistance + threshold*stddev)."))
    , normalThreshold(initData(&normalThreshold,(Real)0,"normalThreshold","suppress outliers when normal.closestPointNormal < threshold."))
    , projectToPlane(initData(&projectToPlane,false,"projectToPlane","project closest points in the plane defined by the normal."))
    , rejectBorders(initData(&rejectBorders,false,"rejectBorders","ignore border vertices."))
    , rejectOutsideBbox(initData(&rejectOutsideBbox,false,"rejectOutsideBbox","ignore source points outside bounding box of target points."))
    , springs(initData(&springs,"spring","index, stiffness, damping"))
	, sourceSurfacePositions(initData(&sourceSurfacePositions,"sourceSurface","Points of the surface of the source mesh."))
	, sourcePositions(initData(&sourcePositions,"sourcePositions","Points of the mesh."))
	, targetPositions(initData(&targetPositions,"targetPositions","Points of the surface of the source mesh."))
	, sourceTriangles(initData(&sourceTriangles,"sourceTriangles","Triangles of the source mesh."))
    , sourceNormals(initData(&sourceNormals,"sourceNormals","Normals of the source mesh."))
	, sourceSurfaceNormals(initData(&sourceSurfaceNormals,"sourceSurfaceNormals","Normals of the surface of the source mesh."))
	, showStrainsPerElement(initData(&showStrainsPerElement, false, "showStrainsPerElement", "  "))
	, plasticStrainsN(initData(&plasticStrainsN, "plasticStrains", "plastic strain per element"))
	, elasticStrainsN(initData(&elasticStrainsN, "elasticStrains", "elastic strain per element"))
	, totalStrainsN(initData(&totalStrainsN, "totalStrains", "total strain per element"))
	, plasticStrainsPerNode(initData(&plasticStrainsPerNode, "plasticStrainsPerNode", "plastic strain per node"))
	, elasticStrainsPerNode(initData(&elasticStrainsPerNode, "elasticStrainsPerNode", "elastic strain per node"))
	, totalStrainsPerNode(initData(&totalStrainsPerNode, "totalStrainsPerNode", "total strain per node"))
	, vonMisesStress(initData(&vonMisesStress, "vonMisesStress", "vonmisesstress per element"))
    , barycenter(initData(&barycenter,"barycenter","Barycenter of the mesh."))
	, drawSource(initData(&drawSource,false,"drawSource"," "))
	, drawTarget(initData(&drawTarget,false,"drawTarget"," "))
	, drawContour(initData(&drawContour,false,"drawContour"," "))
	, showArrowSize(initData(&showArrowSize,0.01f,"showArrowSize","size of the axis."))
	, drawMode(initData(&drawMode,0,"drawMode","The way springs will be drawn:\n- 0: Line\n- 1:Cylinder\n- 2: Arrow."))
    , drawColorMap(initData(&drawColorMap,true,"drawColorMap","Hue mapping of distances to closest point"))
    , theCloserTheStiffer(initData(&theCloserTheStiffer,false,"theCloserTheStiffer","Modify stiffness according to distance"))
	, useContour(initData(&useContour,false,"useContour","Emphasize forces close to the target contours"))
	, useVisible(initData(&useVisible,true,"useVisible","Use the vertices of the viisible surface of the source mesh"))
	, visibilityThreshold(initData(&visibilityThreshold,(Real)0.001,"visibilityThreshold","Threshold to determine visible vertices"))
	, useIntensity(initData(&useIntensity, false,"useIntensity","Use intensity features"))
	, alphaIntensity(initData(&alphaIntensity,(Real)0.00004,"alphaIntensity","Weight of intensity features"))
	, useKLTPoints(initData(&useKLTPoints, false,"useKLTPoints","Use KLT Points"))
	, useCCD(initData(&useCCD, false,"useCCD","Use CCD"))
	, useRealData(initData(&useRealData,true,"useRealData","Use real data"))
	, useGroundTruth(initData(&useGroundTruth,false,"useGroundTruth","Use the vertices of the visible surface of the source mesh"))
	, useSensor(initData(&useSensor,false,"useSensor","Use the sensor"))
	//, useKalman(initData(&useKalman,false,"useKalman","Use the Kalman filter"))
	, sensorType(initData(&sensorType, 0,"sensorType","Type of the sensor"))
	, generateSynthData(initData(&generateSynthData, false,"generateSynthData","Generate synthetic data"))
	, niterations(initData(&niterations,3,"niterations","Number of iterations in the tracking process"))
	, nimages(initData(&nimages,1500,"nimages","Number of images to read"))
	, windowKLT(initData(&windowKLT,5,"windowKLT","window for the KLT tracker"))
	, useDistContourNormal(initData(&useDistContourNormal,false,"useDistContourNormal","use normals to contours"))
	, ipad(initData(&ipad,"ipad"," ip address",false))
	, dataPath(initData(&dataPath,"dataPath","Path for data writings",false))
	, useMassSpring(initData(&useMassSpring,false,"useMassSpring","Use mass spring model"))
	, sigmaWeight(initData(&sigmaWeight,(Real)10,"sigmaWeight","Number of iterations in the tracking process"))
	,translation(initData(&translation,"translation", "translation parameters"))
	,rotation(initData(&rotation,"rotation", "rotation parameters"))
	,errorfunction(initData(&errorfunction,"errorfunction", "error"))
        , viewportWidth(initData(&viewportWidth,640,"viewportWidth","Width of the viewport"))
        , viewportHeight(initData(&viewportHeight,480,"viewportHeight","Height of the viewport"))
{
	nimages = 1500;
	pcl = false;
	disp = false;
	iter_im = 0;
	timeTotal = 0;

	timef = 0;//(double)getTickCount();
	timei = 0;
	timeOverall = 0;
	timeResolution = 0;
	timer = 0;
	timeSecondPass = 0;
	timeFirstPass = 0;
	
	rectRtt.x = 0;
	rectRtt.y = 0;
	rectRtt.height = 480;
	rectRtt.width = 640;
	timeOverall1 = 0;
    // Tracker parameters
    tracker.setTrackerId(1);
    //tracker.setOnMeasureFeature(&modifyFeature);
    tracker.setMaxFeatures(200);
    tracker.setWindowSize(10);
    tracker.setQuality(0.01);
    tracker.setMinDistance(10);
    tracker.setHarrisFreeParameter(0.04);
    tracker.setBlockSize(9);
    tracker.setUseHarris(1);
    tracker.setPyramidLevels(3); 
	
	tracker1.setTrackerId(1);
    tracker1.setMaxFeatures(200);
    tracker1.setWindowSize(10);
    tracker1.setQuality(0.01);
    tracker1.setMinDistance(10);
    tracker1.setHarrisFreeParameter(0.04);
    tracker1.setBlockSize(9);
    tracker1.setUseHarris(1);
    tracker1.setPyramidLevels(3);

}

template <class DataTypes>
RegistrationForceFieldCam<DataTypes>::~RegistrationForceFieldCam()
{
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::initCamera()
{
	
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::reinit()
{

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue(); 	//RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
	this->clearSprings(x.size());
	
    for(unsigned int i=0;i<x.size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());	

}

bool g_bQATest = false;
int  g_nDevice = 0;

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::init()
{

	std::cout << " ipad " << ipad.getValue().c_str() << std::endl;
	
	std::string ip_recipient;	
	ip_recipient = ipad.getValue();

// Socket creation *****************************************************************

	err = p_helper_socket::socket_create_UDP_Sender(sock_sender, 6999, /*ip_recipient*/"192.168.3.29");

	if(err < 0)
	cout<<"errore sender"<<endl;
	
	//**********************************************************************************

    this->Inherit::init();
    core::objectmodel::BaseContext* context = this->getContext();

    if(!(this->mstate)) this->mstate = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(context->getMechanicalState());

	    // Get source triangles
    /*if(!sourceTriangles.getValue().size()) {
        sofa::component::loader::MeshObjLoader *meshobjLoader;
        this->getContext()->get( meshobjLoader, core::objectmodel::BaseContext::Local);
        if (meshobjLoader) {sourceTriangles.virtualSetLink(meshobjLoader->triangles); sout<<"imported triangles from "<<meshobjLoader->getName()<<sendl;
		}
    }*/
    // Get source normals
    if(!sourceNormals.getValue().size()) serr<<"normals of the source model not found"<<sendl;

    // add a spring for every input point
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue(); 			//RDataRefVecCoord x(*this->getMState()->read(core::ConstVecCoordId::position()));
	this->clearSprings(x.size());
	
	npoints = x.size();
	
    for(unsigned int i=0;i<x.size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());	
	
			//initCamera();
	bool opt_device = false;
    bool opt_display = true;
    bool use_cuda = true;
    bool opt_click_allowed = true;
    int start_image = 0;
	
	//if (!useRealData.getValue()) useGroundTruth.setValue(true);
        cv::Rect ROI(160, 120, 320, 240);
	
	Vector4 camParam = cameraIntrinsicParameters.getValue();
	
	rgbIntrinsicMatrix(0,0) = camParam[0];
	rgbIntrinsicMatrix(1,1) = camParam[1];
	rgbIntrinsicMatrix(0,2) = camParam[2];
	rgbIntrinsicMatrix(1,2) = camParam[3];
	
	std::cout << " camparam  " << camParam[0]<< " " << camParam[1] << " " << camParam[2]<< " " << camParam[3] << std::endl;
	
	/*rgbIntrinsicMatrix(0,0) = 275.34;
	rgbIntrinsicMatrix(1,1) = 275.34;
	//rgbIntrinsicMatrix(0,2) = 157.25;
	//rgbIntrinsicMatrix(1,2) = 117.75;
	rgbIntrinsicMatrix(0,2) = 160;
	rgbIntrinsicMatrix(1,2) = 120;*/

        glEnable(GL_BLEND);

     glEnable(GL_DEPTH_TEST);

     // Request Stencil Buffer support

	sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());

	root->get(rgbddataprocessing);
	root->get(meshprocessing);	
	root->get(closestpoint);
	root->get(rendertexturear);
	root->get(dataio);	
		
	if (showStrainsPerElement.getValue())
		npasses = 2;
	else npasses = 1;
	
		
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::resetSprings()
{
	
this->clearSprings(sourceVisiblePositions.getValue().size());	
for(unsigned int i=0;i<sourceVisiblePositions.getValue().size();i++) this->addSpring(i, (Real) ks.getValue(),(Real) kd.getValue());	
	
}

template<class DataTypes>
void RegistrationForceFieldCam<DataTypes>::computeTargetNormals()
{
	cv::Mat normals, normals0;
	cv::Sobel(depth,normals,-1,1,0,3);
	cv::namedWindow("normals");
	cv::imshow("normals", normals);
	
}

template<class DataTypes>
void RegistrationForceFieldCam<DataTypes>::setViewPoint()
{
	Eigen::Affine3f scene_sensor_pose (Eigen::Affine3f::Identity ());
	pcl::PointCloud<pcl::PointXYZRGB>& point_cloud = *rgbddataprocessing->target;
		
	scene_sensor_pose = Eigen::Affine3f (Eigen::Translation3f (point_cloud.sensor_origin_[0],
                                                             point_cloud.sensor_origin_[1],
                                                             point_cloud.sensor_origin_[2])) * 
															 Eigen::Affine3f (point_cloud.sensor_orientation_);
															 
  Eigen::Affine3f viewer_pose = scene_sensor_pose;
  Eigen::Vector3f pos_vector = viewer_pose * Eigen::Vector3f(0, 0, 0);  
  Eigen::Vector3f look_at_vector = viewer_pose.rotation () * Eigen::Vector3f(0, 0, 1) + pos_vector;
  Eigen::Vector3f up_vector = viewer_pose.rotation () * Eigen::Vector3f(0, -1, 0);
  /*viewer.setCameraPosition (pos_vector[0], pos_vector[1], pos_vector[2],
                            look_at_vector[0], look_at_vector[1], look_at_vector[2],
                            up_vector[0], up_vector[1], up_vector[2]);*/
							
        int hght = viewportHeight.getValue();
        int wdth = viewportWidth.getValue();
    sofa::gui::GUIManager::SetDimension(wdth,hght);
	sofa::gui::BaseGUI *gui = sofa::gui::GUIManager::getGUI();
        sofa::gui::BaseViewer * viewer = gui->getViewer();

			{	
			Vec3d position;
            Quat orientation;
			
			orientation[0] = point_cloud.sensor_orientation_.w ();
			orientation[1] = point_cloud.sensor_orientation_.x ();
			orientation[2] = point_cloud.sensor_orientation_.y ();
			orientation[3] = point_cloud.sensor_orientation_.z ();
			
			position[0] = point_cloud.sensor_origin_[0];
			position[1] = point_cloud.sensor_origin_[1];
			position[2] = point_cloud.sensor_origin_[2];

    //if (currentCamera)

                        //sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
                        //sofa::component::visualmodel::BaseCamera::SPtr currentCamera;// = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();
                        //root->get(currentCamera);

                        //double fov = (double)currentCamera->p_fieldOfView;

                        //std:cout << " FoV " << currentCamera->getFieldOfView() << std::endl;

                        //currentCamera->setViewport(480,640);
                        //currentCamera->p_fieldOfView.setValue(10);

                        //currentCamera->setView(position, orientation);


                        viewer->setView(position, orientation);

                        /*sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
                        sofa::component::visualmodel::BaseCamera::SPtr currentCamera;// = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();
                        root->get(currentCamera);

                        double znear = currentCamera->getZNear();
                        double zfar = currentCamera->getZFar();

                        std::cout << " znear0 " << znear << " zfar0 " << zfar << std::endl;

                        gui->redraw();*/

                    //glutPostRedisplay();
                                //glPushAttrib( GL_LIGHTING_BIT | GL_ENABLE_BIT | GL_LINE_BIT | GL_CURRENT_BIT);
        //glPopAttrib();
                        }
	
}

void getEllipses(std::vector<std::vector<cv::Point> >& contours, std::vector<cv::RotatedRect>& ellipses) {
    ellipses.clear();
    cv::Mat img0(cv::Size(320,240), CV_8UC3);
    for (unsigned i = 0; i<contours.size(); i++) {
        if (contours[i].size() >= 5) {
            cv::RotatedRect temp = cv::fitEllipse(cv::Mat(contours[i]));
			{
                //cout << "Reject ellipse " << i << endl;
                cv::drawContours(img0, contours, i, cv::Scalar(0,255,0), -1, 8);
                cv::ellipse(img0, temp, cv::Scalar(255,255,0), 2, 8);
                cv::imshow("Ellipses", img0);
                cv::waitKey(1);
            }
        }
    }
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::addForce(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
		int t = (int)this->getContext()->getTime();	
//if( (t<41) || ( t%20 == 0 || (t+1)%20==0  || (t+2)%20==0 || (t+3)%20==0 || (t+4)%20==0 ) )
                         // for (int i = 0; i <5; i++)
			addForceMesh(mparams, _f, _x, _v);
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::mapKLTPointsTriangles ( sofa::helper::vector< tri > &triangles)
{
    int outside = 0;
        const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    //const sofa::core::topology::BaseMeshTopology::SeqTriangles& triangles = this->fromTopology->getTriangles();
    sofa::helper::vector<Mat3x3d> bases;
    sofa::helper::vector<Vector3> centers;
	
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
			if (sourceVisible[triangles[t][2]] && sourceVisible[triangles[t][1]] && sourceVisible[triangles[t][0]])
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

    if (tt < 0.0000001 || tt!=tt || vv!=vv || uu!=uu)
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
int id;

	 for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
                Mat3x3d m,mt;
				double xt, yt;
			     Vector3 xim0,xim1,xim2;
				 const int kk = k;
//	tracker.getFeature(kk, id, xp, yp);
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
			if (sourceVisible[triangles[t][2]] && sourceVisible[triangles[t][1]] && sourceVisible[triangles[t][0]])
				 {
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
                    double d = max ( max ( -v[0],-v[1] ),max ( ( v[2]<0?-v[2]:v[2] )-0.01,v[0]+v[1]-1 ) );
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
					std::cout << " mapping " << kk << " id " << id << " index " << index << " coefs " << coefs[0] << " " << coefs[1] << " " << coefs[2] << std::endl;
					mappingkltcoef[id] = coefs;
					mappingkltind[id] = index;
			int x_u_0 = (int)(x[triangles[index][0]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[index][0]][2] + rgbIntrinsicMatrix(0,2));
			int x_v_0 = (int)(x[triangles[index][0]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[index][0]][2] + rgbIntrinsicMatrix(1,2));
				}
			
            }
        }
	}
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::KLTPointsTo3D()
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
	
	 znear = 0.0716081;
	 zfar  = 72.8184;
	//std::cout << " znear " << znear << " zfar " << zfar << std::endl;

            Vector3 pos;
            Vector3 col;
float xp, yp;
int id;

//cv::imwrite("depthpp.png", depth);

	 for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
                Mat3x3d m,mt;
				double xt, yt;
			     Vector3 xim0,xim1,xim2;
				 const int kk = k;
	//tracker.getFeature(kk, id, xp, yp);
				int n0 = (int)yp;
				 int m0 = (int)xp;
				 float depthValue;
							if (!useRealData.getValue())
				 			depthValue = (float)depth.at<float>(2*yp,2*xp);
							else depthValue = (float)depth.at<float>(yp,xp);

			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)color.at<Vec4b>(yp,xp)[3];
			if ( depthValue>0 && depthValue < 1)                // if depthValue is not NaN
			{				
				double clip_z = (depthValue - 0.5) * 2.0;
			//double clip_z = (depths1[j-rectRtt.x+(i-rectRtt.y)*(rectRtt.width)] - 0.5) * 2.0;
                if (!useRealData.getValue()) pos[2] = -2*znear*zfar/(clip_z*(zfar-znear)-(zfar+znear));
				else pos[2] = depthValue;
				pos[0] = (xp - rgbIntrinsicMatrix(0,2)) * pos[2] * rgbFocalInvertedX;
				pos[1] = (yp - rgbIntrinsicMatrix(1,2)) * pos[2] * rgbFocalInvertedY;
				targetpos[id]=pos;
				//std::cout << " id " <<pos[2] << std::endl;
				//std::cout << " id " << id << " size targetpos " << targetpos.size() << " pos " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
				
			}
		
            }
					//getchar();	
    const VecCoord&  p = targetpos;
	targetKLTPositions.setValue(p);

}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::CCDPointsTo3D()
{
	
    int outside = 0;
	
	float rgbFocalInvertedX = 1/rgbIntrinsicMatrix(0,0);	// 1/fx
	float rgbFocalInvertedY = 1/rgbIntrinsicMatrix(1,1);	// 1/fy

	VecCoord targetpos;
	targetpos.resize(pointsCCD.size());
	
	sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
	sofa::component::visualmodel::BaseCamera::SPtr currentCamera;// = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();
	root->get(currentCamera);
	
	double znear = currentCamera->getZNear();
	double zfar = currentCamera->getZFar();
	
	 /*znear = 0.0716081;
	 zfar  = 72.8184;*/
	 
            Vector3 pos;
            Vector3 col;
			float xp, yp;
			int id;

	 for (unsigned int k = 0; k < pointsCCD.size(); k++){
                Mat3x3d m,mt;
				double xt, yt;
			     Vector3 xim0,xim1,xim2;
				 const int kk = k;
				 				//std::cout << " pos " << std::endl;
	            xp = ccd.pointsccdmin[k].xu;
				yp = ccd.pointsccdmin[k].xv;
				
				int n0 = (int)yp;
				 int m0 = (int)xp;
				 
				//std::cout << " xp " <<xp << " yp " << yp << std::endl;
				 float depthValue;
							if (!useRealData.getValue())
				 			depthValue = (float)depth.at<float>(2*yp,2*xp);
							else depthValue = (float)depth.at<float>(yp,xp);

			//depthValue =  1.0 / (depthValue*-3.0711016 + 3.3309495161);;
			int avalue = (int)color.at<Vec4b>(yp,xp)[3];
			if ( depthValue>0 && depthValue < 1)                // if depthValue is not NaN
			{	
				//std::cout << " pos " << std::endl;
			
				double clip_z = (depthValue - 0.5) * 2.0;
                if (!useRealData.getValue()) pos[2] = -2*znear*zfar/(clip_z*(zfar-znear)-(zfar+znear));
				else pos[2] = depthValue;
				
				/*pos[2] = pointsCCD[k].Z;
				pos[1] = pointsCCD[k].Y;
				pos[0] = pointsCCD[k].X;*/
				//pos[2] = pointsCCD[k].Z;
				
				pos[0] = (xp - rgbIntrinsicMatrix(0,2)) * pos[2] * rgbFocalInvertedX;
				pos[1] = (yp - rgbIntrinsicMatrix(1,2)) * pos[2] * rgbFocalInvertedY;
				targetpos[k]=pos;
				
				//std::cout << " pos " << rgbIntrinsicMatrix(0,2) << " " << rgbIntrinsicMatrix(1,2) << " "  << rgbIntrinsicMatrix(0,0) << std::endl;
				
			}
					
            }
    const VecCoord&  p = targetpos;
	targetCCDPositions.setValue(p);
}


/*template <class DataTypes, class DepthTypes>
void RegistrationForceFieldCam<DataTypes, DepthTypes>::addPointInTriangle ( const int triangleIndex, const Real* baryCoords )
{
    map2d.resize ( map2d.size() +1 );
    MappingData2D& data = *map2d.rbegin();
    data.in_index = triangleIndex;
    data.baryCoords[0] = ( Real ) baryCoords[0];
    data.baryCoords[1] = ( Real ) baryCoords[1];
    return map2d.size()-1;
}*/

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::initCCD()
{
	ReadAccessor< Data< VecCoord > > xcp0(sourceContourPositions);	
	pointsCCD.resize(xcp0.size());
	pointCCD pCCD;
	std::cout << " pointsCCD size " << pointsCCD.size() << std::endl;
		
		for (unsigned int i = 0; i  <xcp0.size() ; i++)
			{
			pCCD.X = xcp0[i][0];
			pCCD.Y = xcp0[i][1];
			pCCD.Z = xcp0[i][2];
			pCCD.x = xcp0[i][0]/xcp0[i][2];
			pCCD.y = xcp0[i][1]/xcp0[i][2];
			int x_u = (int)(xcp0[i][0]*rgbIntrinsicMatrix(0,0)/xcp0[i][2] + rgbIntrinsicMatrix(0,2));
			int x_v = (int)(xcp0[i][1]*rgbIntrinsicMatrix(1,1)/xcp0[i][2] + rgbIntrinsicMatrix(1,2));
			pCCD.xu = x_u;
			pCCD.xv = x_v;
			//double norm = sqrt(Sx.at<float>(x_u,x_v)*Sx.at<float>(x_u,x_v) + Sy.at<float>(x_u,x_v)*Sy.at<float>(x_u,x_v));
			{/*pCCD.nx = cos((double)ori.at<float>(x_v,x_u) + 3.1416/2.0);
			pCCD.ny = sin((double)ori.at<float>(x_v,x_u) + 3.1416/2.0);*/
						pCCD.nx = normalsContour[i].x;
						pCCD.ny = normalsContour[i].y;
						
			}
			pointsCCD[i] = pCCD;
			//std::cout << " pccd " << x_u << "  " << x_v << " " << pCCD.nx << " " << pCCD.ny << " "/* << (double)gray.at<uchar>(x_v,x_u) << " " << (double)Sx.at<float>(x_v,x_u) << " ori " << (double)ori.at<float>(x_v,x_u)*/ << std::endl;
			}
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::updateCCD()
{			
int kcp = 0;
const VecCoord x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
	for (unsigned int i=0; i<x.size(); i++)
		    {
			if(meshprocessing->sourceBorder[i]){
			pointsCCD[kcp].X = x[i][0];
			pointsCCD[kcp].Y = x[i][1];
			pointsCCD[kcp].Z = x[i][2];
			pointsCCD[kcp].x = x[i][0]/x[i][2];
			pointsCCD[kcp].y = x[i][1]/x[i][2];
			int x_u = (int)(x[i][0]*rgbIntrinsicMatrix(0,0)/x[i][2] + rgbIntrinsicMatrix(0,2));
			int x_v = (int)(x[i][1]*rgbIntrinsicMatrix(1,1)/x[i][2] + rgbIntrinsicMatrix(1,2));
			pointsCCD[kcp].xu = x_u;
			pointsCCD[kcp].xv = x_v;
			//std::cout << " xi " << x[i][0] << " "  << x[i][1] << " " << x[i][2] << std::endl;
			kcp++;
			}
			}
			
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::addForceMesh(const core::MechanicalParams* mparams,DataVecDeriv& _f , const DataVecCoord& _x , const DataVecDeriv& _v )
{
		
	sofa::helper::vector< tri > triangles;
	triangles = sourceTriangles.getValue();

	int t = (int)this->getContext()->getTime();
	
	double t00 = this->getContext()->getTime();
		
	//double timef = 0;
	timeii = timef;
	timef = (double)getTickCount();
	//if (t > 0 && t%niterations.getValue() == niterations.getValue()-1)
	timeOverall += (timef - timei)/getTickFrequency();
	timeResolution += (timef - timei)/getTickFrequency();
	
	timeSecondPass = (timef - timeii)/getTickFrequency() ;
		
	if ((t+1)%niterations.getValue() == 0)
	{timeFirstPass = timeResolution;
	timeOverall1 = timeOverall;
	}
	
	{
    cout <<" t " << t00 << " Time resolution " << timeResolution << " time overall " << timeOverall << " time total " << timeSecondPass << endl;
	}
	
	double timeT = (double)getTickCount();
		
	if (t > 0 && t%niterations.getValue() == 0)//niterations.getValue()-1)
	{
	std::cout << " t " << t << " " << t/niterations.getValue() <<  " time overall " << timeOverall << std::endl;
	timeFile << t/niterations.getValue();	
	timeFile << "\t";
	timeFile << timeOverall;
	timeFile << "\t";
	timeFile << timeSeg;
	timeFile << "\t";
	timeFile << timeRigid;
	timeFile << "\t";
	timeFile << timeSourceContour;
	timeFile << "\t";
	timeFile << timeAddforce;
	timeFile << "\t";
	timeFile << timeResolution;
	timeFile << "\t";
	timeFile << timeFirstPass;
	timeFile << "\t";
	timeFile << timeSecondPass;
	timeFile << "\t";
	timeFile << timeOverall - timeSecondPass;
	timeFile << "\t";
	timeFile << timeOverall1;
	timeFile << "\n";
	}

	timeFile.close();

	
	cout << "Time t " <<  ((double)getTickCount() - timeT)/getTickFrequency() << endl;	
			
		//cv::Mat color11 = imconv->color;
		//cv::imwrite("color.png", color11)
	double timeAcq0 = (double)getTickCount();
	if (useRealData.getValue())
	{
	if (useSensor.getValue()){
		sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
		typename sofa::core::objectmodel::ImageConverter<DataTypes,DepthTypes>::SPtr imconv;// = root->getNodeObject<sofa::component::visualmodel::InteractiveCamera>();
		root->get(imconv);
		color = imconv->color;
                color_1 = imconv->color_1;
		depth = imconv->depth;
		depth00 = depth.clone();
		//cv::imwrite("depth00.png", depth00);
	}
	 else {
		color = dataio->color;
		depth = dataio->depth;
		depth00 = depth.clone();
		color_1 = dataio->color_1;
		color_5 = dataio->color_5.clone();
		color_4 = dataio->color_4.clone();
		color_3 = dataio->color_3.clone();
		color_2 = dataio->color_2.clone();
	 }
	}
	else dataio->readData();
	
	double timeAcq1 = (double)getTickCount();
   cout <<"time acq " << (timeAcq1 - timeAcq0)/getTickFrequency()<< " t " << t << endl;
	
	bool reinitv = false;
	if (t == 0)
	{
        timeFile.open ("timesFract.txt");
		timeFile.clear();
		errorFile.open ("errorsGeom.txt");
		errorFile.clear();
		if(useGroundTruth.getValue())
		{
		std::string datafile = dataPath.getValue() + "Gt.txt";
		errorGtFile.open(datafile.c_str());
		errorGtFile.clear();
		if(useMassSpring.getValue())
		correspMassSpringFile.open ("correspMassSpring.txt");
		//correspMassSpringFile.clear();
		}
	if (useRealData.getValue())
		targetPositions.setValue(rgbddataprocessing->getTargetPositions());
	
		timeOverall = 0;
		timeTotal = 0;
		timeAddforce = 0;
		timeSourceContour = 0;
		timeResolution = 0;
		
		color_init = color;
		
		if (useKLTPoints.getValue())
		{
		
	   cvtColor(color,gray,CV_BGR2GRAY);	
	   cv::Mat gray1;
       cv::flip(gray,gray1,0);
	   vpImageConvert::convert(gray,vpI);
		
	  display.init(vpI, 100, 100,"Display...") ;
      // Display the image
      vpDisplay::display(vpI) ;
      vpDisplay::flush(vpI) ;
	  
    // Point detection using Harris. In input we have an OpenCV image
    tracker.initTracking(gray);
	tracker1.initTracking(gray);	
	tracker1.display(vpI, vpColor::red);
		
	    }
	}
	else
	{
		 if (t > 0 && t%niterations.getValue() == 0){
		
		timeOverall = 0;
		timeAddforce = 0;
		timeSourceContour = 0;
		timeResolution = 0;
		
		if (t%(niterations.getValue()) == 0)
		{color_init = color;
		//cv::imwrite("color_init.png", color_init);
		}
		
		//std::cout << " ssize 0 " << (this->springs.getValue()).size() << std::endl;
		
			double time0 = (double)getTickCount();
			timeT = (double)getTickCount();

                if (npoints != (this->mstate->read(core::ConstVecCoordId::position())->getValue()).size())
		{
			reinit();
			reinitv = true;
		}
					
                npoints = (this->mstate->read(core::ConstVecCoordId::position())->getValue()).size();
		
		if(useRealData.getValue())
		{	
		if(!useContour.getValue()){
		targetPositions.setValue(rgbddataprocessing->getTargetPositions());
		}
		else {
		targetPositions.setValue(rgbddataprocessing->getTargetPositions());
		targetContourPositions.setValue(rgbddataprocessing->getTargetContourPositions());
		//targetBorder = rgbddataprocessing->targetBorder;
		targetWeights = rgbddataprocessing->targetWeights;
		}
		
                }

		//cv::Mat oriMap = orientationMap(mag, ori, 1.0);
		time0 = ((double)getTickCount() - time0)/getTickFrequency();
		cout << "Time segment + extract PCD " << time0 << endl;
		timeSeg = time0;
		timeOverall += time0;
		
		double time1 = (double)getTickCount();

		if (useVisible.getValue()) 			
		{
				time1 = (double)getTickCount();
			//if (t%(npasses + niterations.getValue() - 1) ==0 )
				{
				sourceVisible = meshprocessing->sourceVisible;
				sourceVisiblePositions.setValue(meshprocessing->getSourceVisiblePositions());
				indicesVisible = meshprocessing->indicesVisible;
				depthMap = meshprocessing->depthMap.clone();
				std::cout << " " << indicesVisible.size() << std::endl;

			time1 = ((double)getTickCount() - time1)/getTickFrequency();
			cout << "Time get source visible " << time1 << endl;
			
		if(useContour.getValue()){
		   sourceWeights = meshprocessing->sourceWeights;
		   sourceContourPositions.setValue(meshprocessing->getSourceContourPositions());
           cout << "size contour positions " << meshprocessing->getSourceContourPositions().size() << endl;
			}
				}
		   
		         double time = vpTime::measureTimeMs();
			  	   cv::Mat gray1;
				   cv::flip(gray,gray1,0);
				   
		   		if (useKLTPoints.getValue() && t >= 3){
				if (t == 3 || t%(windowKLT.getValue()*niterations.getValue()) ==0 )
				{
			tracker.initTracking(gray);
			tracker1.initTracking(gray);
			mappingkltcoef.resize(tracker.getMaxFeatures());
			mappingkltind.resize(tracker.getMaxFeatures());
		   	mapKLTPointsTriangles(triangles);
				}
			
			  cvtColor(color,gray,CV_BGR2GRAY);	
			  
			cv::flip(gray,gray1,0);
			vpImageConvert::convert(gray,vpI);
		      //vpImageConvert::convert(gray,vpI)
        // Display the image
        vpDisplay::display(vpI) ;
      // Tracking of the detected points
        tracker.track(gray);
        tracker1.track(gray);

        // Display the tracked points
        tracker1.display(vpI, vpColor::red);
        vpDisplay::flush(vpI) ;
		vpImage<vpRGBa> Icol;
		cv::Mat icol;
		display.getImage(Icol);
		vpImageConvert::convert(Icol,icol);
		imgklt = new cv::Mat;
       //*imgl = color;
	    *imgklt = icol;		

		if (t%npasses == 0)
        dataio->listimgklt.push_back(imgklt);		  

		KLTPointsTo3D();
				}
				
				}
				
			timeRigid = time1;
			timeOverall += time1;
			
		cv::Mat gray0,gray;
		//cvtColor( rtd, gray0, CV_BGR2GRAY );
		//cv::GaussianBlur( rtd, gray, Size( 3, 3), 0, 0 );
                gray = meshprocessing->depthMap;
		//cv::imwrite("gray.png", gray);
		
		cv::Mat Sx;
		//cv::Sobel(gray, Sx, CV_32F, 1, 0, 7);

		cv::Mat Sy;
		//cv::Sobel(gray, Sy, CV_32F, 0, 1, 7);

		cv::Mat mag, ori;
		
		ReadAccessor< Data< VecCoord > > xcp0(sourceContourPositions);	
		normalsContour.resize(xcp0.size());
		
		cv::Mat gray1;
                //cv::GaussianBlur( gray, gray1, Size( 9, 9), 0, 0 );
                //cv::imwrite("gray1.png",gray1);

		if (useContour.getValue()){
		cout << "Rigid 0 " << xcp0.size() << endl;
	    Eigen::Matrix<float,3,2> gradient;
		for (unsigned int i = 0; i  < xcp0.size() ; i++)
			{
			int x_u = (int)(xcp0[i][0]*2*rgbIntrinsicMatrix(0,0)/xcp0[i][2] + 2*rgbIntrinsicMatrix(0,2));
			int x_v = (int)(xcp0[i][1]*2*rgbIntrinsicMatrix(1,1)/xcp0[i][2] + 2*rgbIntrinsicMatrix(1,2));

			//double norm = sqrt(Sx.at<float>(x_u,x_v)*Sx.at<float>(x_u,x_v) + Sy.at<float>(x_u,x_v)*Sy.at<float>(x_u,x_v));
			{/*pCCD.nx = cos((double)ori.at<float>(x_v,x_u) + 3.1416/2.0);
			pCCD.ny = sin((double)ori.at<float>(x_v,x_u) + 3.1416/2.0);*/
			
			gradient(0,0) = (2047.0 *(gray1.at<uchar>(x_v,x_u+1) - gray1.at<uchar>(x_v,x_u-1)) + 913.0*(gray1.at<uchar>(x_v,x_u+2) - gray1.at<uchar>(x_v,x_u-2))+112.0 *(gray1.at<uchar>(x_v,x_u+3) - gray1.at<uchar>(x_v,x_u-3)))/8418.0;
            gradient(0,1) = (2047.0 *(gray1.at<uchar>(x_v+1,x_u) - gray1.at<uchar>(x_v-1,x_u)) + 913.0*(gray1.at<uchar>(x_v+2,x_u) - gray1.at<uchar>(x_v-2,x_u))+112.0 *(gray1.at<uchar>(x_v+3,x_u) - gray1.at<uchar>(x_v-3,x_u)))/8418.0;
							cv::Point2f normal;
						if(gradient(0,0) == 0) {
						normal.x = cos(CV_PI/2.0);
												
						normal.y = sin((CV_PI/2.0));
						}
						else{normal.x = cos((atan(gradient(0,1) / gradient(0,0))) );
									
						normal.y = sin((atan(gradient(0,1) / gradient(0,0))));}
						normalsContour[i] = normal;
						
			}
			}
			
		if (useCCD.getValue())
		{
			initCCD();
			ccd.init(rgbIntrinsicMatrix, pointsCCD);
		
		}
		}			
	}
		std::cout << " source size 1 " << sourceVisiblePositions.getValue().size() << std::endl;
		if (useVisible.getValue() && t >= 3 && t%niterations.getValue()!= 0) 
		{
		sourceVisiblePositions.setValue(meshprocessing->getSourceVisiblePositions());
		std::cout << " source size 2 " << sourceVisiblePositions.getValue().size() << std::endl;

		}
}

	//closestpoint->initSourceSurface();
		
        if (t < 2)
        setViewPoint();
	
		sofa::gui::BaseGUI *gui = sofa::gui::GUIManager::getGUI();
		sofa::gui::BaseViewer * viewer = gui->getViewer();
		std::string opath00 = "out/images38y/img2%06d.png";
		int iterm;
	if (niterations.getValue() == 2)
		iterm = 0;
		else iterm = 0;

        if (t >= 3 ){
    if ( t%npasses == iterm)
	{

        std::cout << " source size 3 " << std::endl;
	double timertt = (double)getTickCount();
	rtt = new cv::Mat;
	cv::Mat rtt_,rtt_2,foreground2;
        //rendertexturear->renderToTexture(rtt_);

        //std::cout << " write color1 " << std::endl;
        //cv::imwrite("color_10.png",color_1);

        rendertexturear->renderToTextureD(rtt_, color_1);
        *rtt = rtt_.clone();
	cv::cvtColor(rtt_,rtt_2,CV_BGR2RGB);
	cv::cvtColor(rgbddataprocessing->foreground,foreground2,CV_RGBA2RGB);
        dataio->listrtt.push_back(rtt);
	//dataio->listrttstress.push_back(rtt);

	timertt = ((double)getTickCount() - timertt)/getTickFrequency();
    cout << "Time RTT " << timertt << endl;
	
	double time3 = (double)getTickCount();
	
	/*if(useContour.getValue())
	extractSourceContour();	*/
    time3 = ((double)getTickCount() - time3)/getTickFrequency();
	cout << "Rigid extractcont " << time3 << endl;
	timeOverall += time3;
	timeSourceContour = time3;
	
	}
	else if (t%npasses == 1){
		
			rtt = new cv::Mat;
			cv::Mat rtt_,rtt_2,foreground2;	
	//setViewPoint1();
	
	rendertexturear->renderToTexture(rtt_);
	setViewPoint();

	*rtt = rtt_;
	cv::cvtColor(rtt_,rtt_2,CV_BGR2RGB);
	cv::cvtColor(rgbddataprocessing->foreground,foreground2,CV_RGBA2RGB);
				
	dataio->listrttstress.push_back(rtt);
	//dataio->listrttstressplast.push_back(rtt);
	
		
	}
	else{
		
	rtt = new cv::Mat;
	cv::Mat rtt_,rtt_2,foreground2;
    rendertexturear->renderToTexture(rtt_);
	setViewPoint();

	*rtt = rtt_;
	cv::cvtColor(rtt_,rtt_2,CV_BGR2RGB);
	cv::cvtColor(rgbddataprocessing->foreground,foreground2,CV_RGBA2RGB);

	dataio->listrttstressplast.push_back(rtt);
		
	}

	}

	//computeTargetNormals();
	
	double time = (double)getTickCount();

    if(ks.getValue()==0) return;
	
    VecDeriv&        f = *_f.beginEdit();       //WDataRefVecDeriv f(_f);
    const VecCoord&  x = _x.getValue();			//RDataRefVecCoord x(_x);
    const VecDeriv&  v = _v.getValue();			//RDataRefVecDeriv v(_v);
    ReadAccessor< Data< VecCoord > > tn(targetNormals);
    ReadAccessor< Data< VecCoord > > tp(targetPositions);
	ReadAccessor< Data< VecCoord > > tcp(targetContourPositions);
    ReadAccessor< Data< VecCoord > > xcp(sourceContourPositions);
	ReadAccessor< Data< VecCoord > > ssn(closestpoint->sourceSurfaceNormalsM);
	ReadAccessor< Data< VecCoord > > sn(sourceNormals);
	
	//VecDeriv&        f1; //= *_f.beginEdit();       //WDataRefVecDeriv f(_f);
    //const VecCoord&  x1;//= _x.getValue();			//RDataRefVecCoord x(_x);
    //const VecDeriv&  v1; // = _v.getValue();
	
			if (t%niterations.getValue() == 0) {
				f_.resize(f.size());
				x_.resize(x.size());
				v_.resize(v.size());
				/*for (int jj = 0; jj < x1.size(); jj++)
				{
					f_[jj] = f1[jj];					
					x_[jj] = x1[jj];
					v_[jj] = v1[jj];
				}*/
			}
	
	if (useCCD.getValue() && t >= niterations.getValue()){
		//updateCCD();
		
	//std::cout << " ok updateccd" << pointsCCD.size() << std::endl;
	//std::cout << " pointsCCD size " << pointsCCD.size() << " " << pointsCCD[0].X << std::endl;
	//ccd.local_statistics(pointsCCD,color);
	//ccd.refine_parameters(pointsCCD,color);
	
	double timeccd0 = (double)getTickCount();
    //ccd.local_statistics(pointsCCD,color);
	ccd.local_statistics_all(pointsCCD,color);
	pointsCCDmin = ccd.pointsccdmin;
	CCDPointsTo3D();
	
	double timeccd1 = ((double)getTickCount() - timeccd0)/getTickFrequency();
	
	}

    const vector<Spring>& s = this->springs.getValue();
    this->dfdx.resize(s.size());
    this->closestPos.resize(s.size());
	
	dfdx1.resize(s.size());
	
		   std::cout <<" ok ok 1  " << endl;

	
	//closestpoint->updateClosestPointsGt();

        if (useVisible.getValue() && t)
	closestpoint->sourceVisiblePositions.setValue(sourceVisiblePositions.getValue());
	
	closestpoint->timer = t;
	closestpoint->targetPositions.setValue(targetPositions.getValue());
	closestpoint->sourceSurfacePositions.setValue(sourceSurfacePositions.getValue());
    closestpoint->sourceBorder = meshprocessing->sourceBorder;
	closestpoint->targetBorder = rgbddataprocessing->targetBorder;

        std::cout <<" ok ok 2  "<< (sourceVisiblePositions.getValue()).size() << endl;
	
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
		//closestpoint->updateClosestPointsSoft();
	}

	
	/*closestSource = closestpoint->getClosestSource();
	closestTarget = closestpoint->getClosestTarget();
	
	sourceIgnored = closestpoint->getSourceIgnored();
	targetIgnored = closestpoint->getTargetIgnored();*/
	indices = closestpoint->getIndices();
			
    m_potentialEnergy = 0;

    // get attraction/ projection factors
    Real attrF=(Real) blendingFactor.getValue();
    if(attrF<(Real)0.) attrF=(Real)0.;
    if(attrF>(Real)1.) attrF=(Real)1.;
    Real projF=((Real)1.-attrF);
	
	std::cout << " tp size " << tp.size() << std::endl;
	

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
			int kkt=0;
			if (useContour.getValue()){
					if (t >2 ) 
					{
					for (unsigned int i=0; i<tp.size(); i++) 
					{					
						if(!closestpoint->targetIgnored[i])	
						cnt[indicesVisible[closestpoint->closestTarget[i].begin()->second]]++;
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
												
						//std::cout << " ind " << indicesVisible[closestpoint->closestTarget[i].begin()->second] << " " << closestpoint->closestTarget[i].begin()->second << std::endl;
							if(!closestpoint->targetIgnored[i])// && !targetBackground[i])
							cnt[indicesVisible[closestpoint->closestTarget[i].begin()->second]]++;
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

				std::cout << " tp size1 " << tp.size() << std::endl;

        // compute targetpos = projF*closestto + attrF* sum closestfrom / count
		
        // projection to point or plane
		
		double error = 0;
		int nerror = 0;
				
		int ivis=0;
		int kk = 0;
		unsigned int id;
		//if (t%niterations.getValue() == 0) 
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
						if(!meshprocessing->sourceBorder[i])						
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

				std::cout << " tp size2 " << tp.size() << std::endl;
						
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
               // if(/*!closestpoint->sourceIgnored[i] &&*/ sourceVisible[i])
					{
						
					if (sourceVisible[i]){
					if(!meshprocessing->sourceBorder[i])
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
                        std::cout << " tp size2 " << tcp.size() << " " << id << " " << kk << " " << indices.size()<< std::endl;

						closestPos[i]=tcp[id]*projF;
						
						if(!cnt[i]) closestPos[i]+=x[i]*attrF;
                        kk++;
                        std::cout << " tp size20 " << id << std::endl;

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
               // if(/*!closestpoint->sourceIgnored[i] &&*/ sourceVisible[i])
					{
						
					if (sourceVisible[i]){
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
               // if(/*!closestpoint->sourceIgnored[i] &&*/ sourceVisible[i])
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
			
		double errorICP = 0;
		// attraction

				std::cout << " tp size3 " << tp.size() << std::endl;
																
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
					//if (meshprocessing->sourceBorder[id] && rgbddataprocessing->targetBorder[i])
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
					unsigned int id1 = indicesVisible[id];
					closestPos[id1]+=tp[i]*attrF/(Real)cnt[id1];
					}
					

					}

					}
				//errorICP = determineErrorICP();
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
					id1 = indicesVisible[id];
					/*else {
					id1 = indicesTarget[kkt];
					kkt++;	
					}*/
					//if (meshprocessing->sourceBorder[id1] && rgbddataprocessing->targetBorder[i])
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
                //errorMatching = (double)error/nerror;
				//if ((t-46)%50==0)

        //errorfunction.setValue(errorMatching);
            //std::cout << " error tracking " << error << " norm " << (double)error/nerror << std::endl;
    }
		//if(useMassSpring.getValue()) correspMassSpringFile.open ("correspMassSpring.txt");

	std::cout << " ok 0 ok 0 " << std::endl;

	if(useGroundTruth.getValue())
		{
		double error_gt = 0;
		double error_vm_gt = 0;
		double error_ps_gt = 0;
		double error_es_gt = 0;
		double error_ts_gt = 0;


		double error_esn_gt = 0;
		double error_psn_gt = 0;
		double error_tsn_gt = 0;

		unsigned int indGt;
		int indS;
		ReadAccessor< Data< VecCoord > > tpGt(targetGtPositions);
		ReadAccessor< Data<helper::vector<Real> > > vM(vonMisesStress);
		ReadAccessor< Data<helper::vector<Real> > > eS(elasticStrainsN);
		ReadAccessor< Data<helper::vector<Real> > > pS(plasticStrainsN);
		ReadAccessor< Data<helper::vector<Real> > > tS(totalStrainsN);
		
		ReadAccessor< Data<helper::vector<Real> > > eSNode(elasticStrainsPerNode);
		ReadAccessor< Data<helper::vector<Real> > > pSNode(plasticStrainsPerNode);
		ReadAccessor< Data<helper::vector<Real> > > tSNode(totalStrainsPerNode);

		for (unsigned int i=0; i<tpGt.size(); i++)
		error_gt += this->computeError(x[i],tpGt[i]);
		
		double meanvm = 0;
		double meanps = 0;
		double meants = 0;
		double meanes = 0;

		
		double meanesn = 0;
		double meanpsn = 0;
		double meantsn = 0;
				
		for (unsigned int i=0; i<vM.size(); i++)
		{
			//if ((vonmisesstressGt[i]) > 0 )
				{
			error_vm_gt += abs((sqrt(vM[i]) /*- sqrt(vonmisesstressGt[i])*/));///sqrt(vonmisesstressGt[i]);
			error_es_gt += abs(sqrt(eS[i]) /*- sqrt(elasticstrainsGt[i])*/);
			error_ps_gt += abs(sqrt(pS[i]) /*- sqrt(plasticstrainsGt[i])*/);
			error_ts_gt += abs(sqrt(tS[i]) /*- sqrt(totalstrainsGt[i])*/);
			if (vonmisesstressGt.size()>0){
			meanvm += sqrt(vonmisesstressGt[i]);
			meanes += sqrt(elasticstrainsGt[i]);
			meanps += sqrt(plasticstrainsGt[i]);
			meants += sqrt(totalstrainsGt[i]);}

			}
		}
				
		for (unsigned int i=0; i<eSNode.size(); i++)
		{
			//if ((elasticstrainsnodeGt[i]) > 0 )
				{
			error_esn_gt += abs((sqrt(eSNode[i])/* - sqrt(elasticstrainsnodeGt[i])*/));///sqrt(vonmisesstressGt[i]);
			}
			//if ((plasticstrainsnodeGt[i]) > 0 )
				{
			error_psn_gt += abs(sqrt(pSNode[i]) /*- sqrt(plasticstrainsnodeGt[i])*/);///sqrt(vonmisesstressGt[i]);
			}
			error_tsn_gt += abs(sqrt(tSNode[i]) /*- sqrt(plasticstrainsnodeGt[i])*/);///sqrt(vonmisesstressGt[i]);

			if (elasticstrainsnodeGt.size()>0){
			meanesn += sqrt(elasticstrainsnodeGt[i]);
			meanpsn += sqrt(plasticstrainsnodeGt[i]);
			meantsn += sqrt(totalstrainsnodeGt[i]);
			}

		}

		std::cout << " Error Gt " << error_gt << " vm size "<< vM.size() << " esnode size " << eSNode.size() << std::endl;
		errorGroundTruth = (double)error_gt/tpGt.size();
		errorGroundTruth_vM = (double)((error_vm_gt)/vM.size());
		meanGt_vM = (double)meanvm/vM.size();
		errorGroundTruth_eS = (double)((error_es_gt)/vM.size());
		meanGt_eS = (double)meanes/vM.size();
		errorGroundTruth_pS = (double)((error_ps_gt)/vM.size());
		meanGt_pS = (double)meanps/vM.size();
		errorGroundTruth_tS = (double)((error_ts_gt)/vM.size());
		meanGt_tS = (double)meants/vM.size();

		errorGroundTruth_eSN = (double)((error_esn_gt)/eSNode.size());
		errorGroundTruth_pSN = (double)((error_psn_gt)/eSNode.size());
		errorGroundTruth_tSN = (double)((error_tsn_gt)/eSNode.size());
		meanGt_eSN = (double)meanesn/eSNode.size();
		meanGt_pSN = (double)meanpsn/eSNode.size();
		meanGt_tSN = (double)meantsn/eSNode.size();

		}
	
	ind = 0;
	int kc = 0;
	
	VecCoord  targetKLTPos;
	if (useKLTPoints.getValue() && t%niterations.getValue() == 0 )
	targetKLTPos = targetKLTPositions.getValue();
	
	VecCoord  targetCCDPos;
	if (useCCD.getValue() && t%niterations.getValue() == 0 )
	targetCCDPos = targetCCDPositions.getValue();
	
		 	std::cout << " Error " << std::endl;	
	
    for (unsigned int i=0; i<s.size(); i++)
    {
        //serr<<"addForce() between "<<springs[i].m1<<" and "<<closestPos[springs[i].m1]<<sendl;
				if (t > 2*niterations.getValue() && t%(niterations.getValue()) == 0)	
					{ 
					
                                        if( !useContour.getValue())
					this->addSpringForce(m_potentialEnergy,f,x,v, i, s[i]);
                                        else this->addSpringForceWeight(m_potentialEnergy,f,x,v, i, s[i]);
									
					if (useContour.getValue())
						{
						if (useCCD.getValue() && t > niterations.getValue()*3 && meshprocessing->sourceBorder[i]) 
							{//this->addSpringForceColorContour(m_potentialEnergy,f,x,v, kc, s[i]);
			//this->addSpringForceKLTA(m_potentialEnergy,f,x,v, targetCCDPos[kc], i, s[i], 1);
										kc++;
							}
						}
				
					}
    }
	
	 Vector3 coefs;
	 int index, id;
	 float xp0, yp0;
	 
	 cv::Mat distimg = rgbddataprocessing->seg.distImage;
	
	 if (useKLTPoints.getValue() && t >= 3 && t%(niterations.getValue()) == 0){
	 for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
		 //std::cout << " k " << k << std::endl;
					/*std::cout << " index " << index << std::endl;
					std::cout << " triangleindex " << triangles[index][0] << std::endl;
					std::cout << " targetKLTPos " << triangles[index][0] << std::endl;*/
					
			int x_u = (int)(x[triangles[index][0]][0]*rgbIntrinsicMatrix(0,0)/x[triangles[index][0]][2] + rgbIntrinsicMatrix(0,2));
			int x_v = (int)(x[triangles[index][0]][1]*rgbIntrinsicMatrix(1,1)/x[triangles[index][0]][2] + rgbIntrinsicMatrix(1,2));
	        //tracker.getFeature(k, id, xp0, yp0);
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
				
				
            if( depthValue < 1 && depthValue > 0 && avalue > 2/* && foreground.at<Vec4b>(yp0,xp0)[3] > 0*/){
			//std::cout << " x_u " <<targetKLTPos[id][2] << " kltfeat " << x[triangles[index][0]][2] << std::endl;
			//if (!meshprocessing->sourceBorder[triangles[index][0]]) this->addSpringForceKLTA(m_potentialEnergy,f,x,v, targetKLTPos[id], triangles[index][0], s[triangles[index][0]], coefs[0]);
			//if (!meshprocessing->sourceBorder[triangles[index][1]]) this->addSpringForceKLTA(m_potentialEnergy,f,x,v, targetKLTPos[id], triangles[index][1], s[triangles[index][1]], coefs[1]);
			//if (!meshprocessing->sourceBorder[triangles[index][2]]) this->addSpringForceKLTA(m_potentialEnergy,f,x,v, targetKLTPos[id], triangles[index][2], s[triangles[index][2]], coefs[2]);
			}
			
            }
			//getchar();
			 }
			//vpDisplay::getClick(vpI) ;
	
    _f.endEdit();
	
	time = ((double)getTickCount() - time)/getTickFrequency();
    cout << "Time external force addforce " << time << endl;
	timeAddforce += time;
	timeOverall += time;
	
	timei = (double)getTickCount();
	timeTotal =((double)getTickCount() - timeT)/getTickFrequency();
	cout << "Time total " << timeTotal << endl;
		
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
void RegistrationForceFieldCam<DataTypes>::addSpringForceKLT(double& potentialEnergy, VecDeriv& f, const  VecCoord& p,const VecDeriv& v, Coord& KLTtarget, int i, const Spring& spring, double coef)
{
    int a = spring.m1;
    Coord u = KLTtarget-p[a];
	//std::cout << " u " << u.norm() << std::endl;
    Real d = u.norm();
    if( d>1.0e-4 )
    {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        potentialEnergy += coef*elongation * elongation * spring.ks / 2;
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
            for( int k=0; k<N; ++k )
            {
                m[j][k] = 0;
            }
        }
    }
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::addSpringForceKLTA(double& potentialEnergy, VecDeriv& f, const  VecCoord& p,const VecDeriv& v, Coord& KLTtarget, int i, const Spring& spring, double coef)
{
    int a = spring.m1;
    Coord u = KLTtarget-p[a];
	//std::cout << " u " << p[a][2] << " t " << KLTtarget[2] << std::endl;
    Real d = u.norm();
    if( d>1.0e-4 )
    {
        Real inverseLength = 1.0f/d;
        u *= inverseLength;
        Real elongation = (Real)d;
        potentialEnergy += coef*elongation * elongation * spring.ks / 2;
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


template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::addStoredSpringForce(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
{
    int a = spring.m1;

        f[a]=f_[a];
        Mat& m = this->dfdx[i];
		m = dfdx1[i];
}

template <class DataTypes>
void RegistrationForceFieldCam<DataTypes>::addSpringForceColorContour(double& potentialEnergy, VecDeriv& f,const  VecCoord& p,const VecDeriv& v,int i, const Spring& spring)
{
    int a = spring.m1;
	Real alpha0 = alphaIntensity.getValue();
	Eigen::Matrix<float,3,1> displacement;
	//displacement = ccd.displacements[i];
	Coord disp;
	
	//std::cout << " force contour " << i << std::endl;
	
	disp[0] = ccd.displacements[i].x;
	disp[1] = ccd.displacements[i].y;
	disp[2] = ccd.displacements[i].z;
	/*disp[0] = 0;
	disp[1] = 0;
	disp[2] = 0;*/
	//std::cout << " jacob mat " << IntensityCurrent[i][0] << std::endl;
	//std::cout << " error intensity " << errorIntensity[0] << std::endl;
	Coord u0 = this->closestPos[i]-p[a];
	//std::cout << " displacement " << jacobian_matrix[0] << " " << jacobian_matrix[1] << " " << jacobian_matrix[2] << std::endl;
	//std::cout << " closest point " << disp[0] << " " << disp[1] << " " << disp[2] << std::endl;
	
	//getchar();
	
    Coord u =/* this->closestPos[i]-p[a]*/ -alpha0*disp;
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
            for( int k=0; k<N; ++k ) m[j][k] += ((Real)spring.ks-tgt) * u[j] * u[k];
            m[j][j] += tgt;
        }
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
void RegistrationForceFieldCam<DataTypes>::normalizeWeights()
{
            const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
		double totalweights = 0;
		combinedWeights.resize(0);
	for (int i = 0; i < x.size(); i++)	
	for (int k = 0; k < targetPositions.getValue().size(); k++)
			if(closestpoint->closestTarget[k].begin()->second == i || closestpoint->closestSource[i].begin()->second == k)
			{ 
			combinedWeights.push_back((double)sourceWeights[i]*targetWeights[k]);
			totalweights += (double)sourceWeights[i]*targetWeights[k];
			}
			
			std::cout << " combined weights " << combinedWeights.size() << std::endl;
	for (int i = 0; i < combinedWeights.size(); i++)
	{
          combinedWeights[i] *= (double)(combinedWeights.size()/totalweights);
		  //std::cout << " combined weights " << combinedWeights[i] << std::endl;
	}			
	
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
		
		for (int k = 0; k < targetPositions.getValue().size(); k++){
			if(closestpoint->closestTarget[k].begin()->second == i || closestpoint->closestSource[i].begin()->second == k)
			{
				if(meshprocessing->sourceBorder[i])
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
			
		if (meshprocessing->sourceBorder[i]) stiffweight*=1;
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

/*template<class DataTypes, class DepthTypes>
void RegistrationForceFieldCam<DataTypes, DepthTypes>::addKToMatrix(const core::MechanicalParams* mparams,const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    if(ks.getValue()==0) return;

    double kFact = mparams->kFactor();

    sofa::core::behavior::MultiMatrixAccessor::MatrixRef mat = matrix->getMatrix(this->mstate);
	
    if (!mat) return;
    const vector<Spring >& ss = this->springs.getValue();
    const unsigned int n = ss.size() < this->dfdx.size() ? ss.size() : this->dfdx.size();
    for (unsigned int e=0; e<n; e++)
    {
        const Spring& s = ss[e];
        unsigned p1 = mat.offset+Deriv::total_size*s.m1;
        const Mat& m = this->dfdx[e];
        for(int i=0; i<N; i++)
            for (int j=0; j<N; j++)
            {
                Real k = (Real)(m[i][j]*kFact);
                mat.matrix->add(p1+i,p1+j, -k);
            }
    }
}*/

            
template<class DataTypes>
void RegistrationForceFieldCam<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

        int t = (int)this->getContext()->getTime();

        double timef = 0;
        timef = (double)getTickCount();
        //if (t > 0 && t%niterations.getValue() == niterations.getValue()-1)
        timeOverall += (timef - timei)/getTickFrequency();
        timeResolution += (timef - timei)/getTickFrequency();

        {
    cout <<" t " << t << " Time overall draw " << timei << " " << timef << " " << (timef - timei)/getTickFrequency() << endl;
        }
        /*
        //double timef = 0;
        timeii = timef;
        timef = (double)getTickCount();
        //if (t > 0 && t%niterations.getValue() == niterations.getValue()-1)
        timeOverall += (timef - timei)/getTickFrequency();
        timeResolution += (timef - timei)/getTickFrequency();
        {
    cout <<" t " << t << " Time overall draw " << timei << " " << timef << " " << (timef - timeii)/getTickFrequency() << endl;
        }

        if (t > 0 && t%niterations.getValue() == 0)//niterations.getValue()-1)
        {
        std::cout << " t " << t << " " << t/niterations.getValue() <<  " time overall " << timeOverall << std::endl;
        timeFile << t/niterations.getValue();
        timeFile << "\t";
        timeFile << timeOverall;
        timeFile << "\t";
        timeFile << timeSeg;
        timeFile << "\t";
        timeFile << timeRigid;
        timeFile << "\t";
        timeFile << timeSourceContour;
        timeFile << "\t";
        timeFile << timeAddforce;
        timeFile << "\t";
        timeFile << timeResolution;
        timeFile << "\n";

        }

        if (t == 1480)
        {
                timeFile.close();
        }*/

                ReadAccessor< Data< VecCoord > > x(this->getMState()->read(core::ConstVecCoordId::position()));

        const VecCoord&  p0 = this->mstate->read(core::ConstVecCoordId::position())->getValue();

        sourcePositions.setValue(p0);

        if (t > 0 && t%niterations.getValue() == niterations.getValue()-1)
        {
        std::cout << " t " << t << " " << t/niterations.getValue() <<  " time overall " << timeOverall << std::endl;
        timeFile << t/niterations.getValue();
        timeFile << "\t";
        timeFile << timeOverall;
        timeFile << "\t";
        timeFile << timeSeg;
        timeFile << "\t";
        timeFile << timeRigid;
        timeFile << "\t";
        timeFile << timeSourceContour;
        timeFile << "\t";
        timeFile << timeAddforce;
        timeFile << "\t";
        timeFile << timeResolution;
        timeFile << "\n";

        errorFile << t/niterations.getValue();
        errorFile << "\t";
        errorFile << errorMatching;
        errorFile << "\n";

        if(useGroundTruth.getValue())
        {
        double errorICP = 0; //determineErrorICP();
                        double error1 = 0;
                int nerror1 = 0;

                VecCoord  tcp = targetContourPositions.getValue();

                                                if (t >2)
                                                {
                                //updateClosestPointsVisible();
                                //updateClosestPointsVisibleContours();
                int kc = 0;
                for(int i=0;i<(int)x.size();i++)
        {

                        //int id = indicesVisible[i];
                        bool useP = true;
            /*if (useContour.getValue())
                        {if(!meshprocessing->sourceBorder[i]) useP = false;}// && t%niterations.getValue() == 0)
                        else
                        {
                        if (useP)
            {
                                double distmin = 1000;
                                double distmin1;
                                double dist, dist1,dist2;
                                int kmin2,kmin1;

                                for (int k = 0; k < tcp.size(); k++)
                                {
                                        dist = (x[i][0] - tcp[k][0])*(x[i][0] - tcp[k][0]) + (x[i][1] - tcp[k][1])*(x[i][1] - tcp[k][1]) + (x[i][2] - tcp[k][2])*(x[i][2] - tcp[k][2]);

                                        //if (dist < distmin)
                                        if (dist < distmin)
                                        {
                                                distmin = dist;
                                                kmin1 = k;
                                        }
                                }

                                error1 += computeError(x[i], tcp[kmin1]);
                                nerror1++;
                                kc++;
            }

            }*/
        }

                                                }

        std::cout << " error1 " << (double)error1/nerror1 << " error matching " << errorMatching << std::endl;
        errorGtFile << t/niterations.getValue();
        errorGtFile << "\t";
        errorGtFile << (double)error1/nerror1;
        errorGtFile << "\t";
        errorGtFile << errorGroundTruth;
        errorGtFile << "\t";
        errorGtFile << errorGroundTruth_vM;
        errorGtFile << "\t";
        errorGtFile << meanGt_vM;
        errorGtFile << "\t";
        errorGtFile << errorGroundTruth_eS;
        errorGtFile << "\t";
        errorGtFile << meanGt_eS;
        errorGtFile << "\t";
        errorGtFile << errorGroundTruth_pS;
        errorGtFile << "\t";
        errorGtFile << meanGt_pS;
        errorGtFile << "\t";
        errorGtFile << errorGroundTruth_tS;
        errorGtFile << "\t";
        errorGtFile << meanGt_tS;
        errorGtFile << "\t";
        errorGtFile << errorGroundTruth_eSN;
        errorGtFile << "\t";
        errorGtFile << meanGt_eSN;
        errorGtFile << "\t";
        errorGtFile << errorGroundTruth_pSN;
        errorGtFile << "\t";
        errorGtFile << meanGt_pSN;
        errorGtFile << "\t";
        errorGtFile << errorGroundTruth_tSN;
        errorGtFile << "\t";
        errorGtFile << meanGt_tSN;

                for (unsigned int i=0; i < x.size(); i++)
                {
        /*errorGtFile << x[i][0];
        errorGtFile << "\t";
        errorGtFile << x[i][1];
        errorGtFile << "\t";
        errorGtFile << x[i][2];
        errorGtFile << "\n";*/
                }
                        errorGtFile << "\n";
        }

        vpHomogeneousMatrix eMc ;

  /*eMc[0][0] = 0;
  eMc[0][1] = -1;
  eMc[0][2] = 0;
  eMc[0][3] = h;

  eMc[1][0] = 1;
  eMc[1][1] = 0;
  eMc[1][2] = 0;
  eMc[1][3] = -L;

  eMc[2][0] = 0;
  eMc[2][1] = 0;
  eMc[2][2] = 1;
  eMc[2][3] = 0;

  eMc[3][0] = 0;
  eMc[3][1] = 0;
  eMc[3][2] = 0;
  eMc[3][3] = 1;*/

  /*eMc[0][0] = 0.08135292269;
  eMc[0][1] = -0.5450835618;
  eMc[0][2] = -0.8344253188;
  eMc[0][3] = 0.1057791273;

  eMc[1][0] = 0.9929451685;
  eMc[1][1] = -0.02813634394;
  eMc[1][2] = 0.11518784;
  eMc[1][3] = 0.02447798141;

  eMc[2][0] = -0.08626467586;
  eMc[2][1] = -0.8379094562;
  eMc[2][2] = 0.5389491153;
  eMc[2][3] = 0.3080637527;*/

  eMc[0][0] =  1;
  eMc[0][1] = 0;
  eMc[0][2] = 0;
  eMc[0][3] = 0.005947680623;

  eMc[1][0] = 0;
  eMc[1][1] =1;
  eMc[1][2] = 0;
  eMc[1][3] = -0.2665557298 ;

  eMc[2][0] = 0;
  eMc[2][1] = 0;
  eMc[2][2] = 1;
  eMc[2][3] = 0.02643054428;

  /*0.05583840319  -0.9157786066  0.3977833788  0.2560252907
0.9949306037  0.01766200412  -0.09900074455  -0.02609790375
0.08363711222  0.4012949007  0.91212238  -0.04549550763  */

  vpColVector xc,xe;
  xc.resize(4);
  xe.resize(4);

                xc[0] = x[4][0];
                xc[1] = x[4][1];
                xc[2] = x[4][2];
                xc[3] = 1;
                xe = eMc*xc;
                center_pos.x = xe[0];
                center_pos.y = xe[1];
                center_pos.z = xe[2];

                        int ns = write(sock_sender, &center_pos, sizeof(center_pos));
                        if (ns < 0) std::cout<<"Errore Sendto"<< std::endl;
        }

        if (t == nimages.getValue())
        {
                timeFile.close();
                errorFile.close();
                errorGtFile.close();
        }

          glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT | GL_ENABLE_BIT);
          glDisable( GL_LIGHTING);
        if(drawColorMap.getValue())
    {
        ReadAccessor< Data< VecCoord > > xtarget(targetPositions);

                glEnable(GL_BLEND);
        glPointSize( 10);
        glBegin( GL_POINTS);

        VecCoord  targetKLTPos;
        targetKLTPos = targetKLTPositions.getValue();

        VecCoord  targetCCDPos;
        targetCCDPos = targetCCDPositions.getValue();

                                sofa::helper::vector< tri > triangles;
                                triangles = sourceTriangles.getValue();
                                if (t > 5)
                                {
                if ( drawTarget.getValue())
                for (unsigned int i=0; i< targetPositions.getValue().size(); i++)
          {
            //sofa::helper::gl::Color::setHSVA(dists[i]*240./max,1.,.8,1.);
                        //if(targetBackground[i])
                        sofa::helper::gl::Color::setHSVA(80,1.,.8,1.);
            if (useContour.getValue() && drawContour.getValue())
                        {if(rgbddataprocessing->targetBorder[i]) sofa::helper::gl::Color::setHSVA(40,1.,.8,1.);}

                        glVertex3d(xtarget[i][0],xtarget[i][1],xtarget[i][2]);

                  }
                                std::vector< Vector3 > points;
                                const vector<Spring>& springs = this->springs.getValue();

                                int ivis = 0;

        /*for (unsigned int i=0; i<springs.size(); i++){
            if(closestpoint->sourceIgnored[ivis] ){

                        bool dispP = true;
                        if (useContour.getValue() && drawContour.getValue())
                        {if(!meshprocessing->sourceBorder[i]) dispP = false;}// && t%niterations.getValue() == 0)
                Vector3 point1 = DataTypes::getCPos(x[springs[i].m1]);
                Vector3 point2 = DataTypes::getCPos(this->closestPos[i]);
                                //std::cout << " pt " << point2[0] << " " << point2[1] << " " << point2[2] << std::endl;
                                if (dispP){
                points.push_back(point1);
                points.push_back(point2);}
                }
                ivis++;
                }*/
                                         Vector3 coefs;
         int index, id;
         float xp0, yp0;
                                        if (t > 50 && useKLTPoints.getValue())
                                 for (unsigned int k = 0; k < static_cast<unsigned int>(tracker.getNbFeatures()); k++){
                 //std::cout << " k " << k << std::endl;
                                        coefs = mappingkltcoef[k];
                                        index = mappingkltind[k];
                                //tracker.getFeature(k, id, xp0, yp0);
                        float depthValue = (float)depth.at<float>(2*yp0,2*xp0);
            if(1-(coefs[0] + coefs[1]) > 0 && depthValue < 1 && depthValue > 0){
                        //std::cout << " x_u " << 1-(coefs[0] + coefs[1]) << " " << coefs[0] << " kltfeat " << coefs[1] << " " << targetKLTPos[k][2] << std::endl;
                Vector3 point1 = DataTypes::getCPos(x[springs[triangles[index][2]].m1]);
                Vector3 point2 = DataTypes::getCPos(targetKLTPos[k]);
                                //std::cout << " pt " << point2[0] << " " << point2[1] << " " << point2[2] << std::endl;
                points.push_back(point1);
                points.push_back(point2);
                        //this->addSpringForceKLTA(m_potentialEnergy,f,x,v, targetKLTPos[k], triangles[index][1], s[triangles[index][1]],coefs[0]);
                        //this->addSpringForceKLTA(m_potentialEnergy,f,x,v, targetKLTPos[k], triangles[index][2], s[triangles[index][2]],coefs[1]);
                        }

            }

                                if (t > 20 && useCCD.getValue())
                                {
                                        int kc = 0;

                                        for (unsigned int i=0; i<x.size(); i++)
                                        {
                                        bool dispPV = true;
                                        if (drawContour.getValue())
                                                if (meshprocessing->sourceBorder[i])
                                                {
                                                        //std::cout << " pt 0 "<< kc << " " << targetCCDPositions.getValue()[kc][0] << " " << targetCCDPositions.getValue()[kc][1] << " " << targetCCDPositions.getValue()[kc][2] << std::endl;
                                Vector3 point1 = DataTypes::getCPos(x[springs[i].m1]);
                Vector3 point2 = DataTypes::getCPos(targetCCDPos[kc]);
                                //std::cout << " pt 0 " << i << " " << kc << " " << targetCCDPos.size() << std::endl;
                                //std::cout << " pt 0 " << point1[0] << " " << point1[1] << " " << point1[2] << std::endl;
                                //std::cout << " pt 1 " << point2[0] << " " << point2[1] << " " << point2[2] << std::endl;

                //points.push_back(point1);
                                //points.push_back(point2);
                                //if (kc == 30)
                                        {
                                sofa::helper::gl::Color::setHSVA(140,1.,.8,1.5);
                                glVertex3d(targetCCDPos[kc][0],targetCCDPos[kc][1],targetCCDPos[kc][2]);
                                        sofa::helper::gl::Color::setHSVA(200,1.,.8,1.5);}
                                //glVertex3d(x[springs[i].m1][0],x[springs[i].m1][1],x[springs[i].m1][2]);}

                                kc++;
                                }

                                        }
                                }

        const Vec<4,float> c(1,0,0,1);
        //if (showArrowSize.getValue()==0 || drawMode.getValue() == 0)
                        //vparams->drawTool()->drawLines(points, 1, c);
        if (drawMode.getValue() == 1)
                        for (unsigned int i=0;i<points.size()/2;++i) vparams->drawTool()->drawCylinder(points[2*i+1], points[2*i], showArrowSize.getValue(), c);
        else if (drawMode.getValue() == 2)
                        for (unsigned int i=0;i<points.size()/2;++i) vparams->drawTool()->drawArrow(points[2*i+1], points[2*i], 0.002, c);


                                        //}
        //else serr << "No proper drawing mode found!" << sendl;
                  int kcp = 0;

                  if (drawSource.getValue())
                        for (unsigned int i=0; i<x.size(); i++)
                                        {
                                        bool dispPV = true;
                                        if (useVisible.getValue()){if(!sourceVisible[i]) dispPV = false;}
                                                sofa::helper::gl::Color::setHSVA(140,1.,.8,1.5);
                                        if (useContour.getValue() && drawContour.getValue())
                                                {if (meshprocessing->sourceBorder[i]) sofa::helper::gl::Color::setHSVA(200,1.,.8,1.5);}
                                                if (dispPV)
                                                glVertex3d(x[i][0],x[i][1],x[i][2]);
                                        }

                        }

                glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT | GL_ENABLE_BIT);
        glDisable( GL_LIGHTING);
                glBegin( GL_POINTS);

        glEnd();
        glPointSize( 1);

        glPopAttrib();
    }

}

}
}
} // namespace sofa

//#endif  /* SOFA_COMPONENT_INTERACTIONFORCEFIELD_RegistrationForceFieldCam_INL */


