#ifndef SOFA_RGBDTRACKING_ZMQCOMMUNICATION_CPP
#define SOFA_RGBDTRACKING_ZMQCOMMUNICATION_CPP

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/Mapping.inl>
#include <sofa/simulation/Simulation.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/gui/BaseGUI.h>
#include <sofa/gui/BaseViewer.h>
#include <sofa/gui/GUIManager.h>
#include <iostream>
#include <map>


#include <sofa/helper/gl/Color.h>
#include <SofaBaseVisual/BaseCamera.h>
#include <SofaBaseVisual/InteractiveCamera.h>
#include <sofa/core/behavior/ForceField.inl>
#include <pcl/keypoints/sift_keypoint.h>

#include "ZMQCommunication.h"
#include "ImageConverter.h"


namespace sofa
{
namespace core
{
namespace objectmodel
{


using namespace sofa::defaulttype;
using namespace helper;


template<class DataTypes>
ZMQCommunication<DataTypes>::ZMQCommunication()
    : Inherited()
    , d_host(initData(&d_host, std::string("127.0.0.1"), "host",
                      "hostname to connect to")),
      d_port(initData(&d_port, ushort(8888), "port", "port to connect to")),
      d_positions(initData(&d_positions, "positions", "3D Positions to send over the network")),
      d_normals(initData(&d_normals, "normals", "3D normals to send over the network")),
      d_SubORPub(initData(&d_SubORPub, "SubORPub", "publisher = 0, subscriber = 1")),
      displayBackgroundImage(initData(&displayBackgroundImage,false,"displayBackgroundImage"," ")),
      useSensor(initData(&useSensor,true,"useSensor","Use real data")),
      niterations(initData(&niterations,1,"niterations","Number of images"))

{
    this->f_listening.setValue(true);
}

template<class DataTypes>
ZMQCommunication<DataTypes>::~ZMQCommunication()
{
cleanup();
}


template<class DataTypes>
void ZMQCommunication<DataTypes>::init()
{
        if (d_SubORPub.getValue())
        {
                std::stringstream serialized;
                serialize(serialized);
                m_sockPub = new zmq::socket_t(m_context, ZMQ_PUB);
                m_sockPub->bind("tcp://127.0.0.1:6667");
        }
        else
        {


                m_sockSub = new zmq::socket_t(m_context, ZMQ_SUB);
                m_sockSub->setsockopt(ZMQ_SUBSCRIBE, "", 0);
        m_sockSub->connect("tcp://127.0.0.1:6667");
        }
}


template<class DataTypes>
void ZMQCommunication<DataTypes>::serialize(std::stringstream& s)
{
        //double t = (double)this->getContext()->getTime();
        // 3d positions


        for (int k=0; k < d_positions.getValue().size(); k++)
        {
                s << d_positions.getValue()[k];
                s << ";";
                //std::cout << " pos send " << d_positions.getValue()[k][0] << " " << d_positions.getValue()[k][1] << " " << d_positions.getValue()[k][2] << std::endl;
        }
        //std::cout << " ok serialize " << std::endl;
}


template<class DataTypes>
void ZMQCommunication<DataTypes>::desserialize(std::stringstream& s)
{
        sofa::helper::vector<Vec3d>  positions_;

        positions_.resize(0);
        int signe = 0;

        while (!s.eof())
        {

                int npoints = 0;

                while (s.peek()!=';')
                {
                        Vec3d point;
                        double x,y,z;
                        s >> x;
                        s >> y;
                        s >> z;
                        if(signe==1) x=-x;
                        point[0] = x;
                        point[1] = y;
                        point[2] = z;
                        npoints ++;
                        //std::cout << " pos receive " << x << " " << y << " " << z << std::endl;
                        positions_.push_back(point);
                }
                s.ignore();
                if(s.peek()=='-') signe = 1;
                else signe = 0;
                s.ignore();
        }

        d_positions.setValue(positions_);
}

std::string save( const cv::Mat & mat )
{
    std::ostringstream oss;
    boost::archive::text_oarchive toa( oss );
    toa << mat;

    return oss.str();
}

void loadImage( cv::Mat & mat, const char * data_str )
{
    std::stringstream ss;
    ss << data_str;

    boost::archive::text_iarchive tia( ss );
    tia >> mat;
}

template<class DataTypes>
void ZMQCommunication<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{
        if (dynamic_cast<simulation::AnimateBeginEvent*>(event)){

        int t = (int)this->getContext()->getTime();

        if (d_SubORPub.getValue())
        {
            if (t%niterations.getValue()==0)
            {

            sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
            typename sofa::core::objectmodel::ImageConverter<DataTypes,DepthTypes>::SPtr imconv;
            root->get(imconv);

            typename sofa::core::objectmodel::DataIO<DataTypes>::SPtr dataio;
            root->get(dataio);

            if (useSensor.getValue()){
                if(!((imconv->depth).empty()) && !((imconv->color).empty()))
                {
                    color = imconv->color;
                }
            }
            else {
                    color = dataio->color;

             }

            if(!(color.empty()))
            {

                std::stringstream serialized;
                this->serialize(serialized);

                zmq::message_t message(serialized.str().length());
                std::memcpy(message.data(), serialized.str().c_str(),
                      serialized.str().length());

                std::string serializedImage = save(color);

                // Send data here

                 zmq::message_t messageImage(serializedImage.length());
                 memcpy(messageImage.data(), serializedImage.c_str(), serializedImage.length());

                bool status = m_sockPub->send(message);
                bool status1 = m_sockPub->send(messageImage);


                if (!status) msg_error(getName() + "::update()") << "could not send message";
                if (!status1) msg_error(getName() + "::update()") << "could not send image";

            }
            }
        }
        else
        {
                zmq::message_t message1,messageImage;

                bool status = m_sockSub->recv(&message1);
                bool status1 = m_sockSub->recv(&messageImage);

                std::cout << " ok receive " << status << " status 1 " << status1 << std::endl;

                if(status1)
                {
                std::string rpl = std::string(static_cast<char*>(messageImage.data()), messageImage.size());

                std::cout << "messageimage size " << messageImage.size() << std::endl;

                const char *cstr = rpl.c_str();
                loadImage(color,cstr);

               // memcpy(img.data, message1.data(), imgSize);
            //  // Assign pixel value to img
            //  for (int i = 0;  i < img1.rows; i++) {
            //   for (int j = 0; j < img1.cols; j++) {
            //    img1.at<Vec4b>(i,j)[0] = img.at<uchar>(0,i*img1.cols+j);
            //    img1.at<Vec4b>(i,j)[1] = img.at<uchar>(0,i*img1.cols+j + 1);
            //    img1.at<Vec4b>(i,j)[2] = img.at<uchar>(0,i*img1.cols+j + 2);
            //    img1.at<Vec4b>(i,j)[3] = img.at<uchar>(0,i*img1.cols+j + 3);
            //    }
            //   }
            //  cv::imwrite("socketimage100.png", img1);
                }
                std::cout << " receive ok " << std::endl;
                char messageChar1[message1.size()];
                memcpy(&messageChar1, message1.data(), message1.size());

                std::stringstream stream1;
                unsigned int nbDataFieldReceived = 0;

                for(unsigned int i=0; i<message1.size(); i++)
                        stream1 << messageChar1[i];

                if (status)
                {
                        desserialize(stream1);
                }
                else
                        msg_error(getName() + "::update()") << "could not retrieve message";
        }
}
}


template<class DataTypes>
void ZMQCommunication<DataTypes>::cleanup()
{
        if (m_sockSub)
        {
                m_sockSub->close();
                delete m_sockSub;
                m_sockSub = NULL;
        }
        //if (m_sockPub->connected())
        {
                m_sockPub->close();
                delete m_sockPub;
                m_sockPub = NULL;
        }
}

template <class DataTypes>
void ZMQCommunication<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

    vparams->drawTool()->saveLastState();


    if (displayBackgroundImage.getValue())
    {
    GLfloat projectionMatrixData[16];
    glGetFloatv(GL_PROJECTION_MATRIX, projectionMatrixData);
    GLfloat modelviewMatrixData[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, modelviewMatrixData);

    cv::Mat colorrgb = color.clone();
    if (!color.empty())
    cv::cvtColor(color, colorrgb, CV_RGB2BGR);

    std::stringstream imageString;
    imageString.write((const char*)colorrgb.data, colorrgb.total()*colorrgb.elemSize());
    // PERSPECTIVE

    glMatrixMode(GL_PROJECTION);	//init the projection matrix
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, 1, 0, 1, -1, 1);  // orthogonal view
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // BACKGROUND TEXTURING
    //glDepthMask (GL_FALSE);		// disable the writing of zBuffer
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_TEXTURE_2D);	// enable the texture
    glDisable(GL_LIGHTING);		// disable the light

    glBindTexture(GL_TEXTURE_2D, 0);  // texture bind
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, colorrgb.cols, colorrgb.rows, 0, GL_RGB, GL_UNSIGNED_BYTE, imageString.str().c_str());

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);	// Linear Filtering

                                                                        // BACKGROUND DRAWING
                                                                        //glEnable(GL_DEPTH_TEST);

    glBegin(GL_QUADS); //we draw a quad on the entire screen (0,0 1,0 1,1 0,1)
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glTexCoord2f(0, 1);		glVertex2f(0, 0);
    glTexCoord2f(1, 1);		glVertex2f(1, 0);
    glTexCoord2f(1, 0);		glVertex2f(1, 1);
    glTexCoord2f(0, 0);		glVertex2f(0, 1);
    glEnd();

    //glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);		// enable light
    glDisable(GL_TEXTURE_2D);	// disable texture 2D
    glEnable(GL_DEPTH_TEST);
    //glDepthMask (GL_TRUE);		// enable zBuffer

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);


    vparams->drawTool()->restoreLastState();
    }

}


////////////////////////////////////////////    FACTORY    ////////////////////////////////////////////
using sofa::core::RegisterObject ;

// Registering the component
SOFA_DECL_CLASS(ZMQCommunication)

int ZMQCommunicationClass = RegisterObject("This component is used to build a communication between two simulations")

//#ifdef SOFA_DOUBLE
.add< ZMQCommunication<defaulttype::Vec3dTypes> >(true)
//#endif
//#ifdef SOFA_FLOAT
.add< ZMQCommunication<defaulttype::Vec3fTypes> >()
//#endif
;

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Force template specialization for the most common sofa floating point related type.
// This goes with the extern template declaration in the .h. Declaring extern template
// avoid the code generation of the template for each compilation unit.
// see: http://www.stroustrup.com/C++11FAQ.html#extern-templates
//#ifdef SOFA_DOUBLE
template class SOFA_RGBDTRACKING_API ZMQCommunication<defaulttype::Vec3dTypes>;
//#endif
//#ifdef SOFA_FLOAT
template class SOFA_RGBDTRACKING_API ZMQCommunication<defaulttype::Vec3fTypes>;
//#endif


//	err = p_helper_socket::socket_create_UDP_Sender(sock_sender, 6999, /*ip_recipient*/"192.168.3.29");

//	if(err < 0)
//	cout<<"errore sender"<<endl;

}  // namespace fusion
}  // namespace sofa
}

#endif
