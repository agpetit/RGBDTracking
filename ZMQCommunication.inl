#ifndef SOFA_RGBDTRACKING_ZMQCOMMUNICATION_INL
#define SOFA_RGBDTRACKING_ZMQCOMMUNICATION_INL

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/core/ObjectFactory.h>
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
#include "ImageConverter.h"
#include "ZMQCommunication.h"

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
        d_SubORPub(initData(&d_SubORPub, "SubORPub", "publisher = 0, subscriber = 1"))
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
        if (d_SubORPub.getValue())
        {

            sofa::simulation::Node::SPtr root = dynamic_cast<simulation::Node*>(this->getContext());
            typename sofa::core::objectmodel::ImageConverter<DataTypes,DepthTypes>::SPtr imconv;
            root->get(imconv);

            if(!((imconv->depth).empty()) && !((imconv->color).empty()))
            {

                color = imconv->color;
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
            }
        }
        else
        {
                zmq::message_t message1,messageImage;
                std::cout << " receive " << std::endl;
                cv::Mat img;

                bool status = m_sockSub->recv(&message1);
                bool status1 = m_sockSub->recv(&messageImage);

                //if(status1)
                {
                std::string rpl = std::string(static_cast<char*>(messageImage.data()), messageImage.size());

                const char *cstr = rpl.c_str();
                loadImage(img,cstr);

               // memcpy(img.data, message1.data(), imgSize);
                std::cout << " ok receive " << std::endl;
                cv::imwrite("img.png", img);

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



                //std::cout << " receive ok " << std::endl;
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


}  // namespace fusion
}  // namespace sofa
}

#endif
