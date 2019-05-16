#ifndef SOFA_RGBDTRACKING_ZMQCOMMUNICATION_H
#define SOFA_RGBDTRACKING_ZMQCOMMUNICATION_H

#include <RGBDTracking/config.h>
#include <boost/thread.hpp>
#include <SofaBaseTopology/TopologyData.h>
#include <SofaUserInteraction/Controller.h>
#include "DataIO.h"


#include <sofa/core/core.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/accessor.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <SofaBaseTopology/TopologyData.h>
#include <sofa/gui/BaseGUI.h>
#include <sofa/gui/BaseViewer.h>
#include <sofa/gui/GUIManager.h>
#include <image/ImageTypes.h>


#include <zmq.hpp>
#include <string>

#include <sofa/helper/OptionsGroup.h>
#include <sofa/helper/vectorData.h>

#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/simulation/AnimateEndEvent.h>

#include <opencv2/core.hpp>


#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <fstream>
#include "serialization.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


#include <sstream>

namespace sofa
{
namespace core
{

namespace objectmodel
{

using core::objectmodel::Event;
using core::objectmodel::BaseObjectDescription;
using std::map;
using std::string;
using sofa::helper::vectorData;
using helper::vector;
using namespace sofa::defaulttype;

template< class DataTypes >
class ZMQCommunication : public virtual sofa::core::objectmodel::BaseObject
{

  typedef sofa::defaulttype::Vec3d Vec3d;
  typedef sofa::defaulttype::Vec2d Vec2d;
  typedef sofa::defaulttype::Vec3f Vec3f;
  typedef defaulttype::ImageF DepthTypes;


 public:
  typedef BaseObject Inherited;
  SOFA_CLASS(SOFA_TEMPLATE(ZMQCommunication,DataTypes), Inherited);


  ZMQCommunication();
  virtual ~ZMQCommunication();

  sofa::Data<std::string> d_host;
  sofa::Data<ushort> d_port;

  sofa::Data<sofa::helper::vector<Vec3d> > d_positions;
  sofa::Data<sofa::helper::vector<Vec3d> > d_normals;

  Data<bool> useSensor;

  cv::Mat color;

  sofa::Data<bool> d_SubORPub;
  Data<bool> displayBackgroundImage;
  Data<int> niterations;

  ////////////////////////// Inherited from BaseObject ////////////////////
  virtual void init() override;
  //virtual void reinit() override;
  //virtual void reset() override;
  /// Parse the given description to assign values to this object's fields and potentially other parameters
  //virtual void parse(BaseObjectDescription *arg) override;
  /// Assign the field values stored in the given map of name -> value pairs
  //virtual void parseFields(const map<string,string*>& str) override;
  void handleEvent(sofa::core::objectmodel::Event *event);
  /////////////////////////////////////////////////////////////////////////

  void serialize(std::stringstream& s);
  void desserialize(std::stringstream& s);
  void cleanup();
  void draw(const core::visual::VisualParams* vparams) ;


 private:
  zmq::context_t m_context{1};
  zmq::socket_t* m_sockSub;
  zmq::socket_t* m_sockPub;

};


}  // namespace controller
}  // namespace component
}  // namespace sofa

#endif  // SOFA_ZMQCOMMUNICATION_H
