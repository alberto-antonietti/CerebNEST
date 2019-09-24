#ifndef CORTEX_NEURON_H
#define CORTEX_NEURON_H

#include <fstream>

// Includes from librandom:
#include "poisson_randomdev.h"

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "mynames.h"

namespace mynest
{
class cortex_neuron : public nest::Archiving_Node
{

public:
  cortex_neuron();
  cortex_neuron( const cortex_neuron& );
  ~cortex_neuron();

  using nest::Node::handle;
  using nest::Node::handles_test_event;

  nest::port send_test_event( nest::Node&, nest::rport, nest::synindex, bool );
  nest::port handles_test_event( nest::SpikeEvent&, nest::rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  void init_state_( const nest::Node& proto );
  void init_buffers_();

  void calibrate();
  void update( nest::Time const&, const long, const long );
  void handle( nest::SpikeEvent& );

  struct Buffers_
  {
    std::map<long, int> out_spikes_;
    std::vector<int> in_spikes_;
  };

  struct Parameters_
  {
    long trial_length_;
    long joint_id_;
    long fiber_id_;
    long fibers_per_joint_;
    double rbf_sdev_;
    double baseline_rate_;
    double gain_rate_;
    bool to_file_;

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dicitonary
  };

  struct Variables_
  {
    double in_rate_;
    librandom::PoissonRandomDev poisson_dev_; //!< Random deviate generator
    std::ofstream out_file_; //!< OutputFile
  };

  Buffers_ B_;
  Parameters_ P_;
  Variables_ V_;
};

inline nest::port
cortex_neuron::send_test_event( nest::Node& target, nest::rport receptor_type, nest::synindex, bool )
{
  nest::SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline nest::port
cortex_neuron::handles_test_event( nest::SpikeEvent&, nest::rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw nest::UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline void
cortex_neuron::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  nest::Archiving_Node::get_status( d );
}

inline void
cortex_neuron::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d );         // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  nest::Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
}

} // namespace

#endif // CORTEX_NEURON_H
