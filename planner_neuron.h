#ifndef PLANNER_NEURON_H
#define PLANNER_NEURON_H

#include <fstream>

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "mynames.h"

namespace mynest
{
class planner_neuron : public nest::Archiving_Node
{

public:
  planner_neuron();
  planner_neuron( const planner_neuron& );
  ~planner_neuron();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding,
   * Overloading, and Hiding
   */
  using nest::Node::handle;
  using nest::Node::handles_test_event;

  nest::port send_test_event( nest::Node&, nest::rport, nest::synindex, bool );
  void handle( nest::SpikeEvent& );
  nest::port handles_test_event( nest::SpikeEvent&, nest::rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  void init_state_( const nest::Node& proto );
  void init_buffers_();
  void calibrate();
  void update( nest::Time const&, const long, const long );

  librandom::RngPtr rng_;

  //! Model parameters
  struct Parameters_
  {
    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dicitonary
  };

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
  };
  Parameters_ P_;
  Variables_ V_;
};

inline nest::port
planner_neuron::send_test_event( nest::Node& target,
  nest::rport receptor_type,
  nest::synindex,
  bool )
{
  nest::SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline nest::port
planner_neuron::handles_test_event( nest::SpikeEvent&, nest::rport receptor_type )
{
  // Allow connections only to port 0
  if ( receptor_type == 0 )
  {
    return receptor_type;
  }
  else
    throw nest::UnknownReceptorType( receptor_type, get_name() );
}

inline void
planner_neuron::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  nest::Archiving_Node::get_status( d );

  //( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
planner_neuron::set_status( const DictionaryDatum& d )
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

#endif // PLANNER_NEURON_H
