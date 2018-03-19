/*
 *  closed_loop_neuron.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


/* BeginDocumentation
Name: closed_loop_neuron - Neuron that repeats incoming spikes.

Description:

The parrot neuron simply emits one spike for every incoming spike.
An important application is to provide identical poisson spike
trains to a group of neurons. The poisson_generator sends a different
spike train to each of its target neurons. By connecting one
poisson_generator to a closed_loop_neuron and then that closed_loop_neuron to
a group of neurons, all target neurons will receive the same poisson
spike train.

Remarks:

- Weights on connection to the closed_loop_neuron are ignored.
- Weights on connections from the closed_loop_neuron are handled as usual.
- Delays are honored on incoming and outgoing connections.

Only spikes arriving on connections to port 0 will be repeated.
Connections onto port 1 will be accepted, but spikes incoming
through port 1 will be ignored. This allows setting exact pre-
and post-synaptic spike times for STDP protocols by connecting
two parrot neurons spiking at desired times by, e.g., a
stdp_synapse onto port 1 on the post-synaptic parrot neuron.

Receives: SpikeEvent

Sends: SpikeEvent

Parameters:
No parameters to be set in the status dictionary.

Author: David Reichert, Abigail Morrison, Alexander Seeholzer, Hans Ekkehard
Plesser
FirstVersion: May 2006
*/


/**
 * The parrot neuron emits one spike for every incoming spike,
 * but may use multiplicity to indicate number of spikes in a single
 * time step.
 * Instead of the accumulated weigths of the incoming spikes, the
 * number of the spikes is stored within a ring buffer.
 */

#ifndef CLOSED_LOOP_NEURON_H
#define CLOSED_LOOP_NEURON_H

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
class closed_loop_neuron : public nest::Archiving_Node
{

public:
  closed_loop_neuron();
  closed_loop_neuron( const closed_loop_neuron& );
  ~closed_loop_neuron();

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

  /**
     Buffers and accumulates the number of incoming spikes per time step;
     RingBuffer stores doubles; for now the numbers are casted.
  */
  struct Buffers_
  {
    std::vector< double > spike_gids_;
  };


  //! Model parameters
  struct Parameters_
  {
    double
      Gain_; //!< Gain to multiply the output for VOR - CR_Threshold for EBCC
    double NumDCN_; //!< Total Number of DCN connected to the Closed_loop_neuron
    double FirstDCN_;      //!< The number (GID) of the First DCN
    bool Positive_;        //!< = True if the output goes to positive IO
    bool ToFile_;          //!< = True if the neuron writes the output files
    double Protocol_;      //!< 1.0 EBCC, 2.0 VOR
    double USOnset_;       //!< in ms the relative onset of US (EBCC Protocol)
    double USDuration_;    //!< in ms the duration of the US (e.g. 100 ms) (EBCC
                           //Protocol)
                           //!< is the total number of trials (VOR Protocol)
    double TrialDuration_; //!< in ms the duration of each trial
    double Phase_; //!< indicates the number of trial when the Extinction begins
                   //(EBCC and VOR Protocols)
    std::string FileDesired_; //!< indicates the Path of the File with the
                              //Desired Trajectory (VOR Protocol)

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dicitonary
  };

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    double DCNAvg_; //!< the actual value of the DCNAvg
    std::vector< double >
      DCNBuffer_; //!< stores the DCNBuffer for the mobile window average
    std::vector< double >
      DesValues_; //!< stores the Desired Values of Trajectory (VOR Protocol)
    double
      OutputVariables_[ 2 ];   //!< actual Positive and Negative DCN Firing Rate
    std::ofstream OutputFile_; //!< OutputFile
    std::ofstream CRFile_;     //!< CRFIle
    bool CRFlag_;              //!< becomes true when a CR is detected
    int Trial_;                //!< counts the number of trials
  };
  Buffers_ B_;
  Parameters_ P_;
  Variables_ V_;
};

inline nest::port
closed_loop_neuron::send_test_event( nest::Node& target,
  nest::rport receptor_type,
  nest::synindex,
  bool )
{
  nest::SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline nest::port
closed_loop_neuron::handles_test_event( nest::SpikeEvent&, nest::rport receptor_type )
{
  // Allow connections only to port 0
  if ( receptor_type == 0 )
    return receptor_type;
  else
    throw nest::UnknownReceptorType( receptor_type, get_name() );
}

inline void
closed_loop_neuron::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  nest::Archiving_Node::get_status( d );

  //( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
closed_loop_neuron::set_status( const DictionaryDatum& d )
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

#endif // CLOSED_LOOP_NEURON_H
