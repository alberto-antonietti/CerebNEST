/*
 *  radial_basis_function_input.h
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
Name: radial_basis_function_input - ...

Description:

...


Receives: Nothing

Sends: SpikeEvent

Parameters:
...

Author: Alberto Antonietti
FirstVersion: October 2017
*/


#ifndef RADIAL_BASIS_FUNCTION_INPUT_H
#define RADIAL_BASIS_FUNCTION_INPUT_H

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
class radial_basis_function_input : public nest::Archiving_Node
{

public:
  radial_basis_function_input();
  radial_basis_function_input( const radial_basis_function_input& );
  ~radial_basis_function_input();

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
    double Rate_;
    double Noise_;
    double Mean_;
    double SDev_;
    std::string FileDesired_; //!< indicates the Path of the File with the
                              //!< Input Variable

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dicitonary
  };

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    std::vector< double >
      DesValues_; //!< stores the Desired Values of Trajectory (VOR Protocol)
  };
  
  Parameters_ P_;
  Variables_ V_;
};

inline nest::port
radial_basis_function_input::send_test_event( nest::Node& target,
  nest::rport receptor_type,
  nest::synindex,
  bool )
{
  nest::SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline nest::port
radial_basis_function_input::handles_test_event( nest::SpikeEvent&, nest::rport receptor_type )
{
  return 0;
}

inline void
radial_basis_function_input::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  nest::Archiving_Node::get_status( d );

}

inline void
radial_basis_function_input::set_status( const DictionaryDatum& d )
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

#endif // RADIAL_BASIS_FUNCTION_INPUT_H
