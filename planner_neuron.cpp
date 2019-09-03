#include "planner_neuron.h"

// C++ includes:
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <fstream>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "exceptions.h"
#include "kernel_manager.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"


/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

mynest::planner_neuron::Parameters_::Parameters_()
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
mynest::planner_neuron::Parameters_::get( DictionaryDatum& d ) const
{
}

void
mynest::planner_neuron::Parameters_::set( const DictionaryDatum& d )
{
  // allow setting the neuron parameters
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::planner_neuron::planner_neuron()
  : Archiving_Node()
  , P_()
  , V_()
{
}

mynest::planner_neuron::planner_neuron( const planner_neuron& n )
  : Archiving_Node( n )
  , P_( n.P_ )
{
}

mynest::planner_neuron::~planner_neuron()
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
mynest::planner_neuron::init_state_( const Node& proto )
{
}


void
mynest::planner_neuron::init_buffers_()
{
  Archiving_Node::clear_history();
}

void
mynest::planner_neuron::calibrate()
{
  rng_ = nest::kernel().rng_manager.get_rng( get_thread() );
}


void
mynest::planner_neuron::update( nest::Time const& origin, const long from, const long to )
{
  assert(
    to >= 0 && ( nest::delay ) from < nest::kernel().connection_manager.get_min_delay() );
  assert( from < to );
}

void
mynest::planner_neuron::handle( nest::SpikeEvent& e )
{
}
