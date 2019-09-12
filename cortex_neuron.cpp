#include "cortex_neuron.h"

// C++ includes:
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <fstream>

#include "poisson_generator.h"

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

mynest::cortex_neuron::Parameters_::Parameters_()
  : trial_length_( 1000 )
  , joint_id_( 0 )
  , fiber_id_( 0 )
  , fibers_per_joint_( 100 )
  , rbf_sdev_( 10.0 )
  , baseline_rate_( 10.0 )
  , gain_rate_( 10.0 )
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
mynest::cortex_neuron::Parameters_::get( DictionaryDatum& d ) const
{
  def< long >( d, mynames::trial_length, trial_length_ );
  def< long >( d, mynames::joint_id, joint_id_ );
  def< long >( d, mynames::fiber_id, fiber_id_ );
  def< long >( d, mynames::fibers_per_joint, fibers_per_joint_ );
  def< double >( d, mynames::rbf_sdev, rbf_sdev_ );
  def< double >( d, mynames::baseline_rate, baseline_rate_ );
  def< double >( d, mynames::gain_rate, gain_rate_ );
}

void
mynest::cortex_neuron::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< long >( d, mynames::trial_length, trial_length_ );
  if ( trial_length_ <= 0 )
  {
    throw nest::BadProperty( "The trial length cannot be zero or negative." );
  }
  updateValue< long >( d, mynames::joint_id, joint_id_ );
  updateValue< long >( d, mynames::fiber_id, fiber_id_ );
  updateValue< long >( d, mynames::fibers_per_joint, fibers_per_joint_ );
  updateValue< double >( d, mynames::rbf_sdev, rbf_sdev_ );
  updateValue< double >( d, mynames::baseline_rate, baseline_rate_ );
  updateValue< double >( d, mynames::gain_rate, gain_rate_ );
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::cortex_neuron::cortex_neuron()
  : Archiving_Node()
  , P_()
  , V_()
{
}

mynest::cortex_neuron::cortex_neuron( const cortex_neuron& n )
  : Archiving_Node( n )
  , P_( n.P_ )
{
}

mynest::cortex_neuron::~cortex_neuron()
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
mynest::cortex_neuron::init_state_( const Node& proto )
{
}


void
mynest::cortex_neuron::init_buffers_()
{
  B_.spike_mult.clear();
  B_.spike_lag.clear();
  Archiving_Node::clear_history();
}

void
mynest::cortex_neuron::calibrate()
{
  double time_res = nest::Time::get_resolution().get_ms();  // 0.1
  long from = 0;
  long to = (double)P_.trial_length_ / time_res;

  librandom::RngPtr rng = nest::kernel().rng_manager.get_rng( get_thread() );

  for (long lag = from; lag < to; ++lag )
  {
    double t = ( from + lag ) * time_res * 1e-3;  // [s]
    double sdev = P_.rbf_sdev_;
    double mean = P_.fiber_id_;
    double desired = P_.fibers_per_joint_ * ( 0.5 + 0.5*std::sin( 2*M_PI * t ) );
    double rate = P_.baseline_rate_ * exp(-pow(((desired - mean) / sdev), 2 ));

    V_.poisson_dev_.set_lambda( time_res * rate * 1e-3 );

    long n_spikes = V_.poisson_dev_.ldev( rng );

    if ( n_spikes > 0 ) // we must not send events with multiplicity 0
    {
      B_.spike_mult.push_back( n_spikes );
      B_.spike_lag.push_back( lag );
    }
  }
}


void
mynest::cortex_neuron::update( nest::Time const& origin, const long from, const long to )
{
  assert( to >= 0 );
  assert( static_cast<nest::delay>(from) < nest::kernel().connection_manager.get_min_delay() );
  assert( from < to );

  double time_res = nest::Time::get_resolution().get_ms();  // 0.1
  long trial_len_ticks = (double)P_.trial_length_ / time_res;

  long const T = origin.get_steps() + from;

  if ( T % trial_len_ticks == 0 )
  {
    for ( size_t i = 0; i < B_.spike_lag.size(); i++ )
    {
      nest::SpikeEvent se;
      int n_spikes = B_.spike_mult[i];
      long lag = B_.spike_lag[i];  // it works fine without adding any offset

      se.set_multiplicity( n_spikes );
      nest::kernel().event_delivery_manager.send( *this, se, lag );
    }
  }
}

void
mynest::cortex_neuron::handle( nest::SpikeEvent& e )
{
}
