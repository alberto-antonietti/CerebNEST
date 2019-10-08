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
  , to_file_( false )
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
  def< bool >( d, nest::names::to_file, to_file_ );
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
  updateValue< bool >( d, nest::names::to_file, to_file_ );
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
  Archiving_Node::clear_history();

  double time_res = nest::Time::get_resolution().get_ms();  // 0.1
  long ticks = 50.0 / time_res;  // 50ms
  // long ticks = 100.0 / time_res;  // 100ms
  for ( long i = 0; i < ticks; i++ )
  {
    B_.in_spikes_.push_back(0);
  }

  B_.traj_.resize(4, std::vector<double>(P_.trial_length_));
  std::ifstream traj_file("JointTorques.dat");

  V_.joint_scale_factors_.resize(4);
  traj_file >> V_.joint_scale_factors_[0]
            >> V_.joint_scale_factors_[1]
            >> V_.joint_scale_factors_[2]
            >> V_.joint_scale_factors_[3];

  for (int i = 0; i < P_.trial_length_; i++)
  {
    traj_file >> B_.traj_[0][i]
              >> B_.traj_[1][i]
              >> B_.traj_[2][i]
              >> B_.traj_[3][i];
  }

  for (int i = 0; i < P_.trial_length_; i++)
  for (int j = 0; j < 4; j++)
  {
    B_.traj_[j][i] = 0.5 + 0.5 * B_.traj_[j][i];
  }
}

void
mynest::cortex_neuron::calibrate()
{
  if ( P_.to_file_ )
  {
    V_.out_file_.open( "cortex_out.dat" );
  }
}


void
mynest::cortex_neuron::update( nest::Time const& origin, const long from, const long to )
{
  assert( to >= 0 );
  assert( static_cast<nest::delay>(from) < nest::kernel().connection_manager.get_min_delay() );
  assert( from < to );

  double time_res = nest::Time::get_resolution().get_ms();  // 0.1
  long trial_ticks = (double)P_.trial_length_ / time_res;
  librandom::RngPtr rng = nest::kernel().rng_manager.get_rng( get_thread() );

  for ( long lag = from; lag < to; ++lag )
  {
    long tick = origin.get_steps() + lag;
    // int n_spikes = B_.out_spikes_[tick % trial_ticks];
    double t = tick * time_res * 1e-3;  // [s]
    double sdev = P_.rbf_sdev_;
    double mean = P_.fiber_id_;
    // double desired = P_.fibers_per_joint_ * ( 0.5 + 0.5*std::sin( 2*M_PI * t ) );
    double desired = P_.fibers_per_joint_ * B_.traj_[P_.joint_id_][(int)(tick * time_res) % P_.trial_length_];

    double baseline_rate;
    int j_id = P_.joint_id_;

    baseline_rate = P_.baseline_rate_;

    if ( j_id == 1 )  // Second joint
    {
      baseline_rate = std::max( 0.0, V_.in_rate_ );
    }

    double rate = baseline_rate * exp(-pow(((desired - mean) / sdev), 2 ));

    V_.poisson_dev_.set_lambda( time_res * rate * 1e-3 );

    long n_spikes = V_.poisson_dev_.ldev( rng );

    if ( n_spikes > 0 )
    {
      nest::SpikeEvent se;
      se.set_multiplicity( n_spikes );
      nest::kernel().event_delivery_manager.send( *this, se, lag );

      // set the spike times, respecting the multiplicity
      for ( unsigned long i = 0; i < n_spikes; i++ )
      {
        set_spiketime( nest::Time::step( tick ) );
      }
    }

    long buf_size = B_.in_spikes_.size();
    B_.in_spikes_[ tick % buf_size ] = 0;
  }
}

void
mynest::cortex_neuron::handle( nest::SpikeEvent& e )
{
  nest::Time stamp = e.get_stamp();
  long t = stamp.get_steps();

  long buf_size = B_.in_spikes_.size();

  B_.in_spikes_[ t % buf_size ] += e.get_weight() * e.get_multiplicity();

  int spike_count = 0;
  for ( long i = 0; i < buf_size; i++ )
  {
    spike_count += B_.in_spikes_[ i ];
  }
  double time_res = nest::Time::get_resolution().get_ms();  // 0.1
  V_.in_rate_ = std::max( 0.0, 1000.0 * spike_count / (buf_size * time_res) );

  if ( P_.to_file_ )
  {
    V_.out_file_ << "in_rate\t" << V_.in_rate_ << std::endl;
  }
  // std::cout << "Rate: " << V_.in_rate_ << std::endl;
}
