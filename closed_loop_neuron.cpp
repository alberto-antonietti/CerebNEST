/*
 *  closed_loop_neuron.cpp
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

#include "closed_loop_neuron.h"

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

mynest::closed_loop_neuron::Parameters_::Parameters_()
  : Gain_( 1.0 )      // adimensional
  , NumDCN_( 1.0 )    // Number of DCN as output
  , Positive_( true ) // Positive or Negative Neuron
  , ToFile_( false )  // If it writes the OutputFile.dat
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
mynest::closed_loop_neuron::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, mynames::gain, Gain_ );
  def< double >( d, mynames::num_dcn, NumDCN_ );
  def< double >( d, mynames::first_dcn, FirstDCN_ );
  def< bool >( d, mynames::positive, Positive_ );
  def< bool >( d, nest::names::to_file, ToFile_ );
  def< double >( d, mynames::protocol, Protocol_ );
  def< double >( d, nest::names::Tstart, USOnset_ );
  def< double >( d, nest::names::Tstop, USDuration_ );
  def< double >( d, mynames::Tduration, TrialDuration_ );
  def< double >( d, nest::names::phase, Phase_ );
  def< std::string >( d, nest::names::filenames, FileDesired_ );
}

void
mynest::closed_loop_neuron::Parameters_::set( const DictionaryDatum& d )
{
  // allow setting the neuron parameters
  updateValue< double >( d, mynames::gain, Gain_ );
  updateValue< double >( d, mynames::num_dcn, NumDCN_ );
  updateValue< double >( d, mynames::first_dcn, FirstDCN_ );
  updateValue< bool >( d, mynames::positive, Positive_ );
  updateValue< bool >( d, nest::names::to_file, ToFile_ );
  updateValue< double >( d, mynames::protocol, Protocol_ );
  updateValue< double >( d, nest::names::Tstart, USOnset_ );
  updateValue< double >( d, nest::names::Tstop, USDuration_ );
  updateValue< double >( d, mynames::Tduration, TrialDuration_ );
  updateValue< double >( d, nest::names::phase, Phase_ );
  updateValue< std::string >( d, nest::names::filenames, FileDesired_ );
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::closed_loop_neuron::closed_loop_neuron()
  : Archiving_Node()
  , P_()
{
}

mynest::closed_loop_neuron::closed_loop_neuron( const closed_loop_neuron& n )
  : Archiving_Node( n )
  , P_( n.P_ )
{
}

mynest::closed_loop_neuron::~closed_loop_neuron()
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
mynest::closed_loop_neuron::init_state_( const Node& proto )
{
}


void
mynest::closed_loop_neuron::init_buffers_()
{
  B_.spike_gids_.clear();
  Archiving_Node::clear_history();
}

void
mynest::closed_loop_neuron::calibrate()
{
  V_.DCNAvg_ = 0.0;
  V_.OutputVariables_[ 0 ] = 0.0; // POSITIVE VARIABLE
  V_.OutputVariables_[ 1 ] = 0.0; // NEGATIVE VARIABLE
  for ( int i = 0; i < 50; i++ )
    V_.DCNBuffer_.push_back( 0.0 );

  if ( P_.Protocol_ == 2.0 )
  {
    std::ifstream DesFile;
    double Val = 0.0;
    DesFile.open( ( P_.FileDesired_ ).c_str() );
    if ( DesFile.is_open() == false )
    {
      std::cout << "ERROR! Could not open the Desired Variable File"
                << std::endl;
      return;
    }
    for ( int i = 0; i < P_.TrialDuration_ * P_.USDuration_; i++ )
    {
      DesFile >> Val;
      V_.DesValues_.push_back( Val );
    }
  }
  V_.CRFlag_ = false;
  if ( P_.ToFile_ )
  {
    V_.OutputFile_.open( "OutputFile.dat" );
    if ( P_.Protocol_ == 1.0 )
      V_.CRFile_.open( "CR.dat" );
  }
  V_.Trial_ = 0;
  rng_ = nest::kernel().rng_manager.get_rng( get_thread() );
}


void
mynest::closed_loop_neuron::update( nest::Time const& origin, const long from, const long to )
{
  assert(
    to >= 0 && ( nest::delay ) from < nest::kernel().connection_manager.get_min_delay() );
  assert( from < to );

  double tau_time_constant = 0.060;
  double kernel_amplitude = sqrt( 2.0 / tau_time_constant );
  V_.OutputVariables_[ 0 ] *= exp( -0.001 / tau_time_constant );
  V_.OutputVariables_[ 1 ] *= exp( -0.001 / tau_time_constant );


  for ( long lag = from; lag < to; ++lag )
  {
    const unsigned long current_spikes_n =
      static_cast< unsigned long >( B_.spike_gids_.size() );
    if ( current_spikes_n > 0 )
    {
      for ( unsigned long i = 0; i < current_spikes_n; i++ )
      {
        if ( B_.spike_gids_[ i ] < P_.FirstDCN_ + ( P_.NumDCN_ ) / 2 )
        { // POSITIVE DCN
          V_.OutputVariables_[ 0 ] += kernel_amplitude;
        }
        else
        { // NEGATIVE DCN
          V_.OutputVariables_[ 1 ] += kernel_amplitude;
        }
      }
    }

    V_.DCNBuffer_.erase( V_.DCNBuffer_.begin() );
    V_.DCNBuffer_.push_back(
      V_.OutputVariables_[ 0 ] - V_.OutputVariables_[ 1 ] );
    V_.DCNAvg_ = accumulate( V_.DCNBuffer_.begin(), V_.DCNBuffer_.end(), 0.0 )
      / ( V_.DCNBuffer_.size() * ( ( P_.NumDCN_ ) / 2 ) );
    int t = origin.get_steps() + lag;
    double Error = 0.0;
    if ( t % ( int ) P_.TrialDuration_ == 0 )
    {
      V_.Trial_++;
      if ( P_.ToFile_ )
        std::cout << "Trial : " << V_.Trial_ << std::endl;
      V_.CRFlag_ = false;
    }

    //!< EBCC PROTOCOL - ACQUISITION AND EXTINCTION
    if ( P_.Protocol_ == 1.0 && V_.Trial_ <= P_.Phase_ )
    { // EBCC ACQUISITION
      if ( ( t - ( V_.Trial_ - 1 ) * P_.TrialDuration_ ) >= P_.USOnset_
        && ( t - ( V_.Trial_ - 1 ) * P_.TrialDuration_ ) < P_.USOnset_
            + P_.USDuration_
        && P_.Positive_ )
      {
        if ( !V_.CRFlag_ )
          Error = 1.0; // NO CR before
        else
          Error = 0.5; // CR before
      }
    }
    else if ( P_.Protocol_ == 1.0 && V_.Trial_ > P_.Phase_ )
    {
      // EBCC EXTINCTION
      Error = 0.0;
    }
    if ( P_.Protocol_ == 1.0 )
    {
      // EBCC ACQUISITION AND EXTINCTION
      // CR DETECTION
      if ( ( t - ( V_.Trial_ - 1 ) * P_.TrialDuration_ ) < P_.USOnset_
        && ( t - ( V_.Trial_ - 1 ) * P_.TrialDuration_ )
          > ( P_.USOnset_ - 200.0 )
        && V_.DCNAvg_ > P_.Gain_
        && !V_.CRFlag_ )
      {
        V_.CRFlag_ = true;
        if ( P_.ToFile_ )
        {
          V_.CRFile_ << t << std::endl;
        }
      }
      if ( P_.ToFile_ )
        V_.OutputFile_ << V_.DCNAvg_ << std::endl;
    }

    //!< VOR PROTOCOL - ACQUISITION AND EXTINCTION
    if ( P_.Protocol_ == 2.0 )
    { // VOR ACQUISITION
      double Desired = V_.DesValues_[ t ];
      Error = Desired - V_.DCNAvg_ * P_.Gain_;
      if ( P_.ToFile_ )
        V_.OutputFile_ << Desired << "\t" << V_.DCNAvg_* P_.Gain_ << "\t"
                       << Error << std::endl;
    }


    //!< Check that the Error variable ranges in [-1, 1]
    if ( Error > 1.0 )
    {
      Error = 1.0;
      std::cout << "Warning: Error Exceeds the nominal range (> 1.0)"
                << std::endl;
    }
    if ( Error < -1.0 )
    {
      Error = -1.0;
      std::cout << "Warning: Error Exceeds the nominal range (< -1.0)"
                << std::endl;
    }


    if ( P_.Positive_ && Error > 0.0 )
    { // Positive Neuron and Positive Error
      if ( 0.01 * Error > rng_->drand() )
      {
        // When the error is max == 1, we have a mean firing rate equal to 10 Hz
        // (max freq)
        nest::SpikeEvent se;
        se.set_multiplicity( 1 );
        nest::kernel().event_delivery_manager.send( *this, se, lag );
        set_spiketime( nest::Time::step( t + 1 ) );
      }
    }
    else if ( !P_.Positive_ && Error < 0.0 )
    {
      if ( -0.01 * Error > rng_->drand() )
      {
        // When the error is max == -1, we have a mean firing rate equal to 10
        // Hz (max freq)
        nest::SpikeEvent se;
        se.set_multiplicity( 1 );
        nest::kernel().event_delivery_manager.send( *this, se, lag );
        set_spiketime( nest::Time::step( t + 1 ) );
      }
    }
  }
  B_.spike_gids_.clear();
}

void
mynest::closed_loop_neuron::handle( nest::SpikeEvent& e )
{
  // Repeat only spikes incoming on port 0, port 1 will be ignored
  if ( 0 == e.get_rport() )
  {
    B_.spike_gids_.push_back( e.get_sender_gid() );
  }
}
