/*
 *  radial_basis_function_input.cpp
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


#include "radial_basis_function_input.h"

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

mynest::radial_basis_function_input::Parameters_::Parameters_()
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
mynest::radial_basis_function_input::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, nest::names::rate, Rate_ );
  def< double >( d, nest::names::noise,Noise_ );
  def< double >( d, nest::names::mean, Mean_ );
  def< double >( d, mynames::sdev, SDev_ );
  def< std::string >( d, nest::names::filenames, FileDesired_ );
}

void
mynest::radial_basis_function_input::Parameters_::set( const DictionaryDatum& d )
{
  // allow setting the neuron parameters
  updateValue< double >( d, nest::names::rate, Rate_ );
  updateValue< double >( d, nest::names::noise, Noise_ );
  updateValue< double >( d, nest::names::mean, Mean_ );
  updateValue< double >( d, mynames::sdev, SDev_ );
  updateValue< std::string >( d, nest::names::filenames, FileDesired_ );
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::radial_basis_function_input::radial_basis_function_input()
  : Archiving_Node()
  , P_()
{
}

mynest::radial_basis_function_input::radial_basis_function_input( const radial_basis_function_input& n )
  : Archiving_Node( n )
  , P_( n.P_ )
{
}

mynest::radial_basis_function_input::~radial_basis_function_input()
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
mynest::radial_basis_function_input::init_state_( const Node& proto )
{
}


void
mynest::radial_basis_function_input::init_buffers_()
{
  Archiving_Node::clear_history();
}

void
mynest::radial_basis_function_input::calibrate()
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
    while ( not( DesFile.eof() ) )
    {
      DesFile >> Val;
      V_.DesValues_.push_back( Val );
    }
    
  rng_ = nest::kernel().rng_manager.get_rng( get_thread() );
}


void
mynest::radial_basis_function_input::update( nest::Time const& origin, const long from, const long to )
{
  assert(
    to >= 0 && ( nest::delay ) from < nest::kernel().connection_manager.get_min_delay() );
  assert( from < to );

  for ( long lag = from; lag < to; ++lag )
  {
    int t = origin.get_steps() + lag;
    double Desired = V_.DesValues_[ t ];
    double Value = 0.0; // ranges from 0 (0 Hz) to 1 (Max Frequency)
    double MaxFraction = nest::Time::get_resolution().get_ms() * P_.Rate_ / 1000.0; // Since it is updated each X ms
    double MaxFractionNoise = nest::Time::get_resolution().get_ms() * P_.Noise_ / 1000.0; // Since it is updated each X ms

    Value = exp(-pow(((Desired-P_.Mean_)/P_.SDev_), 2 ));

    if ( MaxFraction * Value > rng_->drand() || MaxFractionNoise > rng_->drand() )
      {
        nest::SpikeEvent se;
        se.set_multiplicity( 1 );
        nest::kernel().event_delivery_manager.send( *this, se, lag );
        set_spiketime( nest::Time::step( t + 1 ) );
      }
  }
}

void
mynest::radial_basis_function_input::handle( nest::SpikeEvent& e )
{

}
