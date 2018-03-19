/*
 *  stdp_connection_sinexp.cpp
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
 
 /*

   Alberto Antonietti
   alberto.antonietti@polimi.it
  
   Cerebellar PF-PC Plasticity with an exp. sin. Kernel LTP and LTD
 
 */

#include "stdp_connection_sinexp.h"

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connector_model.h"
#include "event.h"
#include "kernel_manager.h"

// Includes from sli:
#include "dictdatum.h"

namespace mynest
{
//
// Implementation of class STDPSinExpCommonProperties.
//

STDPSinExpCommonProperties::STDPSinExpCommonProperties()
  : nest::CommonSynapseProperties()
  , A_plus_( 1.0 )
  , A_minus_( 1.5 )
  , Wmin_( 0.0 )
  , Wmax_( 200.0 )
  , vtC_ ( 0 )
{
}

void STDPSinExpCommonProperties::get_status( DictionaryDatum& d ) const{
  nest::CommonSynapseProperties::get_status( d );

  def< double >( d, "A_plus", A_plus_ );
  def< double >( d, "A_minus", A_minus_ );
  def< double >( d, "Wmin", Wmin_ );
  def< double >( d, "Wmax", Wmax_ );
}

void STDPSinExpCommonProperties::set_status( const DictionaryDatum& d, nest::ConnectorModel& cm ){
  nest::CommonSynapseProperties::set_status( d, cm );

  updateValue< double >( d, "A_plus", A_plus_ );
  updateValue< double >( d, "A_minus", A_minus_ );
  updateValue< double >( d, "Wmin", Wmin_ );
  updateValue< double >( d, "Wmax", Wmax_ );
}

} // End of namespace mynest
