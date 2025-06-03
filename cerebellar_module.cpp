
/*
*  cerebellar_module.cpp
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
*  2024-11-05 14:48:37.050316
*/

// Include from NEST
#include "nest_extension_interface.h"

// include headers with your own stuff
#include "eglif_io_nestml.h"
#include "eglif_pc_nestml.h"
#include "eglif_mli.h"
#include "eglif_dcnp.h"
#include "stdp_synapse_sinexp.h"
#include "stdp_synapse_alpha.h"

#include "tracking_neuron_nestml.h"
#include  "basic_neuron_nestml.h"
#include  "diff_neuron_nestml.h"
#include  "rb_neuron_nestml.h"
#include  "tracking_neuron_planner_nestml.h"
#include "state_neuron_nestml.h"

#include "eglif_cond_alpha_multisyn.h"


class cerebellar_module : public nest::NESTExtensionInterface
{
  public:
    cerebellar_module() {}
    ~cerebellar_module() {}

    void initialize() override;
};

// Define module instance outside of namespace to avoid name-mangling problems
cerebellar_module cerebellar_module_LTX_module;

void cerebellar_module::initialize()
{
    // register neurons, synapse or device models
    register_eglif_cond_alpha_multisyn("eglif_cond_alpha_multisyn");
    nest::register_eglif_io_nestml("eglif_io_nestml");
    nest::register_eglif_pc_nestml("eglif_pc_nestml");
    nest::register_eglif_mli("eglif_mli");
    nest::register_eglif_dcnp("eglif_dcnp");
    
    register_tracking_neuron_nestml("tracking_neuron_nestml");
    register_basic_neuron_nestml("basic_neuron_nestml");
    register_diff_neuron_nestml("diff_neuron_nestml");
    register_rb_neuron_nestml("rb_neuron_nestml");
    register_tracking_neuron_planner_nestml("tracking_neuron_planner_nestml");
    register_state_neuron_nestml("state_neuron_nestml");


    nest::register_stdp_synapse_sinexp( "stdp_synapse_sinexp" );
    nest::register_stdp_synapse_alpha( "stdp_synapse_alpha" );
}
