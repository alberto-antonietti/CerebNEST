import nest
import sys

nest.Install("albertomodule")

VT = nest.Create("volume_transmitter_alberto", 1)

CLOSED = nest.Create("closed_loop_neuron", 1)
RBF = nest.Create("radial_basis_function_input", 1)
PRE = nest.Create("iaf_cond_exp", 1)
POST = nest.Create("iaf_cond_exp", 1)

nest.SetDefaults('stdp_synapse_sinexp',{"A_minus":  -0.10,   # double - Amplitude of weight change for depression
                                        "A_plus":    0.01,   # double - Amplitude of weight change for facilitation 
                                        "Wmin":      0.0,    # double - Minimal synaptic weight 
                                        "Wmax":      4.0,    # double - Maximal synaptic weight
                                        "vt":        VT[0] })

nest.Connect(PRE,POST,{'rule': 'one_to_one'},'stdp_synapse_sinexp')
A=nest.GetConnections(PRE,POST)
nest.SetStatus(A,{'vt_num': 0})

nest.Connect(POST,PRE,{'rule': 'one_to_one'},'stdp_synapse_cosexp')
A=nest.GetConnections(POST,PRE)
nest.SetStatus(A,{'vt_num': 0})

sys.exit(0) #Everything went fine
