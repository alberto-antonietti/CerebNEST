import nest
import numpy as np
import json
from validation_synapse_alpha import plot_synaptic_weight, reformat_weights, plot_synaptic_matrix, correct_spike_times
import sys



f = open('eglif_params.json')
params = json.load(f)
f.close()

IO_params = params["EGLIF_model"]["IO"]
GR_params = params["EGLIF_model"]["Granule"]
MLI_params = params["EGLIF_model"]["Stellate"]

N = 3
M = 5


Aminus_values = [-0.005, -0.003, -0.007]
Aplus = 0.0

weights_dict = {}

for idx, value in enumerate(Aminus_values):
    TEST_NUMBER = 7 + idx/10
    nest.ResetKernel()
    nest.Install("cerebellar_module")
    nest.SetKernelStatus({"resolution": 0.1})

    # Create the network and tune parameters
    IO = nest.Create("eglif_io_nestml", N) # IO
    nest.SetStatus(IO, IO_params)
    MLI = nest.Create("eglif_mli", N) # PC
    nest.SetStatus(MLI, MLI_params)
    Gr = nest.Create("eglif_cond_alpha_multisyn", M) # Gr


    # Connections
    nest.Connect(IO, MLI, "one_to_one", syn_spec={"receptor_type": 3})


    wr = nest.Create("weight_recorder")
    nest.CopyModel("stdp_synapse_alpha", "my_stdp_synapse_rec", {"weight_recorder": wr, "Aplus": Aplus, "Aminus": value})
    
    conn_spec_dict = {'rule': 'fixed_indegree', 'indegree': 4}
    syn_spec_dict = {'synapse_model': "my_stdp_synapse_rec", "receptor_type": 1}

    nest.Connect(Gr, MLI, conn_spec = conn_spec_dict, syn_spec = syn_spec_dict)

    pf_pc_conns = nest.GetConnections(synapse_model="my_stdp_synapse_rec")
    IO_mli_conns = nest.GetConnections(source = IO, target = MLI)

    # Devices   
    cf_recorder = nest.Create("spike_recorder", {"record_to": "memory"})
    MLI_recorder = nest.Create("spike_recorder",  {"record_to": "memory"})
    Gr_recorder = nest.Create("spike_recorder",  {"record_to": "memory"})

    nest.Connect(IO, cf_recorder)
    nest.Connect(MLI, MLI_recorder)
    nest.Connect(Gr, Gr_recorder)

    # SIMULATION
    weights = []
    for i in range(1000):
        nest.Simulate(1.0)
        weights.append(nest.GetStatus(pf_pc_conns, "weight"))

    weight_matrix = reformat_weights(weights)

    weights_dict[idx] = weight_matrix
    corrected_cf_tms, cf_evs, corrected_MLI_tms, MLI_evs, MLI_complex_tms, MLI_complex_evs, corrected_gr_spikes, Gr_evs = correct_spike_times(cf_recorder, MLI_recorder, Gr_recorder)
    
    data = nest.GetStatus(wr)[0]

    plot_synaptic_weight(pf_pc_conns, data, IO, MLI, Gr, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_MLI_tms, MLI_evs, TEST_NUMBER) 

# Check across the trials if, when LTD happens, the ratio of weights is the same of Aminus

weights_1 = weights_dict.get(0)
weights_2 = weights_dict.get(1)
weights_3 = weights_dict.get(2)

LTD_1 = []
LTD_2 = []
LTD_3 = []

for n, conn in enumerate(pf_pc_conns):
    # Retrieve weights
    wgh_1 = weights_1[n, :]
    wgh_2 = weights_2[n, :]
    wgh_3 = weights_3[n, :]

    # Get the indexes of where LTD happens
    diff_wgh_1 = np.diff(wgh_1)
    diff_wgh_2 = np.diff(wgh_2)
    diff_wgh_3 = np.diff(wgh_3)

    LTD_index_1 = np.where(diff_wgh_1 < 0) # Non sto considerando la situazione dove il peso è già 0 e rimane 0
    LTD_index_2 = np.where(diff_wgh_2 < 0)
    LTD_index_3 = np.where(diff_wgh_3 <  0)


    # First assert LTD happens at the same times
    assert(np.array_equal(LTD_index_1, LTD_index_2))
    assert(np.array_equal(LTD_index_2, LTD_index_3))
    
    LTD_amounts_1 = diff_wgh_1[LTD_index_1]
    LTD_amounts_2 = diff_wgh_2[LTD_index_2]
    LTD_amounts_3 = diff_wgh_3[LTD_index_3]

    # Check that the ratio corresponds to the ratio between the Aminus values
    ratio1 = np.round(np.average(LTD_amounts_2) / np.average(LTD_amounts_1), 1)
    ratio2 = np.round(np.average(LTD_amounts_3 )/ np.average(LTD_amounts_1), 1)
    ratio3 = np.round(np.average(LTD_amounts_3)/np.average(LTD_amounts_2), 1)

    target_value = round(Aminus_values[1] / Aminus_values[0], 1)
    #print(ratio1)
    #print(target_value)
    assert np.allclose(ratio1, target_value, atol=0.1)   
    target_value = round(Aminus_values[2] / Aminus_values[0], 1)
    #print(ratio2)
    #print(target_value)
    assert np.allclose(ratio2, target_value, atol=0.1)
 
    target_value = round(Aminus_values[2] / Aminus_values[1], 1)
    #print(ratio3)
    #print(target_value)
    assert np.allclose(ratio3, target_value, atol=0.1)

sys.exit(0)