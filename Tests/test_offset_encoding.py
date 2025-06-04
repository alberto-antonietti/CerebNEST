import nest
import numpy as np
import json
from validation_synapse import plot_synaptic_weight, reformat_weights, plot_synaptic_matrix, correct_spike_times
import sys



f = open('eglif_params.json')
params = json.load(f)
f.close()

IO_params = params["EGLIF_model"]["IO"]
GR_params = params["EGLIF_model"]["Granule"]
PC_params = params["EGLIF_model"]["Purkinje"]

N = 3
M = 5

for trial in range(N):
    TEST_NUMBER = 5 + trial/10
    nest.ResetKernel()
    nest.Install("cerebellar_module")
    nest.SetKernelStatus({"resolution": 0.1})

    # Create the network and tune parameters
    IO = nest.Create("eglif_io_nestml", N) # IO
    nest.SetStatus(IO, IO_params)
    PC = nest.Create("eglif_pc_nestml", N) # PC
    nest.SetStatus(PC, PC_params)
    Gr = nest.Create("eglif_cond_alpha_multisyn", M) # Gr
    #nest.SetStatus(Gr, GR_params) # If used, Gr does not spike

    nest.SetStatus(Gr, {"V_th": -71.0}) # Increase firing of GrCs

    # "Deactivate" two out of three PCs
    for idx, cell in enumerate(IO):
        if idx != trial:
            nest.SetStatus(cell, {"V_th":0.0})
    
    # Connections
    nest.Connect(IO, PC, "one_to_one", syn_spec={"receptor_type": 5})


    wr = nest.Create("weight_recorder")
    nest.CopyModel("stdp_synapse_sinexp", "my_stdp_synapse_rec", {"weight_recorder": wr,"t_0":150,  "Aplus":0.001})

    conn_spec_dict = {'rule': 'fixed_indegree', 'indegree': 4}
    syn_spec_dict = {'synapse_model': "my_stdp_synapse_rec", "receptor_type": 1}
    nest.Connect(Gr, PC, conn_spec = conn_spec_dict, syn_spec = syn_spec_dict)

    pf_pc_conns = nest.GetConnections(synapse_model="my_stdp_synapse_rec")
    IO_pc_conns = nest.GetConnections(source = IO, target = PC)

    # Devices   
    cf_recorder = nest.Create("spike_recorder", {"record_to": "memory"})
    PC_recorder = nest.Create("spike_recorder",  {"record_to": "memory"})
    Gr_recorder = nest.Create("spike_recorder",  {"record_to": "memory"})

    nest.Connect(IO, cf_recorder)
    nest.Connect(PC, PC_recorder)
    nest.Connect(Gr, Gr_recorder)

    # SIMULATION
    weights = []
    for i in range(1000):
        nest.Simulate(1.0)
        weights.append(nest.GetStatus(pf_pc_conns, "weight"))

    weight_matrix = reformat_weights(weights)
    
    corrected_cf_tms, cf_evs, corrected_PC_tms, PC_evs, PC_complex_tms, PC_complex_evs, corrected_gr_spikes, Gr_evs = correct_spike_times(cf_recorder, PC_recorder, Gr_recorder)

    data = nest.GetStatus(wr)[0]

    plot_synaptic_weight(pf_pc_conns, data, IO, PC, Gr, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_PC_tms, PC_evs, TEST_NUMBER) 

    # Filter connections
    matching_conn = nest.GetConnections(source = IO[trial], target = PC)
    corresponding_PC = matching_conn.get("target")

    # Find the plastic connections involving the "stimulated" PC
    plastic_connections = nest.GetConnections(synapse_model = "my_stdp_synapse_rec", source = Gr, target = PC[corresponding_PC - N - 1])

    for idx, connection in enumerate(pf_pc_conns):
        if connection not in plastic_connections:
            assert np.all(np.diff(weight_matrix[idx, :]) >= 0) # Make sure LTD is not happening 

sys.exit(0)