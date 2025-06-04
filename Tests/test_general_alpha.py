import nest
import json
from validation_synapse_alpha import plot_synaptic_weight, reformat_weights, plot_synaptic_matrix, correct_spike_times, plot_LTP
import sys

TEST_NUMBER = 1
f = open('eglif_params.json')
params = json.load(f)
f.close()

IO_params = params["EGLIF_model"]["IO"]
GR_params = params["EGLIF_model"]["Granule"]
MLI_params = params["EGLIF_model"]["Stellate"]

nest.ResetKernel()
nest.Install("cerebellar_module")
nest.SetKernelStatus({"resolution": 0.1})

N = 3
M = 5
IO = nest.Create("eglif_io_nestml", N) # IO
nest.SetStatus(IO, IO_params)
MLI = nest.Create("eglif_mli", N) # PC
nest.SetStatus(MLI, MLI_params)
Gr = nest.Create("eglif_cond_alpha_multisyn", M) # Gr
#nest.SetStatus(Gr, GR_params) # If used, Gr does not spike


# Lower threshold of GrCs to increase their firing rate

nest.SetStatus(Gr, {"V_th": -71.0})


# Connections
nest.Connect(IO, MLI, "one_to_one", syn_spec={"receptor_type": 3})


wr = nest.Create("weight_recorder")
nest.CopyModel("stdp_synapse_alpha", "my_stdp_synapse_rec", {"weight_recorder": wr,"weight":1, "Aplus" : 0.05,"Aminus": -0.1})
conn_spec_dict = {'rule': 'fixed_indegree', 'indegree': 4}
syn_spec_dict = {'synapse_model': "my_stdp_synapse_rec", "receptor_type": 1}
nest.Connect(Gr, MLI, conn_spec = conn_spec_dict, syn_spec = syn_spec_dict)

# Devices
cf_recorder = nest.Create("spike_recorder", {"record_to": "memory"})
MLI_recorder = nest.Create("spike_recorder",  {"record_to": "memory"})
Gr_recorder = nest.Create("spike_recorder",  {"record_to": "memory"})

nest.Connect(IO, cf_recorder)
nest.Connect(MLI, MLI_recorder)
nest.Connect(Gr, Gr_recorder)

pf_mli_conns = nest.GetConnections(synapse_model="my_stdp_synapse_rec")


weights = []
for i in range(1000):
    nest.Simulate(1.0)
    weights.append(nest.GetStatus(pf_mli_conns, "weight"))
    #print(nest.GetStatus(pf_pc_conns, "weight"))
weight_matrix = reformat_weights(weights)
#print(pf_pc_conns)
corrected_cf_tms, cf_evs, corrected_MLI_tms, MLI_evs, MLI_complex_tms, MLI_complex_evs, corrected_gr_spikes, Gr_evs = correct_spike_times(cf_recorder, MLI_recorder, Gr_recorder)

data = nest.GetStatus(wr)[0]

plot_synaptic_weight(pf_mli_conns, data, IO, MLI, Gr, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_MLI_tms, MLI_evs, TEST_NUMBER)
sys.exit(0)