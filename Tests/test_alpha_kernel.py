import nest
import json
from validation_synapse_alpha import plot_synaptic_weight, reformat_weights, plot_synaptic_matrix, correct_spike_times, plot_LTD
import sys

TEST_NUMBER = 9
f = open('eglif_params.json')
params = json.load(f)
f.close()

IO_params = params["EGLIF_model"]["IO"]
GR_params = params["EGLIF_model"]["Granule"]
MLI_params = params["EGLIF_model"]["Stellate"]

nest.ResetKernel()
nest.Install("cerebellar_module")
nest.SetKernelStatus({"resolution": 0.1})

N = 1
M = 1
IO = nest.Create("eglif_io_nestml", N) # IO
MLI = nest.Create("eglif_mli", N) # MLI
nest.SetStatus(MLI, MLI_params)
Gr = nest.Create("eglif_cond_alpha_multisyn", M)

simple_spike_times = [
    1.0, 210.0, 420.0, 630.0, 840.0, 1050.0, 1260.0, 1470.0, 1680.0, 1890.0,
    2100.0, 2310.0, 2520.0, 2730.0, 2940.0, 3150.0, 3360.0, 3570.0, 3780.0, 3990.0, 4200.0
]

complex_spike_times = [
    201.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0,
    2200.0, 2400.0, 2600.0, 2800.0, 3000.0, 3200.0, 3400.0, 3600.0, 3800.0, 4000.0, 4200.0
]

simple_spike_gen = nest.Create("spike_generator", {"spike_times": simple_spike_times})
complex_spike_gen = nest.Create("spike_generator", {"spike_times": complex_spike_times})

nest.Connect(simple_spike_gen, Gr, "all_to_all", {"weight": 100, "receptor_type": 1})
nest.Connect(complex_spike_gen, IO, "all_to_all", {"weight": 100, "receptor_type": 1})

# Set thresholds of IOs and GrCs to 0.0 to avoid autonomous firing

nest.SetStatus(Gr, {"V_th": 0.0})
nest.SetStatus(IO, {"V_th": 0.0})



# Connections
nest.Connect(IO, MLI, "one_to_one", syn_spec={"receptor_type": 3})


wr = nest.Create("weight_recorder")
nest.CopyModel("stdp_synapse_alpha", "my_stdp_synapse_rec", {"weight_recorder": wr ,"receptor_type": 1,"weight":50, "Aplus":10, "Aminus":0.0, "tau":50.0})
syn_spec_dict = {'synapse_model': "my_stdp_synapse_rec"}
nest.Connect(Gr, MLI, syn_spec = syn_spec_dict)

# Devices
cf_recorder = nest.Create("spike_recorder", {"record_to": "memory"})
MLI_recorder = nest.Create("spike_recorder",  {"record_to": "memory"})
Gr_recorder = nest.Create("spike_recorder",  {"record_to": "memory"})


nest.Connect(IO, cf_recorder)
nest.Connect(MLI, MLI_recorder)
nest.Connect(Gr, Gr_recorder)




pf_mli_conns = nest.GetConnections(synapse_model="my_stdp_synapse_rec")

weights = []
for i in range(4300):
    nest.Simulate(1.0)
    weights.append(nest.GetStatus(pf_mli_conns, "weight"))
    

weight_matrix = reformat_weights(weights)
print(weight_matrix)


corrected_cf_tms, cf_evs, corrected_PC_tms, MLI_evs, MLI_complex_tms, MLI_complex_evs, corrected_gr_spikes, Gr_evs = correct_spike_times(cf_recorder, MLI_recorder, Gr_recorder)

data = nest.GetStatus(wr)[0]

plot_synaptic_weight(pf_mli_conns, data, IO, MLI, Gr, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_PC_tms, MLI_evs, TEST_NUMBER)
sys.exit(0)