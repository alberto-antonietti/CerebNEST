import nest
import json
from validation_synapse import plot_synaptic_weight, reformat_weights, plot_synaptic_matrix, correct_spike_times, plot_LTD
import sys

TEST_NUMBER = 9
f = open('eglif_params.json')
params = json.load(f)
f.close()

IO_params = params["EGLIF_model"]["IO"]
GR_params = params["EGLIF_model"]["Granule"]
PC_params = params["EGLIF_model"]["Purkinje"]

nest.ResetKernel()
nest.Install("cerebellar_module")
nest.SetKernelStatus({"resolution": 1.0})

N = 1
M = 1
IO = nest.Create("eglif_io_nestml", N) # IO
#nest.SetStatus(IO, IO_params)
PC = nest.Create("eglif_pc_nestml", N) # PC
nest.SetStatus(PC, PC_params)
Gr = nest.Create("eglif_cond_alpha_multisyn", M) # Gr
#nest.SetStatus(Gr, GR_params) # If used, Gr does not spike

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

nest.Connect(simple_spike_gen, Gr, "all_to_all", {"weight": 1000, "receptor_type": 1})
nest.Connect(complex_spike_gen, IO, "all_to_all", {"weight": 1000, "receptor_type": 1})

# Set thresholds of IOs and GrCs to 0.0 to avoid autonomous firing

nest.SetStatus(Gr, {"V_th": 0.0})
nest.SetStatus(IO, {"V_th": 0.0})


# Connections
nest.Connect(IO, PC, "one_to_one", syn_spec={"receptor_type": 5})


wr = nest.Create("weight_recorder")
nest.CopyModel("stdp_synapse_sinexp", "my_stdp_synapse_rec", {"weight_recorder": wr, "Wmax":150,"receptor_type": 1,"t_0":150, "Aplus":0.0, "Aminus":-0.01})
syn_spec_dict = {'synapse_model': "my_stdp_synapse_rec"}
nest.Connect(Gr, PC, syn_spec = syn_spec_dict)

# Devices
cf_recorder = nest.Create("spike_recorder", {"record_to": "memory"})
PC_recorder = nest.Create("spike_recorder",  {"record_to": "memory"})
Gr_recorder = nest.Create("spike_recorder",  {"record_to": "memory"})


nest.Connect(IO, cf_recorder)
nest.Connect(PC, PC_recorder)
nest.Connect(Gr, Gr_recorder)



pf_pc_conns = nest.GetConnections(synapse_model="my_stdp_synapse_rec")

weights = []
for i in range(4300):
    nest.Simulate(1.0)
    weights.append(nest.GetStatus(pf_pc_conns, "weight"))

weight_matrix = reformat_weights(weights)


corrected_cf_tms, cf_evs, corrected_PC_tms, PC_evs, PC_complex_tms, PC_complex_evs, corrected_gr_spikes, Gr_evs = correct_spike_times(cf_recorder, PC_recorder, Gr_recorder)

data = nest.GetStatus(wr)[0]

plot_synaptic_weight(pf_pc_conns, data, IO, PC, Gr, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_PC_tms, PC_evs, TEST_NUMBER)

sys.exit(0)