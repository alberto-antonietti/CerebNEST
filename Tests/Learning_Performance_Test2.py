import sys
import numpy as np
import auxiliary_functions as aux
import nest
import os
import errno
import glob

nest.Install("albertomodule")

# Cell numbers
MF_number = 12
PC_number = 40
DCN_number = 2  

''' SIMULATION PROPERTIES '''
CORES = int(sys.argv[1])
RECORDING_CELLS = True
RECORDING_WEIGHTS = True
PLAST2 = True

''' PROTOCOL SETUP '''
LTP2 =  1.0e-3
LTD2 = -1.0e-2

Init_PCDCN = -0.005
Init_MFDCN = 0.4

aux.tic()
''' VARIABLES INITIALIZATION '''
nest.set_verbosity('M_ERROR')
nest.ResetKernel()
nest.SetKernelStatus({'local_num_threads' : CORES,
                      'total_num_virtual_procs' : CORES,
                      'resolution' : 1.0,
                      'overwrite_files' : True})
if CORES > 1:
    nest.SetNumRecProcesses(0)
msd = 1000 # master seed
msdrange1 = range(msd, msd+CORES )
pyrngs = [np.random.RandomState(s) for s in msdrange1]
msdrange2=range(msd+CORES+1, msd+1+2*CORES)
nest.SetKernelStatus({'grng_seed' : msd+CORES,
                      'rng_seeds' : msdrange2})

# Define the new kinds of neuron
nest.CopyModel('parrot_neuron','mossy_fiber')
nest.CopyModel('parrot_neuron','purkinje_neuron')
nest.CopyModel('iaf_cond_exp','deep_cerebellar_nuclei')


nest.SetDefaults('deep_cerebellar_nuclei',{'t_ref' : 1.0,
                                           'C_m' : 2.0,
                                           'V_th' : -40.0,
                                           'V_reset' : -70.0,
                                           'g_L' : 0.2,
                                           #'tau_m' : 10.0,
                                           'tau_syn_ex' : 0.5,
                                           'tau_syn_in' : 10.0})

MF = nest.Create("mossy_fiber", MF_number)
PC = nest.Create("purkinje_neuron", PC_number)
DCN = nest.Create("deep_cerebellar_nuclei", DCN_number)

InputGen = nest.Create("spike_generator", MF_number)
ErrorGen = nest.Create("spike_generator", PC_number)
if RECORDING_WEIGHTS:
    recdict2 = {"to_memory": False,
         	"to_file":    True,
		"label":     "MFDCN",
		"senders":    MF,
		"targets":    DCN
	       }
    WeightMFDCN = nest.Create('weight_recorder',params=recdict2)
    
if PLAST2:
    vt2=nest.Create("volume_transmitter_alberto",DCN_number)
    
# Connection 0
nest.Connect(InputGen,MF,"one_to_one",{"model": "static_synapse",
                                    "weight": 1.0,
                                    "delay": 1.0 
                                    })
                                    
nest.Connect(ErrorGen,PC,"one_to_one",{"model": "static_synapse",
                                       "weight": 1.0,
                                       "delay": 1.0 
                                      })

# Connection 6
if RECORDING_WEIGHTS:
    nest.SetDefaults('stdp_synapse_cosexp',{"A_minus":   LTD2,   # double - Amplitude of weight change for depression
    					    "A_plus":    LTP2,   # double - Amplitude of weight change for facilitation 
       					    "Wmin":      0.0,    # double - Minimal synaptic weight 
					    "Wmax":      1.0,    # double - Maximal synaptic weight
                                            "vt":        vt2[0],
              				    "weight_recorder": WeightMFDCN[0]
					   })
else:
    nest.SetDefaults('stdp_synapse_cosexp',{"A_minus":   LTD2,   # double - Amplitude of weight change for depression
 					    "A_plus":    LTP2,   # double - Amplitude of weight change for facilitation 
					    "Wmin":      0.0,    # double - Minimal synaptic weight 
					    "Wmax":      1.0,    # double - Maximal synaptic weight
					    "vt":        vt2[0]})	
    MFDCN_conn_param = {"model":  'stdp_synapse_cosexp',
                       "weight": Init_MFDCN,
                       "delay":  1.0,
                       }
    nest.Connect(MF,DCN,{'rule': "fixed_indegree", 'indegree':12, "multapses": False, "autapses": False},MFDCN_conn_param)

MFDCN_conn = nest.GetConnections(MF,DCN)
print("Number of GR-PC synapses: " + str(len(PFPC_conn)))

# Connection 7
if PLAST1:
    # IO-PC teaching connections
    nest.Connect(IO,vt,'all_to_all',{"model": "static_synapse",
                                     "weight": 1.0,
                                     "delay": 1.0})
    IOPC_conn = nest.GetConnections(IO,vt)
    print("Number of IO-PC (volume_transmitter) synapses: " + str(len(IOPC_conn)))

if RECORDING_CELLS:    
    # Create Auxiliary tools
    recdict = [{"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_GR"},
               {"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_PC"},
               {"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_IO"}]
    spikedetector = nest.Create("spike_detector",3,params=recdict)
    nest.Connect(GR, [spikedetector[0]])
    nest.Connect(PC, [spikedetector[1]])
    nest.Connect(IO, [spikedetector[2]])
    
# Load input activity on GRs
GRinput_file = open("GR_65600.dat",'r')
for InputGeni in InputGen:
    Spikes_s = GRinput_file.readline()
    Spikes_s = Spikes_s.split()
    Spikes_f = []
    for elements in Spikes_s:
        Spikes_f.append(float(elements))
    nest.SetStatus([InputGeni],{'spike_times' : Spikes_f})
    
nest.SetStatus(ErrorGen,{'spike_times' : [8.0, 98.0, 298.0, 308.0, 318.0, 498.0, 598.0, 698.0, 798.0, 997.0]})

aux.toc()

msd = 1000 # master seed
n_vp = nest.GetKernelStatus('total_num_virtual_procs')
msdrange1 = range(msd, msd+n_vp )
pyrngs = [np.random.RandomState(s) for s in msdrange1]
msdrange2=range(msd+n_vp+1, msd+1+2*n_vp)
nest.SetKernelStatus({'grng_seed' : msd+n_vp,
                      'rng_seeds' : msdrange2})

print("### SIMULATION STARTS ###") 
aux.tic()
nest.Simulate(1001)
aux.toc()

sys.exit(0) #Everything went fine 
