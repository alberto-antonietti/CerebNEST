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
LTD2 = -1.0e-3

Init_MFDCN = 0.4
Init_PCDCN = -0.005

aux.tic()
''' VARIABLES INITIALIZATION '''
nest.set_verbosity('M_ERROR')
nest.ResetKernel()
nest.SetKernelStatus({'local_num_threads' : CORES,
                      'total_num_virtual_procs' : CORES,
                      'resolution' : 1.0,
                      'overwrite_files' : True})

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
    for n,vti in enumerate(vt2):
        nest.SetStatus([vti],{"vt_num" : n})
    
nest.Connect(InputGen,MF,"one_to_one",{"model": "static_synapse",
                                       "weight": 1.0,
                                       "delay": 1.0 
                                      })
                                    
nest.Connect(ErrorGen,PC,"one_to_one",{"model": "static_synapse",
                                       "weight": 1.0,
                                       "delay": 1.0 
                                      })

PCDCN_conn_param = {"model": "static_synapse",
                    "weight": Init_PCDCN,
                    "delay": 1.0}

nest.Connect(PC,DCN,{"rule" : "fixed_indegree", "indegree" : 40, "multapses" : False, "autapses" : False},PCDCN_conn_param)


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
                                            "Wmin":	 0.0,    # double - Minimal synaptic weight 
                                            "Wmax":      1.0,    # double - Maximal synaptic weight
                                            "vt":        vt2[0]
                                           })

MFDCN_conn_param = {"model": 'stdp_synapse_cosexp',
                    "weight": Init_MFDCN,
                    "delay": 1.0}

for i,DCNind in enumerate(DCN):
    nest.Connect(MF,[DCNind],{'rule': 'fixed_indegree', 'indegree': 12, "multapses": False, "autapses": False},MFDCN_conn_param)

    A=nest.GetConnections(MF,[DCNind])
    nest.SetStatus(A,{"vt_num": i})
    B=nest.GetConnections(PC,[DCNind])
    B=np.array(B)
    if(len(B)>0):
       source=B[:,0]
       nest.Connect(source.tolist(),[(vt2[i])],{"rule" : "all_to_all"},{'model':'static_synapse',
                                                                             'delay':   1.0,
                                                                             'weight':  1.0})
PCVT2_conn = nest.GetConnections(PC,vt2)
print "Number of PC-Vt2 synapses: " + str(len(PCVT2_conn))


MFDCN_conn = nest.GetConnections(MF,DCN)
print("Number of MF-DCN synapses: " + str(len(MFDCN_conn)))


if RECORDING_CELLS:    
    # Create Auxiliary tools
    recdict = [{"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_MF"},
               {"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_PC"},
               {"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_DCN"}]
    spikedetector = nest.Create("spike_detector",3,params=recdict)
    nest.Connect(MF, [spikedetector[0]])
    nest.Connect(PC, [spikedetector[1]])
    nest.Connect(DCN, [spikedetector[2]])
    
# Load input activity on MFs and PCs
MFinput_file = open("MF_12.dat",'r')
for InputGeni in InputGen:
    Spikes_s = MFinput_file.readline()
    Spikes_s = Spikes_s.split()
    Spikes_f = []
    for elements in Spikes_s:
        Spikes_f.append(float(elements))
    nest.SetStatus([InputGeni],{'spike_times' : Spikes_f})
    
PCinput_file = open("PC_40.dat",'r')
for InputGeni in ErrorGen:
    Spikes_s = PCinput_file.readline()
    Spikes_s = Spikes_s.split()
    Spikes_f = []
    for elements in Spikes_s:
        Spikes_f.append(float(elements))
    nest.SetStatus([InputGeni],{'spike_times' : Spikes_f})

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
