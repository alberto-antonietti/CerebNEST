import auxiliary_functions as aux
import nest
nest.Install("albertomodule")
import os
import sys
import numpy as np

''' SIMULATION PROPERTIES '''
CORES = int(sys.argv[1])
RECORDING_CELLS = True
PLAST1 = True
PLAST2 = True
PLAST3 = True

''' PROTOCOL SETUP '''
NumTrial = 100
NumAcq1 = 80.0
TrialDuration = 500.0 # ms
US_Onset = 400.0 # ISI = 400 ms
US_Duration = 100.0 # US Duration = 100 ms
CR_Threshold = 5.0 # DCNAvg threshold that has to be reached to produce a CR

LTP1 =  1.5e-2
LTD1 = -1e-1
LTP2 =  1e-5
LTD2 = -1e-6
LTP3 =  1e-7
LTD3 =  1e-6

Init_PFPC = 2.0
Init_MFDCN = 0.05
Init_PCDCN = -0.25

''' VARIABLES INITIALIZATION '''

nest.set_verbosity('M_WARNING')
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

#aux.tic()
# Define the new kinds of neuron
# http://www.nest-simulator.org/cc/iaf_cond_exp/
# http://www.nest-simulator.org/cc/iaf_psc_alpha/
# http://www.nest-simulator.org/cc/iaf_psc_exp/

nest.CopyModel('iaf_cond_exp','granular_neuron')
nest.CopyModel('iaf_cond_exp','purkinje_neuron')
nest.CopyModel('iaf_cond_exp','olivary_neuron')
nest.CopyModel('iaf_cond_exp','nuclear_neuron')

nest.SetDefaults('granular_neuron',{'t_ref' : 1.0,
                                    'C_m' : 2.0,
                                    'V_th' : -40.0,
                                    'V_reset' : -70.0,
                                    'g_L' : 0.2,
                                    'tau_syn_ex' : 0.5,
                                    'tau_syn_in' : 10.0})

nest.SetDefaults('purkinje_neuron',{'t_ref' : 2.0,
                                    'C_m' : 400.0,
                                    'V_th' : -52.0,
                                    'V_reset' : -70.0,
                                    'g_L' : 16.0,
                                    'tau_syn_ex' : 0.5,
                                    'tau_syn_in' : 1.6})

nest.SetDefaults('olivary_neuron',{'t_ref' : 1.0,
                                    'C_m' : 2.0,
                                    'V_th' : -40.0,
                                    'V_reset' : -70.0,
                                    'g_L' : 0.2,
                                    'tau_syn_ex' : 0.5,
                                    'tau_syn_in' : 10.0})

nest.SetDefaults('nuclear_neuron',{'t_ref' : 1.0,
                                    'C_m' : 2.0,
                                    'V_th' : -40.0,
                                    'V_reset' : -70.0,
                                    'g_L' : 0.2,
                                    'tau_syn_ex' : 0.5,
                                    'tau_syn_in' : 10.0})
if PLAST3:
    nest.SetDefaults('nuclear_neuron',{'tau_minus': 30.0})

# Cell numbers
MF_number = 300
GR_number = MF_number*20
PC_number = 72
IO_number = PC_number
DCN_number = int(PC_number/2)


MF = nest.Create("granular_neuron", MF_number)
GR = nest.Create("granular_neuron", GR_number)
PC = nest.Create("purkinje_neuron", PC_number)
IO = nest.Create("olivary_neuron", IO_number)
DCN = nest.Create("nuclear_neuron", DCN_number)

if PLAST1:
    vt=nest.Create("volume_transmitter_alberto",PC_number)
    for n,vti in enumerate(vt):
        nest.SetStatus([vti],{"vt_num" : n})
if PLAST2:
    vt2=nest.Create("volume_transmitter_alberto",DCN_number)
    for n,vti in enumerate(vt2):
        nest.SetStatus([vti],{"vt_num" : n})

recdict2 = {"to_memory": False,
            "to_file":    True,
            "label":     "PFPC_" + str(CORES),
            "senders":    GR,
            "targets":    PC
           }
WeightPFPC = nest.Create('weight_recorder',params=recdict2)

# Define the new kinds of synapses
MFGR_conn_param = {"model": "static_synapse",
                    "weight": {'distribution' : 'uniform', 'low': 0.55, 'high': 0.7}, # -> 0.75 GR fire at 7 Hz
                    "delay": 1.0}

if PLAST3:
    nest.SetDefaults('stdp_synapse',{"tau_plus": 30.0,
                                     "lambda": LTP3,
                                     "alpha": LTD3/LTP3,
                                     "mu_plus": 0.0,  # Additive STDP
                                     "mu_minus": 0.0, # Additive STDP
                                     "Wmax": -0.5,
                                     "weight": Init_PCDCN,
                                     "delay": 1.0})
    PCDCN_conn_param = {"model": "stdp_synapse"}
else:
    PCDCN_conn_param = {"model": "static_synapse",
                        "weight": Init_PCDCN,
                        "delay": 1.0}

# MF-GR excitatory fixed connections - each GR receives 4 connections from 4 random granule cells
nest.Connect(MF,GR,{'rule': 'fixed_indegree', 'indegree': 4, "multapses": False},MFGR_conn_param)
MFGR_conn = nest.GetConnections(MF,GR)

if PLAST1:
    nest.SetDefaults('stdp_synapse_sinexp',{"A_minus":   LTD1,   # double - Amplitude of weight change for depression
                                            "A_plus":    LTP1,   # double - Amplitude of weight change for facilitation 
                                            "Wmin":      0.0,    # double - Minimal synaptic weight 
                                            "Wmax":      4.0,    # double - Maximal synaptic weight
                                            "vt":        vt[0],
                                            "weight_recorder" : WeightPFPC[0]
                                            })
    PFPC_conn_param = {"model":  'stdp_synapse_sinexp',
                       "weight": Init_PFPC,
                       "delay":  1.0}
    # PF-PC excitatory plastic connections - each PC receives the random 80% of the GR
    for i,PCi in enumerate(PC):
        nest.Connect(GR,[PCi],{'rule': 'fixed_indegree', 'indegree': int(0.8*GR_number), "multapses": False},PFPC_conn_param)
        A=nest.GetConnections(GR,[PCi])
        nest.SetStatus(A,{'vt_num': float(i)})
else:
    PFPC_conn_param = {"model":  "static_synapse",
                       "weight": Init_PFPC,
                       "delay":  1.0}
    nest.Connect(GR,PC,{'rule': 'fixed_indegree', 'indegree': int(0.8*GR_number), "multapses": False},PFPC_conn_param)
PFPC_conn = nest.GetConnections(GR,PC)

if PLAST1:
    # IO-PC teaching connections - Each IO is one-to-one connected with each PC
    nest.Connect(IO, vt, {'rule': 'one_to_one'}, {"model": "static_synapse",
                                                  "weight": 1.0,
                                                  "delay": 1.0})
    IOPC_conn = nest.GetConnections(IO,vt)

if PLAST2:
    # MF-DCN excitatory plastic connections - every MF is connected with every DCN
    nest.SetDefaults('stdp_synapse_cosexp',{"A_minus":   LTD2,   # double - Amplitude of weight change for depression
                                            "A_plus":    LTP2,   # double - Amplitude of weight change for facilitation 
                                            "Wmin":      0.0,    # double - Minimal synaptic weight
                                            "Wmax":      0.25,     # double - Maximal synaptic weight
                                            "vt":        vt2[0]
                                            })

    MFDCN_conn_param = {"model": 'stdp_synapse_cosexp',
                        "weight": Init_MFDCN,
                        "delay": 10.0}
    for i,DCNi in enumerate(DCN):
        nest.Connect(MF,[DCNi],'all_to_all',MFDCN_conn_param)
        A=nest.GetConnections(MF,[DCNi])
        nest.SetStatus(A,{'vt_num': float(i)})
else:
    MFDCN_conn_param = {"model":  "static_synapse",
                        "weight": Init_MFDCN,
                        "delay":  10.0}
    nest.Connect(MF,DCN,'all_to_all',MFDCN_conn_param)

MFDCN_conn = nest.GetConnections(MF,DCN)

# PC-DCN inhibitory plastic connections - each DCN receives 2 connections from 2 contiguous PC
count_DCN=0
for P in range(PC_number):
    nest.Connect([PC[P]],[DCN[count_DCN]],'one_to_one', PCDCN_conn_param)
    if PLAST2:
        nest.Connect([PC[P]],[vt2[count_DCN]],'one_to_one', {"model":  "static_synapse",
                                                             "weight": 1.0,
                                                             "delay":  1.0})
    if P%2==1:
        count_DCN+=1
PCDCN_conn = nest.GetConnections(PC,DCN)

CLOSED_LOOP_P = nest.Create("closed_loop_neuron",int(IO_number/2),{'gain' : CR_Threshold,
                                                              'num_dcn' : float(DCN_number),
                                                              'first_dcn' : float(DCN[0]),
                                                              'positive' : True,
                                                              'to_file' : False,
                                                              'protocol' : 1.0,
                                                              'Tstart' : US_Onset,
                                                              'Tstop' : US_Duration,
                                                              'Tduration' : TrialDuration,
                                                              'phase' : NumAcq1 })
CLOSED_LOOP_N = nest.Create("closed_loop_neuron",int(IO_number/2),{'gain' : 0.0,
                                                              'num_dcn' : float(DCN_number),
                                                              'first_dcn' : float(DCN[0]),
                                                              'positive' : False,
                                                              'to_file' : False,
                                                              'protocol' : 1.0,
                                                              'Tstart' : US_Onset,
                                                              'Tstop' : US_Duration,
                                                              'Tduration' : TrialDuration,
                                                              'phase' : NumAcq1 })

nest.SetStatus([CLOSED_LOOP_P[0]],{'to_file' : True})

nest.Connect(DCN,CLOSED_LOOP_P,'all_to_all')
nest.Connect(DCN,CLOSED_LOOP_N,'all_to_all')
nest.Connect(CLOSED_LOOP_P,IO[:int(IO_number/2)],'one_to_one',{'weight' : 100.0})
nest.Connect(CLOSED_LOOP_N,IO[int(IO_number/2):],'one_to_one',{'weight' : 100.0})

nest.PrintNetwork(2)

if RECORDING_CELLS:
    # Create Auxiliary tools
    recdict = [{"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_MF" + str(CORES)},
               {"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_GR" + str(CORES)},
               {"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_PC" + str(CORES)},
               {"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_IO" + str(CORES)},
               {"to_memory": True, "to_file":True, "withgid":True, "withtime":  True, "label":"Spike_Detector_DCN"+ str(CORES)}]
    spikedetector = nest.Create("spike_detector",5,params=recdict)
    nest.Connect(MF,  [spikedetector[0]])
    nest.Connect(GR,  [spikedetector[1]])
    nest.Connect(PC,  [spikedetector[2]])
    nest.Connect(IO,  [spikedetector[3]])
    nest.Connect(DCN, [spikedetector[4]])


Input_generation = nest.Create("spike_generator", MF_number)
nest.Connect(Input_generation,MF,'one_to_one',{'weight': 100.0})

MFinput_file = open("MF_100Trial_EBCC.dat",'r')

for MFi in Input_generation:
    Spikes_s = MFinput_file.readline()
    Spikes_s = Spikes_s.split()
    Spikes_f = []
    for elements in Spikes_s:
        Spikes_f.append(float(elements))
    nest.SetStatus([MFi],{'spike_times' : Spikes_f})

aux.tic()
msd = 1000 # master seed
n_vp = nest.GetKernelStatus('total_num_virtual_procs')
msdrange1 = range(msd, msd+n_vp )
pyrngs = [np.random.RandomState(s) for s in msdrange1]
msdrange2=range(msd+n_vp+1, msd+1+2*n_vp)
nest.SetKernelStatus({'grng_seed' : msd+n_vp,
                      'rng_seeds' : msdrange2})

print("### SIMULATION STARTS ###")
nest.Simulate(TrialDuration*NumTrial)

aux.toc()

sys.exit(0) #Everything went fine
