import nest
nest.Install("albertomodule")

VT = nest.Create("volume_transmitter_alberto", 1)
PRE = nest.Create("iaf_cond_exp", 1)
POST = nest.Create("iaf_cond_exp", 1)

conn_param1 = {"model":    'stdp_synapse_sinexp',
               "A_minus": -0.01,   # double - Amplitude of weight change for depression
               "A_plus":   0.01,   # double - Amplitude of weight change for facilitation 
               "Wmin":     0.0,    # double - Minimal synaptic weight 
               "Wmax":     4.0,    # double - Maximal synaptic weight,              
               "weight":   1.0,
               "delay":    1.0}

nest.Connect(PRE,POST,{'rule': 'one_to_one'},conn_param1)
A=nest.GetConnections(PRE,POST)
nest.SetStatus(A,{'vt': VT[0]})

### 1 CORE

CORES = 1
nest.set_verbosity('M_WARNING')
nest.ResetKernel()
nest.SetKernelStatus({'local_num_threads' : CORES,
                      'total_num_virtual_procs' : CORES,
                      'resolution' : 1.0,
                      'overwrite_files' : True})
msd = 1000 # master seed
n_vp = nest.GetKernelStatus('total_num_virtual_procs')
msdrange1 = range(msd, msd+n_vp )
pyrngs = [np.random.RandomState(s) for s in msdrange1]
msdrange2=range(msd+n_vp+1, msd+1+2*n_vp)
nest.SetKernelStatus({'grng_seed' : msd+n_vp,
                      'rng_seeds' : msdrange2})

SPIKES1 = [50.0, 100.0, 110.0, 120.0, 130.0,370.0]
SPIKES2 = [1.0, 200.0, 220.0, 230.0]

PoissonGen1 = nest.Create('spike_generator',params={'spike_times': SPIKES1})
PoissonGen2 = nest.Create('spike_generator',params={'spike_times': SPIKES2})

nest.Connect(PoissonGen1,PRE,'one_to_one',{'weight': 1000.0})
nest.Connect(PoissonGen2,VT,'one_to_one',{'weight': 1000.0})

WEIGHT_1=[]
for i in range(500):
        nest.Simulate(1.0)
        WEIGHT_1.append(nest.GetStatus(A,'weight')[0]))

print(WEIGHT_1)

### 4 CORES

CORES = 4
nest.set_verbosity('M_WARNING')
nest.ResetKernel()
nest.SetKernelStatus({'local_num_threads' : CORES,
                      'total_num_virtual_procs' : CORES,
                      'resolution' : 1.0,
                      'overwrite_files' : True})
msd = 1000 # master seed
n_vp = nest.GetKernelStatus('total_num_virtual_procs')
msdrange1 = range(msd, msd+n_vp )
pyrngs = [np.random.RandomState(s) for s in msdrange1]
msdrange2=range(msd+n_vp+1, msd+1+2*n_vp)
nest.SetKernelStatus({'grng_seed' : msd+n_vp,
                      'rng_seeds' : msdrange2})

SPIKES1 = [50.0, 100.0, 110.0, 120.0, 130.0,370.0]
SPIKES2 = [1.0, 200.0, 220.0, 230.0]

PoissonGen1 = nest.Create('spike_generator',params={'spike_times': SPIKES1})
PoissonGen2 = nest.Create('spike_generator',params={'spike_times': SPIKES2})

nest.Connect(PoissonGen1,PRE,'one_to_one',{'weight': 1000.0})
nest.Connect(PoissonGen2,VT,'one_to_one',{'weight': 1000.0})

WEIGHT_4=[]
for i in range(500):
        nest.Simulate(1.0)
        WEIGHT_4.append(nest.GetStatus(A,'weight')[0]))

print(WEIGHT_4)

        
sys.exit(0) #Everything went fine
