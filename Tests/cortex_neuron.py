from __future__ import print_function
import nest
import numpy as np
import pylab
import trajectories

nest.Install("albertomodule")


# nest.SetKernelStatus({'resolution': 1.0})


# 300ms to comply with JointTorques_ideal.dat
def run_simulation(trial_len=300, sim_len=900, target=0.0, prism=0.0, n=1):
    nest.ResetKernel()
    trajectories.save_file(prism, trial_len)

    planner = nest.Create(
        "planner_neuron",
        n=n,
        params={
            "trial_length": trial_len,
            "target": target,
            "prism_deviation": prism,
            "baseline_rate": 10.0,
            "gain_rate": 10.0,
            }
        )

    cortex = nest.Create(
        "cortex_neuron",
        n=n,
        params={
            "trial_length": trial_len,
            "fibers_per_joint": n//4,
            "rbf_sdev": 20.0
            }
        )

    for i, neuron in enumerate(cortex):
        nest.SetStatus([neuron], {"joint_id": i // (n//4),
                                  "fiber_id": i % (n//4)})

    # nest.SetStatus(cortex[:1], {"to_file": True})

    # parrot = nest.Create("parrot_neuron", n)
    # poisson = nest.Create('spike_generator')
    # nest.Connect(planner, cortex, "one_to_one")
    # nest.Connect(poisson, parrot, 'all_to_all')
    nest.Connect(planner, cortex, 'one_to_one')

    spikedetector = nest.Create("spike_detector")
    nest.Connect(cortex, spikedetector)

    nest.Simulate(sim_len)

    dSD = nest.GetStatus(spikedetector, keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]

    return evs, ts


def poisson_input(trial_len=500, sim_len=1000, n=400):
    nest.ResetKernel()
    trajectories.save_file(0.0, trial_len)

    cortex = nest.Create(
        "cortex_neuron",
        n=n,
        params={
            "trial_length": trial_len,
            "fibers_per_joint": n//4,
            "rbf_sdev": 20.0
            }
        )

    for i, neuron in enumerate(cortex):
        nest.SetStatus([neuron], {"joint_id": i // (n//4),
                                  "fiber_id": i % (n//4)})

    nest.SetStatus(cortex[:1], {"to_file": True})

    poisson = nest.Create(
        'poisson_generator',
        n=n,
        params={"rate": 20.0}
    )
    nest.Connect(
        poisson, cortex,
        'one_to_one', syn_spec={'weight': 1.0}
    )

    poisson_minus = nest.Create(
        'poisson_generator',
        n=n,
        params={"rate": 10.0}
    )
    nest.Connect(
        poisson_minus, cortex,
        'one_to_one', syn_spec={'weight': -1.0}
    )

    spikedetector = nest.Create("spike_detector")
    nest.Connect(cortex, spikedetector)

    nest.Simulate(sim_len)

    dSD = nest.GetStatus(spikedetector, keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]

    pylab.scatter(ts, evs, marker='.')
    pylab.show()

    return evs, ts


run_simulation(n=400)

poisson_input()
