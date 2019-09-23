from __future__ import print_function
import nest
import numpy as np
import pylab

nest.Install("albertomodule")


def run_simulation(trial_len=1000, sim_len=1000, target=0.0, prism=0.0, n=1):
    nest.ResetKernel()

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


evs, ts = run_simulation(1000, 2000, n=400)

pylab.scatter(ts, evs, marker='.')
pylab.show()
