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
            }
        )

    for i, neuron in enumerate(cortex):
        nest.SetStatus([neuron], {"fiber_id": i})

    spikedetector = nest.Create("spike_detector")
    nest.Connect(cortex, spikedetector)

    nest.Simulate(sim_len)

    dSD = nest.GetStatus(spikedetector, keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]

    return evs, ts


evs, ts = run_simulation(n=100)

pylab.scatter(ts, evs, marker='.')
pylab.show()
