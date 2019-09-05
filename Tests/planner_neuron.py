from __future__ import print_function
import nest
# import pylab

nest.Install("albertomodule")


def test_poisson():
    nest.ResetKernel()

    planner = nest.Create("planner_neuron", params={"trial_length": 1000})

    spikedetector = nest.Create("spike_detector")
    nest.Connect(planner, spikedetector)

    nest.Simulate(3000.0)

    dSD = nest.GetStatus(spikedetector, keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]

    print(len(evs))
    print(ts)

    expected = [
        109.7, 193.2, 443.0, 473.7, 538.6, 705.6, 768.2, 978.6,
        1109.7, 1193.2, 1443.0, 1473.7, 1538.6, 1705.6, 1768.2, 1978.6,
        2109.7, 2193.2, 2443.0, 2473.7, 2538.6, 2705.6, 2768.2, 2978.6,
    ]
    assert(all(round(t, 1) in expected for t in ts))
    # pylab.plot(ts, evs, ".")
    # pylab.show()


test_poisson()
