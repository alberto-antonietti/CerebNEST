from __future__ import print_function
import nest
# import pylab

nest.Install("albertomodule")


def test_poisson():
    nest.ResetKernel()

    planner = nest.Create("planner_neuron", params={"rate": 10.0})

    spikedetector = nest.Create("spike_detector")
    nest.Connect(planner, spikedetector)

    nest.Simulate(1000.0)

    dSD = nest.GetStatus(spikedetector, keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]

    expected = [109.7, 193.2, 443.0, 473.7, 538.6, 705.6, 768.2, 978.6]
    assert(all(round(t, 1) in expected for t in ts))

    print(len(evs))
    print(ts)
    # pylab.plot(ts, evs, ".")
    # pylab.show()


test_poisson()
