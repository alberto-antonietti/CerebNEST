from __future__ import print_function
import nest
# import pylab

nest.Install("albertomodule")


def run_simulation(trial_len, sim_len):
    nest.ResetKernel()

    planner = nest.Create("planner_neuron", params={"trial_length": trial_len})

    spikedetector = nest.Create("spike_detector")
    nest.Connect(planner, spikedetector)

    nest.Simulate(sim_len)

    dSD = nest.GetStatus(spikedetector, keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]

    return evs, ts


def test_poisson():
    evs, ts = run_simulation(1000, 1000)

    print(len(evs))
    print(ts)

    expected = [
        109.7, 193.2, 443.0, 473.7, 538.6, 705.6, 768.2, 978.6,
    ]
    assert(all(round(t, 1) in expected for t in ts))
    # pylab.plot(ts, evs, ".")
    # pylab.show()


def test_periodicity(trial_len=1000, sim_len=3000):
    _, ts = run_simulation(trial_len, sim_len)

    first_trial_ts = [
        round(t, 1) for t in ts
        if t < trial_len
    ]
    assert(all(
        round(t % trial_len, 1) in first_trial_ts
        for t in ts
    ))


test_periodicity(1000, 3000)
test_periodicity(500, 1000)
test_periodicity(100, 500)
test_poisson()
