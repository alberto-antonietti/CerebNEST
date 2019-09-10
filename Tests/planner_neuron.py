from __future__ import print_function
import nest
import numpy as np
import pylab

nest.Install("albertomodule")


def run_simulation(trial_len, sim_len, target=0.0, prism=0.0, n=1):
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


def test_rate():
    evs, ts = run_simulation(1000, 1000, 0.0, 0.0)
    assert(len(evs) == 8)

    evs, ts = run_simulation(1000, 1000, 10.0, 0.0)
    assert(len(evs) == 91)

    evs, ts = run_simulation(1000, 1000, 10.0, 5.0)
    assert(len(evs) == 150)
    # print(len(evs), "spikes")


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


def test_recurrency():
    n = 1000
    delta_t = 10
    trial_t = 300
    evs, ts = run_simulation(trial_t, trial_t, n=n)

    x = np.zeros([trial_t/delta_t, n])

    for (n_id, t) in zip(evs, ts):
        if t >= trial_t:  # FIXME in planner_neuron
            continue
        time_chunk = int(t // delta_t)
        x[time_chunk, n_id-1] += 1

    cor_mat = np.corrcoef(x)

    pylab.pcolor(cor_mat)
    pylab.show()


test_periodicity(1000, 3000)
test_periodicity(500, 1000)
test_periodicity(100, 500)
test_poisson()
test_rate()

test_recurrency()
