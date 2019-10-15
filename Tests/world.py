import numpy as np
import nest
from itertools import accumulate
import matplotlib.pyplot as plt

import trajectories

nest.Install("albertomodule")


def run_simulation(trial_len=300, n_trials=1, prism=0.0, n=400):
    nest.ResetKernel()
    trajectories.save_file(prism, trial_len)

    planner = nest.Create(
        "planner_neuron",
        n=n,
        params={
            "trial_length": trial_len,
            "target": 0.0,
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
            "rbf_sdev": 15.0,
            "baseline_rate": 10.0,
            }
        )

    for i, neuron in enumerate(cortex):
        nest.SetStatus([neuron], {"joint_id": i // (n//4),
                                  "fiber_id": i % (n//4)})

    nest.Connect(planner, cortex, 'one_to_one')

    spikedetector = nest.Create("spike_detector")
    nest.Connect(cortex, spikedetector)

    nest.Simulate(trial_len * n_trials)

    dSD = nest.GetStatus(spikedetector, keys="events")[0]
    evs = dSD["senders"]
    ts = dSD["times"]

    return evs, ts


def integrate_torque(evs, ts, j_id, pop_size, pop_offset):
    j_evs = []
    j_ts = []
    pop_size = pop_size // 4

    for ev, t in zip(evs, ts):
        fiber_id = ev - pop_offset - 1
        joint_id = fiber_id // pop_size

        if joint_id == j_id:
            j_evs.append(fiber_id % pop_size)
            j_ts.append(t)

    torques = [2.0*ev / pop_size - 1.0 for ev in j_evs]
    vel = np.array(list(accumulate(torques))) / pop_size
    pos = np.array(list(accumulate(vel))) / pop_size

    return j_ts, torques, vel, pos


def cut_trial(evs, ts, trial_len, trial_i, norm_times=False):
    trial_events = [
        (ev, t)
        for (ev, t) in zip(evs, ts)
        if trial_len*trial_i <= t < trial_len*(trial_i+1)
    ]
    trial_evs, trial_ts = zip(*trial_events)
    if norm_times:
        return np.array(trial_evs), np.array(trial_ts) - trial_len*trial_i
    else:
        return np.array(trial_evs), np.array(trial_ts)


def compute_trajectories(evs, ts, n, trial_len, n_trials):
    trjs = []

    for i in range(n_trials):
        trial_evs, trial_ts = cut_trial(evs, ts, trial_len, i)
        q_ts, qdd, qd, q = integrate_torque(trial_evs, trial_ts, 1, n, n)
        trjs.append([q_ts, q])

    return trjs


def get_final_x(trjs):
    xs = [q[-1] for (q_ts, q) in trjs]
    return np.mean(xs), np.std(xs)


def test_integration():
    prism = 0.0
    duration = 300
    q_in = np.array((10.0, -10.0, -90.0, 170.0))
    q_out = np.array((0.0, prism, 0.0,   0.0))

    q, qd, qdd = trajectories.jtraj(q_in, q_out, duration)
    fig, axs = plt.subplots(6, 4)

    for j in range(4):
        axs[0, j].plot(q[:, j])
        axs[1, j].plot(qd[:, j])
        axs[2, j].plot(qdd[:, j])

    n = 400
    trial_len = 300

    evs, ts = run_simulation(n=n, trial_len=trial_len)

    for j in range(4):
        q_ts, qdd, qd, q = integrate_torque(evs, ts, j, n, n)
        axs[3, j].scatter(q_ts, qdd, marker='.')
        axs[4, j].plot(q_ts, qd)
        axs[5, j].plot(q_ts, q)

    plt.show()


def test_trajectories(n_trials):
    n = 400
    trial_len = 300

    evs, ts = run_simulation(n=n, n_trials=n_trials, trial_len=trial_len)

    trjs = compute_trajectories(evs, ts, n, trial_len, n_trials)

    # mean, std = get_final_x(trjs)
    # print("Mean:", mean)
    # print("Std:", std)

    for q_ts, q in trjs:
        plt.plot(q_ts, q)

    plt.show()


def test_prism(n_trials, prism_values):
    n = 400
    trial_len = 300

    evs, ts = run_simulation(n=n, n_trials=n_trials, prism=0.0)
    trjs = compute_trajectories(evs, ts, n, trial_len, n_trials)

    ref_mean, ref_std = get_final_x(trjs)
    deltas = [0.0]
    stds = [ref_std]

    for prism in prism_values:
        evs, ts = run_simulation(n=n, n_trials=n_trials, prism=prism)
        trjs = compute_trajectories(evs, ts, n, trial_len, n_trials)

        mean, std = get_final_x(trjs)

        delta_x = mean - ref_mean
        deltas.append(delta_x)
        stds.append(std)

    plt.errorbar([0.0] + list(prism_values), deltas, stds)
    plt.show()


def main():
    test_integration()
    test_trajectories(10)
    test_prism(4, [25.0, 50.0, 75.0, 100.0])


if __name__ == '__main__':
    main()
