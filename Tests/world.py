import numpy as np
import nest
from itertools import accumulate
import matplotlib.pyplot as plt

import trajectories

nest.Install("albertomodule")


def run_simulation(trial_len=300, sim_len=300, target=0.0, prism=0.0, n=1):
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

    nest.Simulate(sim_len)

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


def plot_integration(evs, ts, n):
    fig, axs = plt.subplots(3, 4, figsize=(18, 10))

    for j in range(4):
        j_ts, qdd, qd, q = integrate_torque(evs, ts, j, n, n)
        print(q[-1])

        axs[0, j].scatter(j_ts, qdd, marker='.')
        axs[1, j].plot(j_ts, qd)
        axs[2, j].plot(j_ts, q)

    plt.show()


def get_final_position(ts, evs, n):
    results = [integrate_torque(evs, ts, j, n, n) for j in range(4)]
    return np.array([q[-1] for (j_ts, qdd, qd, q) in results])


def plot_multiple(trials, n):
    fig, axs = plt.subplots(3, 4, figsize=(18, 10))

    final_pos = []

    for i in range(trials):
        evs, ts = run_simulation(n=n)
        for j in range(4):
            j_ts, qdd, qd, q = integrate_torque(evs, ts, j, n, n)
            print(q[-1])
            final_pos.append(get_final_position(ts, evs, n))

            axs[0, j].scatter(j_ts, qdd, marker='.')
            axs[1, j].plot(j_ts, qd)
            axs[2, j].plot(j_ts, q)

    final_pos = np.array(final_pos)
    print(final_pos)
    print(np.array([
        (np.mean(final_pos[:, i]), np.var(final_pos[:, i]))
        for i in range(4)
    ]))
    plt.show()


def plot_prims():
    n = 400
    n_trials = 5
    # fig, axs = plt.subplots(n_trials)

    start_x = -10.0
    final_qs = []

    for i in range(n_trials):
        evs, ts = run_simulation(n=n, prism=10.0*i)
        j_id = 1
        j_ts, qdd, qd, q = integrate_torque(evs, ts, j_id, n, n)
        # axs[i].plot(q)
        plt.plot(j_ts, q)
        # final_xs.append(q[-1] * 10 + start_x)
        final_qs.append(q[-1])

    final_qs = np.array(final_qs)
    final_xs = abs(start_x) * final_qs / final_qs[0] + start_x
    print("Final xs:", final_xs)
    print(n)
    plt.show()


def prism_deviations(n, n_trials, start_x, prism_values):
    j_id = 1
    avg_qs = []
    std_qs = []

    for prism in [0.0] + list(prism_values):
        final_qs = []
        for i in range(n_trials):
            evs, ts = run_simulation(n=n, prism=float(prism))
            j_ts, qdd, qd, q = integrate_torque(evs, ts, j_id, n, n)
            final_qs.append(q[-1])

        avg_qs.append(np.average(final_qs))
        std_qs.append(np.std(final_qs))

    avg_qs = np.array(avg_qs)
    std_qs = np.array(std_qs)
    final_xs = abs(start_x) * avg_qs / avg_qs[0] + start_x
    std_xs = abs(start_x) * std_qs / avg_qs[0]

    return final_xs, std_xs


if __name__ == '__main__':
    n = 400
    # plot_multiple(10, n)
    # plot_prims()
    # xs = prism_deviations(400, -10.0, [10, 20, 30, 40])
    xs, stds = prism_deviations(400, 5, -10.0, 10*np.arange(4) + 10)
    print("xs:", xs)
    print("stds:", stds)
    plt.errorbar(np.arange(len(xs)), xs, stds)
    # plt.plot(xs)
    plt.show()

    q_in = np.array((10.0, -10.0, -90.0, 170.0))
    q_out = np.array((0.0, 0.0, 0.0, 0.0))
    delta_expected = q_out - q_in
