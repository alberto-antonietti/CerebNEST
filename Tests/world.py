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


n = 400
evs, ts = run_simulation(n=n)

# pylab.scatter(ts, evs, marker='.')
# pylab.show()


fig, axs = plt.subplots(3, 4, figsize=(18, 10))

for j in range(4):
    j_ts, trq, vel, pos = integrate_torque(evs, ts, j, n, n)

    axs[0, j].scatter(j_ts, trq, marker='.')
    axs[1, j].plot(j_ts, vel)
    axs[2, j].plot(j_ts, pos)

plt.show()
