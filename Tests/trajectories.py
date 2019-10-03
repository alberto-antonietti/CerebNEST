from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def jtraj(q0, q1, tv):
    '''
    Compute a joint space trajectory
    :param q_0: Initial joint positions
    :param q_1: Final joint positions
    :oaram tv: Number of time steps
    :return:
    '''

    qd0 = np.zeros(np.size(q0))
    qd1 = np.copy(qd0)

    tscal = 1
    t = (np.arange(0, tv)).T / (tv - 1)

    # compute the polynomial coefficients
    dq = q1 - q0
    A = np.dot(6,   dq) - np.dot(np.dot(3, (qd1 + qd0)),            tscal)
    B = np.dot(-15, dq) + np.dot((np.dot(8, qd0) + np.dot(7, qd1)), tscal)
    C = np.dot(10,  dq) - np.dot((np.dot(6, qd0) + np.dot(4, qd1)), tscal)
    E = np.dot(qd0, tscal)
    F = np.copy(q0)

    tt = np.transpose(np.array([
        t ** 5, t ** 4, t ** 3, t ** 2, t,
        np.ones(np.size(t))
    ]))
    c = np.array([A, B, C, np.zeros(np.size(A)), E, F])
    qt = np.dot(tt, c)

    # compute velocity
    c = np.array([
        np.zeros(np.size(A)),
        np.dot(5, A), np.dot(4, B), np.dot(3, C),
        np.zeros(np.size(A)), E
    ])
    qdt = np.dot(tt, c) / tscal

    # compute acceleration
    c = np.array([
        np.zeros(np.size(A)), np.zeros(np.size(A)),
        np.dot(20, A), np.dot(12, B), np.dot(6, C),
        np.zeros(np.size(A))
    ])
    qddt = np.dot(tt, c) / tscal ** 2

    return [qt, qdt, qddt]


def plot_traj(prism=25.0, duration=300):
    # Joint angles - Initial and final (with and without perturbation)
    q_in = np.array((10.0, -10.0, -90.0, 170.0))  # [Deg] Initial Position

    # [Deg] Final Position (without perturbation)
    q_fin = np.array((0.0,   0.0,   0.0,   0.0))

    # [Deg] Final Position with perturbation
    q_prism = np.array((0.0, prism, 0.0,   0.0))

    q, qd, qdd = jtraj(q_in, q_fin, duration)
    q_p, qd_p, qdd_p = jtraj(q_in, q_prism, duration)

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Joint Trajectories')
    ax1.plot(q, color='black')
    ax1.plot(q_p, '--', color='red')

    ax2.plot(qd, color='black')
    ax2.plot(qd_p, '--', color='red')

    ax3.plot(qdd, color='black')
    ax3.plot(qdd_p, '--', color='red')

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Joint Trajectories')
    ax1.plot(q[:, 1], color='black')
    ax1.plot(q_p[:, 1], '--', color='red')

    ax2.plot(qd[:, 1], color='black')
    ax2.plot(qd_p[:, 1], '--', color='red')

    ax3.plot(qdd[:, 1], color='black')
    ax3.plot(qdd_p[:, 1], '--', color='red')

    plt.show()


def normalize(qdd):
    # Torques
    trs = [qdd[:, i] for i in range(4)]
    amplitudes = np.array([max(abs(tr)) for tr in trs])

    # print(amplitudes)

    for i in range(4):
        qdd[:, i] /= amplitudes[i]

    amplitudes_norm = amplitudes / max(amplitudes)
    return qdd, amplitudes_norm


def save_file(prism=0.25, duration=300, file_name="JointTorques.dat"):
    q_in = np.array((10.0, -10.0, -90.0, 170.0))
    q_out = np.array((0.0, prism, 0.0,   0.0))

    q, qd, qdd = jtraj(q_in, q_out, duration)

    qdd_n, amps = normalize(qdd)

    np.savetxt(file_name, np.vstack([amps, qdd_n]))

    # fig, axs = plt.subplots(4)

    # for i in range(4):
    #     axs[i].plot(qdd[:, i])

    # plt.show()


if __name__ == '__main__':
    prism = 0.25
    duration = 300
    q_in = np.array((10.0, -10.0, -90.0, 170.0))
    q_out = np.array((0.0, prism, 0.0,   0.0))

    q, qd, qdd = jtraj(q_in, q_out, duration)

    qdd_n, amps = normalize(qdd)
    # print(np.vstack([amps, qdd_n]))

    fig, axs = plt.subplots(4)

    for i in range(4):
        axs[i].plot(qdd[:, i])

    plt.show()
