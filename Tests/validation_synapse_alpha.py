import nest
import matplotlib.pyplot as plt
import numpy as np
import os
import math
from scipy.interpolate import interp1d

def plot_synaptic_weight(pf_mli_conns, data, IO, MLI, GrC, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_MLI_tms, MLI_evs, test_n):
    senders = data['events']['senders']
    targets = data['events']['targets']
    weights = data['events']['weights']
    
    
    for conn in pf_mli_conns:
        #print(conn)
        N = len(MLI)
        source_id = conn.get("source")
        target_id = conn.get("target")
        #print("source_id: ", source_id)
        #print("starget_id: ", target_id)

        indices = np.where((senders == source_id) & (targets == target_id))
        result_weights = weights[indices]
        print(result_weights)

        # Find IO corresponding to MLI
        connection = nest.GetConnections(source = IO, target = MLI[target_id - N-1])
        IO_sender = connection.get("source")
        #print("Corresponding IO: ", IO_sender)
        
        # Get the spikes of IO, Gr and MLI
        IO_times = [t for t, e in zip(corrected_cf_tms, cf_evs) if e == IO_sender]
        #print("IO spikes: ", IO_times)
        
        total_MLI_times = [t for t, e in zip(corrected_MLI_tms, MLI_evs) if e == target_id]
        complex_MLI_times = [t for t in total_MLI_times if t-1.0 in IO_times]
        simple_MLI_times = [t for t in total_MLI_times if t not in complex_MLI_times]
        #print("MLI simple spikes: ", simple_MLI_times)
        #print("MLI complex spikes: ", complex_MLI_times)

        total_Gr_times = [t for t, e in zip(corrected_gr_spikes, Gr_evs) if e == source_id]
        #complex_Gr_times = [t for t in total_Gr_times if t in complex_MLI_times]
        simple_Gr_times = [t for t in total_Gr_times ]
        #print("Simple Gr spikes: ", simple_Gr_times)
        #print("Complex Gr spikes: ", complex_Gr_times)
        #assert(simple_Gr_times)

        

    

        fig, ax = plt.subplots(2,1)
        ax[0].plot(data["events"]["times"][indices], result_weights, '--')
        #print(data["events"]["times"][indices])
        #print(result_weights)
        ax[0].set_ylabel("Synaptic weight")

        ax[1].scatter(IO_times, IO_sender*np.ones(len(IO_times)), label='IO')
        ax[1].scatter(simple_MLI_times, target_id*np.ones(len(simple_MLI_times)), label='MLI')
        ax[1].scatter(complex_MLI_times, target_id*np.ones(len(complex_MLI_times)), color="red", label='Complex spikes')
        ax[1].scatter(simple_Gr_times, source_id*np.ones(len(simple_Gr_times)), label='GrCs')
        #ax[1].scatter(complex_MLI_times, source_id*np.ones(len(complex_MLI_times)), color = "white")
        ax[1].legend(fontsize='small')
        ax[1].set_ylim(0, 20)
        ax[1].set_xlabel("Time [ms]")
        ax[1].set_ylabel("Neuron gID")

        subplot_labels = ['A', 'B']
        for i, axs in enumerate(ax):
            axs.text(-0.1, 1.1, subplot_labels[i], transform=axs.transAxes,
                fontsize=16, fontweight='bold', va='top', ha='right')
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            ax[i].spines['bottom'].set_visible(True)
            ax[i].spines['left'].set_visible(True)


        directory = f"figures_plasticity_alpha_normal_synapse/test_{test_n}"
        filename = f"{directory}/synapse_sender_{source_id}_target_{target_id}.png"
        os.makedirs(directory, exist_ok=True)
        plt.tight_layout()
        plt.savefig(filename)
        

        if test_n == 9:
            fig, ax = plt.subplots(3, 1, figsize=(8, 10), sharex=False)
            
            # --- 1. Synaptic weight trace over time ---
            ax[0].plot(data["events"]["times"][indices], result_weights, '*--', label='Synaptic Weight')
            ax[0].set_ylabel("Synaptic weight", fontsize=12)
            ax[0].legend()

            # --- 2. Spike raster ---
            ax[1].scatter(IO_times, IO_sender*np.ones(len(IO_times)), label='IO')
            ax[1].scatter(simple_MLI_times, target_id*np.ones(len(simple_MLI_times)), label='MLI')
            ax[1].scatter(complex_MLI_times, target_id*np.ones(len(complex_MLI_times)), color="red", label='Complex spikes')
            ax[1].scatter(simple_Gr_times, source_id*np.ones(len(simple_Gr_times)), label='GrCs')

            ax[1].set_ylim(0, 5)
            ax[1].set_xlabel("Time [ms]", fontsize=12)
            ax[1].set_ylabel("Neuron gID", fontsize=12)
            ax[1].legend(fontsize='small')

            # --- 3. LTP amount vs spike time difference ---
            diffs = np.diff(result_weights)
            LTP_mask = diffs > 0  

            # Define a Δt window from -200 ms to 0 ms
            spike_time_diffs = np.linspace(-200, 0, len(diffs))  # Assuming uniform mapping
            LTP_values = diffs[LTP_mask]  # Make LTP amounts positive
            LTP_times = spike_time_diffs[LTP_mask]

            ax[2].plot(LTP_times, LTP_values, 'o-', color='black')
            ax[2].set_ylabel('LTP amount', fontsize=12)
            ax[2].set_xlabel('Spike time difference [ms]', fontsize=12)
            ax[2].set_xticks(np.arange(-200, 1, 20))
            ax[2].set_xlim([-200, 0])

            # --- Consistent formatting for all subplots ---
            subplot_labels = ['A', 'B', 'C']
            for i, axs in enumerate(ax):
                axs.text(-0.1, 1.1, subplot_labels[i], transform=axs.transAxes,
                        fontsize=16, fontweight='bold', va='top', ha='right')
                axs.spines['top'].set_visible(False)
                axs.spines['right'].set_visible(False)
                axs.spines['bottom'].set_visible(True)
                axs.spines['left'].set_visible(True)

            # --- Save figure ---
            filename = f"{directory}/LTP_LTD.png"
            os.makedirs(directory, exist_ok=True)
            plt.tight_layout()
            plt.savefig(filename)

            

def plot_LTP(pf_mli_conns, data, weight_matrix, IO, MLI, GrC, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_MLI_tms, MLI_evs, test_n):
    senders = data['events']['senders']
    targets = data['events']['targets']
    weights = data['events']['weights']
    #print(data['events'])
    for idx, conn in enumerate(pf_mli_conns):
        #print(conn)
        N = len(MLI)
        source_id = conn.get("source")
        target_id = conn.get("target")
        #print("source_id: ", source_id)
        #print("starget_id: ", target_id)

        indices = np.where((senders == source_id) & (targets == target_id))
        result_weights = weights[indices]

        # Find IO corresponding to MLI
        connection = nest.GetConnections(source = IO, target = MLI[target_id - N-1])
        IO_sender = connection.get("source")
        #print("Corresponding IO: ", IO_sender)
        
        # Get the spikes of IO, Gr and MLI
        IO_times = [t for t, e in zip(corrected_cf_tms, cf_evs) if e == IO_sender]
        #print("IO spikes: ", IO_times)
        
        total_MLI_times = [t for t, e in zip(corrected_MLI_tms, MLI_evs) if e == target_id]
        complex_MLI_times = [t for t in total_MLI_times if t-1.0 in IO_times]
        simple_MLI_times = [t for t in total_MLI_times if t not in complex_MLI_times]
        #print("MLI simple spikes: ", simple_MLI_times)
        #print("MLI complex spikes: ", complex_MLI_times)

        first_complex_spike = complex_MLI_times[0]
        print(first_complex_spike)
        first_complex_step = int(first_complex_spike/0.1)
        print(first_complex_step)
        total_Gr_times = [t for t, e in zip(corrected_gr_spikes, Gr_evs) if e == source_id]
        complex_Gr_times = [t for t in total_Gr_times if t in complex_MLI_times]
        simple_Gr_times = [t for t in total_Gr_times if t not in complex_Gr_times]
        #print("Simple Gr spikes: ", simple_Gr_times)
        #print("Complex Gr spikes: ", complex_Gr_times)
        #assert(simple_Gr_times)

        result_weights = weight_matrix[idx, :]

        print(type(data["events"]["times"][indices]))
        LTP_indices = np.where(data["events"]["times"][indices] < first_complex_spike)
        print(type(LTP_indices))
 

        print(type(simple_MLI_times))
        simple_MLI_times = np.array(simple_MLI_times)
        simple_MLI_indices = np.where(simple_MLI_times < first_complex_spike)
        print(type(simple_MLI_indices))
        simple_MLI_times = simple_MLI_times[simple_MLI_indices]

        simple_Gr_times = np.array(simple_Gr_times)
        simple_Gr_indices = np.where(simple_Gr_times < first_complex_spike)
        simple_Gr_times = simple_Gr_times[simple_Gr_indices]

        # Plot weight and spikes
        fig, ax = plt.subplots(3,1, sharex = True, figsize =(8,6))

        default_colors = plt.cm.tab10.colors
        green_color = default_colors[2]
        #ax[1].scatter(IO_times, IO_sender*np.ones(len(IO_times)), label='IO')
        ax[0].scatter(simple_MLI_times[:4], target_id*np.ones(len(simple_MLI_times))[:4], label='MLI', color ="orange")
        #ax[1].scatter(complex_MLI_times, target_id*np.ones(len(complex_MLI_times)), color="red", label='Complex spikes')
        ax[0].scatter(simple_Gr_times[:4], source_id*np.ones(len(simple_Gr_times))[:4], label='GrCs', color =green_color)
        #ax[1].scatter(complex_MLI_times, source_id*np.ones(len(complex_MLI_times)), color = "white")
        ax[0].legend(fontsize='small')
        ax[0].set_ylim(0, 20)
        #ax[0].set_xlim(0, first_complex_spike-5)
        ax[0].set_xlabel("Time [ms]")
        ax[0].set_ylabel("Neuron gID")
        ax[0].tick_params(labelbottom=True)
        
        ax[1].plot(np.arange(math.ceil(first_complex_spike-5)), result_weights[:math.ceil(first_complex_spike-5)])
        ax[1].set_ylabel("Synaptic weight")
        ax[1].set_xlabel("Time [ms]")

        ax[1].tick_params(labelbottom=True)
        
        ax[2].plot(np.arange(math.ceil(first_complex_spike-5))[1:], np.diff(result_weights[:math.ceil(first_complex_spike-5)]), 'k--')
        ax[2].set_ylabel("Diff weight")
        ax[2].set_ylim((-1.7,1.2))
        ax[2].set_xlabel("Time [ms]")
        

            
        subplot_labels = ['A', 'B', 'C']
        for i, axs in enumerate(ax):
            axs.text(-0.1, 1.1, subplot_labels[i], transform=axs.transAxes,
                fontsize=16, fontweight='bold', va='top', ha='right')
            ax[i].spines['top'].set_visible(False)
            ax[i].spines['right'].set_visible(False)
            ax[i].spines['bottom'].set_visible(True)
            ax[i].spines['left'].set_visible(True)


        directory = f"figures_plasticity_alpha_normal_synapse/test_{test_n}"
        filename = f"{directory}/synapse_sender_{source_id}_target_{target_id}_LTP.png"
        plt.tight_layout()
        os.makedirs(directory, exist_ok=True)

        plt.savefig(filename)

        if test_n == 9:
            fig, ax = plt.subplots(3,1)

            LTD_idx = np.array(np.where(np.diff(result_weights) < 0.0))
            LTD_values = np.diff(result_weights)[LTD_idx]
            print(type(LTD_idx))
            print(LTD_values)

            ax[0].plot(data["events"]["times"][indices], result_weights, '*--')
            ax[0].set_ylabel("Synaptic weight")
            #ax[1].plot(data["events"]["times"][indices][:-1], np.diff(result_weights))
            
            ax[1].scatter(IO_times, IO_sender*np.ones(len(IO_times)), label='IO')
            ax[1].scatter(simple_MLI_times, target_id*np.ones(len(simple_MLI_times)), label='MLI')
            ax[1].scatter(complex_MLI_times, target_id*np.ones(len(complex_MLI_times)), color="red", label='Complex spikes')
            ax[1].scatter(simple_Gr_times, source_id*np.ones(len(simple_Gr_times)), label='GrCs')
            #ax[1].scatter(complex_MLI_times, source_id*np.ones(len(complex_MLI_times)), color = "white")
            ax[2].legend(fontsize='small')
            ax[1].set_ylim(0, 5)
            ax[1].set_xlabel("Time [ms]")
            ax[1].set_ylabel("Neuron gID")
            ax[2].plot(LTD_idx.flatten(), -LTD_values.flatten())
            ax[2].set_ylabel('LTD amount')

            subplot_labels = ['A', 'B', 'C']
            for i, axs in enumerate(ax):
                axs.text(-0.1, 1.1, subplot_labels[i], transform=axs.transAxes,
                    fontsize=16, fontweight='bold', va='top', ha='right')
                ax[i].spines['top'].set_visible(False)
                ax[i].spines['right'].set_visible(False)
                ax[i].spines['bottom'].set_visible(True)
                ax[i].spines['left'].set_visible(True)

            filename = f"{directory}/LTP_LTD.png"
            os.makedirs(directory, exist_ok=True)
            plt.tight_layout()
            plt.savefig(filename)

def plot_LTD(pf_mli_conns, data, weight_matrix, IO, MLI, GrC, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_MLI_tms, MLI_evs, test_n):
    senders = data['events']['senders']
    targets = data['events']['targets']
    weights = data['events']['weights']
    #print(data['events'])
    for conn in pf_mli_conns:

        #print(conn)
        N = len(MLI)
        source_id = conn.get("source")
        target_id = conn.get("target")
        #print("source_id: ", source_id)
        #print("starget_id: ", target_id)

        indices = np.where((senders == source_id) & (targets == target_id))
        result_weights = weights[indices]

        # Find IO corresponding to MLI
        connection = nest.GetConnections(source = IO, target = MLI[target_id - N-1])
        IO_sender = connection.get("source")
        #print("Corresponding IO: ", IO_sender)
        
        # Get the spikes of IO, Gr and MLI
        IO_times_tot = [t for t, e in zip(corrected_cf_tms, cf_evs) if e == IO_sender]
        #print("IO spikes: ", IO_times)
        
        total_MLI_times = [t for t, e in zip(corrected_MLI_tms, MLI_evs) if e == target_id]
        complex_MLI_times_tot = [t for t in total_MLI_times if t-1.0 in IO_times_tot]
        print('complex times: ', complex_MLI_times_tot)
        simple_MLI_times_tot = [t for t in total_MLI_times if t not in complex_MLI_times_tot]
        #print("MLI simple spikes: ", simple_MLI_times)
        #print("MLI complex spikes: ", complex_MLI_times)

        first_complex_spike = complex_MLI_times_tot[0]
        first_complex_step = int(first_complex_spike/0.1)

        total_Gr_times = [t for t, e in zip(corrected_gr_spikes, Gr_evs) if e == source_id]
        complex_Gr_times = [t for t in total_Gr_times if t in complex_MLI_times_tot]
        simple_Gr_times_tot = [t for t in total_Gr_times if t not in complex_Gr_times]
        #print("Simple Gr spikes: ", simple_Gr_times)
        #print("Complex Gr spikes: ", complex_Gr_times)
        #assert(simple_Gr_times)
        # Plot weight and spikes

        LTD_idx = np.array(np.where(np.diff(result_weights) < 0.0))
        LTD_values = np.diff(result_weights)[LTD_idx]
        #print(type(LTD_idx))
        #print(LTD_values)

        for num_window, cs in enumerate(complex_MLI_times_tot):
            print('num window: ', num_window)
            end = cs + 50
            start = cs - 202.0

            time_window = np.arange(start, end, 0.1)
            print("time_window: ", time_window)
            
            weight_indices = np.where((data["events"]["times"] >= time_window[0]) & (data["events"]["times"] <= time_window[-1]))
            #print(data["events"]["times"])
            weight_times = data["events"]["times"][weight_indices]
            print('weight_times: ', weight_times)
            print('weights: ', result_weights[weight_indices])

            #print('IO_times prima array: ', IO_times)
            IO_times = np.array(IO_times_tot)
            #print('IO_times dopo array: ', IO_times)
            IO_indices = np.where((IO_times >= time_window[0]) & (IO_times <= time_window[-1]))
            IO_times = IO_times[IO_indices]
            print('IO_times: ', IO_times)

            simple_MLI_times = np.array(simple_MLI_times_tot)
            simple_MLI_indices = np.where((simple_MLI_times >= time_window[0]) & (simple_MLI_times <= time_window[-1]))
            simple_MLI_times = simple_MLI_times[simple_MLI_indices]
            print('simple MLI times: ', simple_MLI_times)

            complex_MLI_times = np.array(complex_MLI_times_tot)
            complex_MLI_indices = np.where((complex_MLI_times >= time_window[0]) & (complex_MLI_times <= time_window[-1]))
            complex_MLI_times = complex_MLI_times[complex_MLI_indices]
            print('complex MLI times: ', complex_MLI_times)

            simple_Gr_times = np.array(simple_Gr_times_tot)
            simple_Gr_indices = np.where((simple_Gr_times >= time_window[0]) & (simple_Gr_times <= time_window[-1]))
            simple_Gr_times = simple_Gr_times[simple_Gr_indices]
            print('Gr times: ', simple_Gr_times)

            result_weights = weight_matrix[int(math.floor(start)) : int(math.ceil(end))]
            interpolator = interp1d(np.arange(len(result_weights))+start, result_weights, kind='linear')  # Linear interpolation

            # Compute interpolated values
            values_high_res = interpolator(np.arange(len(result_weights))+start)
            
            fig, ax = plt.subplots(3,1, figsize = (8, 8), sharex=True)
            #ax[0].plot(weight_times, result_weights[weight_indices], linewidth = 4)
            ax[0].scatter(IO_times, IO_sender*np.ones(len(IO_times)), label='IO')
            ax[0].scatter(simple_MLI_times, target_id*np.ones(len(simple_MLI_times)), label='MLI')
            ax[0].scatter(complex_MLI_times, target_id*np.ones(len(complex_MLI_times)), color="red", label='Complex spikes')
            ax[0].scatter(simple_Gr_times, source_id*np.ones(len(simple_Gr_times)), label='GrCs')
            ax[0].set_ylim(0, 8)
            #ax[1].set_xlim(start, end)
            ax[0].set_xlabel("Time [ms]")
            ax[0].set_ylabel("Neuron gID")
            ax[0].legend(fontsize = 'small')
            ax[0].tick_params(labelbottom=True)
            #ax[0].axhline(result_weights[weight_indices][0], xmin=0, xmax=1, color = 'k', linestyle='--')
            #ax[0].axhline(result_weights[weight_indices][-1], xmin=0, xmax=1, color = 'k', linestyle='--')
            #ax[1].plot(data["events"]["times"][indices][:-1], np.diff(result_weights))
            
            ax[1].plot(np.arange(len(result_weights))+start, values_high_res)
            ax[1].set_ylabel("Synaptic weight")
            ax[1].set_xlabel("Time [ms]")
            ax[1].set_ylim((0.14,0.4))
            ax[1].tick_params(labelbottom=True)

            ax[2].plot(np.arange(len(result_weights))[:-1]+start, np.diff(values_high_res), 'k--')
            ax[2].set_ylabel("Diff weight")
            ax[2].set_ylim((-1.7,1.2))
            ax[2].set_xlabel("Time [ms]")
            ax[2].tick_params(labelbottom=True)

            subplot_labels = ['A', 'B', 'C']
            for i, axs in enumerate(ax):
                axs.text(-0.1, 1.1, subplot_labels[i], transform=axs.transAxes,
                    fontsize=16, fontweight='bold', va='top', ha='right')
                ax[i].spines['top'].set_visible(False)
                ax[i].spines['right'].set_visible(False)
                ax[i].spines['bottom'].set_visible(True)
                ax[i].spines['left'].set_visible(True)

            directory = f"figures_plasticity_alpha_normal_synapse/test_{test_n}"
            filename = f"{directory}/window_{num_window}.png"
            os.makedirs(directory, exist_ok=True)
            plt.tight_layout()
            plt.savefig(filename)
        
            

def reformat_weights(weights):
    total_steps = len(weights)
    n_conns = len(weights[0])
    weight_matrix = np.zeros((n_conns, total_steps))
    for timestep in range(total_steps):
        for conn in range(n_conns):
            weight_matrix[conn, timestep] =round(weights[timestep][conn], 3)
    
    return weight_matrix

def plot_synaptic_matrix(pf_mli_conns, data, weight_matrix, IO, MLI, GrC, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_MLI_tms, MLI_evs, test_n):
    senders = data['events']['senders']
    targets = data['events']['targets']
    #weights = data['events']['weights']
    #print(data['events'])
    for idx, conn in enumerate(pf_mli_conns):
        #print(conn)
        N = len(MLI)
        source_id = conn.get("source")
        target_id = conn.get("target")
        #print("source_id: ", source_id)
        #print("starget_id: ", target_id)

        indices = np.where((senders == source_id) & (targets == target_id))
        #print(type(weight_matrix))
        #print(np.shape(weight_matrix))
        result_weights = weight_matrix[idx, :]

        # Find IO corresponding to MLI
        connection = nest.GetConnections(source = IO, target = MLI[target_id - N-1])
        IO_sender = connection.get("source")
        #print("Corresponding IO: ", IO_sender)
        
        # Get the spikes of IO, Gr and MLI
        IO_times = [t for t, e in zip(corrected_cf_tms, cf_evs) if e == IO_sender]
        #print("IO spikes: ", IO_times)
        
        total_MLI_times = [t for t, e in zip(corrected_MLI_tms, MLI_evs) if e == target_id]
        complex_MLI_times = [t for t in total_MLI_times if t-1.0 in IO_times]
        simple_MLI_times = [t for t in total_MLI_times if t not in complex_MLI_times]
        #print("MLI simple spikes: ", simple_MLI_times)
        #print("MLI complex spikes: ", complex_MLI_times)

        total_Gr_times = [t for t, e in zip(corrected_gr_spikes, Gr_evs) if e == source_id]
        complex_Gr_times = [t for t in total_Gr_times if t in complex_MLI_times]
        simple_Gr_times = [t for t in total_Gr_times if t not in complex_Gr_times]
        #print("Simple Gr spikes: ", simple_Gr_times)
        #print("Complex Gr spikes: ", complex_Gr_times)
        #assert(simple_Gr_times)
        # Plot weight and spikes
        directory = f"figures_plasticity_alpha_normal_synapse/test_{test_n}"
        if test_n != 9:
            fig, ax = plt.subplots(2,1)
            ax[0].plot(np.arange(1000), result_weights, '--')
            ax[0].set_ylabel("Synaptic weight")
            ax[1].scatter(IO_times, IO_sender*np.ones(len(IO_times)), label='IO')
            ax[1].scatter(simple_MLI_times, target_id*np.ones(len(simple_MLI_times)), label='MLI')
            ax[1].scatter(complex_MLI_times, target_id*np.ones(len(complex_MLI_times)), color="red", label='Complex spikes')
            ax[1].scatter(simple_Gr_times, source_id*np.ones(len(simple_Gr_times)), label='GrCs')
            ax[1].scatter(complex_MLI_times, source_id*np.ones(len(complex_MLI_times)), color = "white", s = 70)
            ax[1].legend(fontsize='small')
            ax[1].set_ylim(0, 20)
            ax[1].set_xlabel("Time [ms]")
            ax[1].set_ylabel("Neuron gID")

            subplot_labels = ['A', 'B']
            for i, axs in enumerate(ax):
                axs.text(-0.1, 1.1, subplot_labels[i], transform=axs.transAxes,
                    fontsize=16, fontweight='bold', va='top', ha='right')
                ax[i].spines['top'].set_visible(False)
                ax[i].spines['right'].set_visible(False)
                ax[i].spines['bottom'].set_visible(True)
                ax[i].spines['left'].set_visible(True)
                
            filename = f"{directory}/synapse_sender_{source_id}_target_{target_id}.png"
            os.makedirs(directory, exist_ok=True)

        else:
            fig, ax = plt.subplots(3,1)

            ax[0].plot(np.arange(4300), result_weights, '--')
            ax[0].set_ylabel("Synaptic weight")
            ax[1].plot(np.arange(4299), np.diff(result_weights))
            ax[1].set_ylabel('Weight change')
            ax[2].scatter(IO_times, IO_sender*np.ones(len(IO_times)), label='IO')
            ax[2].scatter(simple_MLI_times, target_id*np.ones(len(simple_MLI_times)), label='Simple MLI')
            ax[2].scatter(complex_MLI_times, target_id*np.ones(len(complex_MLI_times)), color="red", label='Complex spikes')
            ax[2].scatter(simple_Gr_times, source_id*np.ones(len(simple_Gr_times)), label='GrCs')
            #ax[2].scatter(complex_MLI_times, source_id*np.ones(len(complex_MLI_times)), color = "white")
            ax[2].legend(fontsize='small')
            ax[2].set_ylim(0, 20)
            ax[2].set_xlabel("Time [ms]")
            ax[2].set_ylabel("Neuron gID")

            subplot_labels = ['A', 'B', 'C']
            for i, axs in enumerate(ax):
                axs.text(-0.1, 1.1, subplot_labels[i], transform=axs.transAxes,
                    fontsize=16, fontweight='bold', va='top', ha='right')
                ax[i].spines['top'].set_visible(False)
                ax[i].spines['right'].set_visible(False)
                ax[i].spines['bottom'].set_visible(True)
                ax[i].spines['left'].set_visible(True)

            filename = f"{directory}/LTP_LTD.png"
            os.makedirs(directory, exist_ok=True)
        plt.tight_layout()
        plt.savefig(filename)
    
def plot_simple_spike(pf_mli_conns, data, weight_matrix, IO, MLI, GrC, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_MLI_tms, MLI_evs, test_n):
        senders = data['events']['senders']
        targets = data['events']['targets']
        #weights = data['events']['weights']
        #print(data['events'])
        for idx, conn in enumerate(pf_mli_conns):
            #print(conn)
            N = len(MLI)
            source_id = conn.get("source")
            target_id = conn.get("target")
            #print("source_id: ", source_id)
            #print("starget_id: ", target_id)

            indices = np.where((senders == source_id) & (targets == target_id))
            #print(type(weight_matrix))
            #print(np.shape(weight_matrix))
            result_weights = weight_matrix[idx, :]

            # Find IO corresponding to MLI
            connection = nest.GetConnections(source = IO, target = MLI[target_id - N-1])
            IO_sender = connection.get("source")
            #print("Corresponding IO: ", IO_sender)
            
            # Get the spikes of IO, Gr and MLI
            IO_times = [t for t, e in zip(corrected_cf_tms, cf_evs) if e == IO_sender]
            #print("IO spikes: ", IO_times)
            
            total_MLI_times = [t for t, e in zip(corrected_MLI_tms, MLI_evs) if e == target_id]
            complex_MLI_times = [t for t in total_MLI_times if t-1.0 in IO_times]
            simple_MLI_times = [t for t in total_MLI_times if t not in complex_MLI_times]
            #print("MLI simple spikes: ", simple_MLI_times)
            #print("MLI complex spikes: ", complex_MLI_times)

            total_Gr_times = [t for t, e in zip(corrected_gr_spikes, Gr_evs) if e == source_id]
            complex_Gr_times = [t for t in total_Gr_times if t in complex_MLI_times]
            simple_Gr_times = [t for t in total_Gr_times if t not in complex_Gr_times]
            #print("Simple Gr spikes: ", simple_Gr_times)
            #print("Complex Gr spikes: ", complex_Gr_times)
            #assert(simple_Gr_times)
            # Plot weight and spikes
            fig, ax = plt.subplots(2,1)
            ax[0].plot(np.arange(1000), result_weights, '--')
            ax[0].set_ylabel("Synaptic weight")
            ax[1].scatter(IO_times, IO_sender*np.ones(len(IO_times)), label='IO')
            ax[1].scatter(simple_MLI_times, target_id*np.ones(len(simple_MLI_times)), label='MLI')
            ax[1].scatter(complex_MLI_times, target_id*np.ones(len(complex_MLI_times)), color="red", label='Complex MLI')
            ax[1].scatter(simple_Gr_times[-1], source_id*np.ones(len(simple_Gr_times))[-1], label='Simple GR')
            ax[1].legend(fontsize='small')
            ax[1].set_ylim(0, 12)
            ax[1].set_xlabel("Time [ms]")
            ax[1].set_ylabel("Neuron gID")

            subplot_labels = ['A', 'B']
            for i, axs in enumerate(ax):
                axs.text(-0.1, 1.1, subplot_labels[i], transform=axs.transAxes,
                    fontsize=16, fontweight='bold', va='top', ha='right')
                ax[i].spines['top'].set_visible(False)
                ax[i].spines['right'].set_visible(False)
                ax[i].spines['bottom'].set_visible(True)
                ax[i].spines['left'].set_visible(True)

            directory = f"figures_plasticity_alpha_normal_synapse/test_{test_n}"
            filename = f"{directory}/synapse_sender_{source_id}_target_{target_id}.png"
            os.makedirs(directory, exist_ok=True)
            plt.tight_layout()
            plt.savefig(filename)

def plot_complex_spike(pf_mli_conns, data, weight_matrix, IO, MLI, GrC, corrected_cf_tms, cf_evs, corrected_gr_spikes, Gr_evs, corrected_MLI_tms, MLI_evs, test_n):
        senders = data['events']['senders']
        targets = data['events']['targets']
        #weights = data['events']['weights']
        #print(data['events'])
        for idx, conn in enumerate(pf_mli_conns):
            #print(conn)
            N = len(MLI)
            source_id = conn.get("source")
            target_id = conn.get("target")
            #print("source_id: ", source_id)
            #print("starget_id: ", target_id)

            indices = np.where((senders == source_id) & (targets == target_id))
            #print(type(weight_matrix))
            #print(np.shape(weight_matrix))
            result_weights = weight_matrix[idx, :]

            # Find IO corresponding to MLI
            connection = nest.GetConnections(source = IO, target = MLI[target_id - N-1])
            IO_sender = connection.get("source")
            #print("Corresponding IO: ", IO_sender)
            
            # Get the spikes of IO, Gr and MLI
            IO_times = [t for t, e in zip(corrected_cf_tms, cf_evs) if e == IO_sender]
            #print("IO spikes: ", IO_times)
            
            total_MLI_times = [t for t, e in zip(corrected_MLI_tms, MLI_evs) if e == target_id]
            complex_MLI_times = [t for t in total_MLI_times if t-1.0 in IO_times]
            simple_MLI_times = [t for t in total_MLI_times if t not in complex_MLI_times]
            #print("MLI simple spikes: ", simple_MLI_times)
            #print("MLI complex spikes: ", complex_MLI_times)

            total_Gr_times = [t for t, e in zip(corrected_gr_spikes, Gr_evs) if e == source_id]
            complex_Gr_times = [t for t in total_Gr_times if t in complex_MLI_times]
            simple_Gr_times = [t for t in total_Gr_times if t not in complex_Gr_times]
            #print("Simple Gr spikes: ", simple_Gr_times)
            #print("Complex Gr spikes: ", complex_Gr_times)
            #assert(simple_Gr_times)
            # Plot weight and spikes
            fig, ax = plt.subplots(2,1)
            ax[0].plot(np.arange(1000), result_weights, '--')
            ax[0].set_ylabel("Synaptic weight")
            ax[1].scatter(IO_times, IO_sender*np.ones(len(IO_times)), label='IO')
            ax[1].scatter(simple_MLI_times, target_id*np.ones(len(simple_MLI_times)), label='Simple MLI')
            ax[1].scatter(complex_MLI_times, target_id*np.ones(len(complex_MLI_times)), color="red", label='Complex spikes')
            #ax[1].scatter(total_Gr_times[-1], source_id*np.ones(len(total_Gr_times))[-1], color = "red")
            ax[1].legend(fontsize='small')
            ax[1].set_ylim(0, 12)
            ax[1].set_xlabel("Time [ms]")
            ax[1].set_ylabel("Neuron gID")

            subplot_labels = ['A', 'B']
            for i, axs in enumerate(ax):
                axs.text(-0.1, 1.1, subplot_labels[i], transform=axs.transAxes,
                    fontsize=16, fontweight='bold', va='top', ha='right')
                ax[i].spines['top'].set_visible(False)
                ax[i].spines['right'].set_visible(False)
                ax[i].spines['bottom'].set_visible(True)
                ax[i].spines['left'].set_visible(True)

            directory = f"figures_plasticity_alpha_normal_synapse/test_{test_n}"
            filename = f"{directory}/synapse_sender_{source_id}_target_{target_id}.png"
            os.makedirs(directory, exist_ok=True)
            plt.tight_layout()
            plt.savefig(filename)

def correct_spike_times(cf_recorder, MLI_recorder, Gr_recorder):
    cf_spikes = nest.GetStatus(cf_recorder)[0]["events"]
    cf_tms = cf_spikes["times"]
    cf_evs = cf_spikes["senders"]
    #print('evs inside function: ', cf_evs)

    MLI_spikes = nest.GetStatus(MLI_recorder)[0]["events"]
    MLI_tms = MLI_spikes["times"]
    MLI_evs = MLI_spikes["senders"]

    Gr_spikes = nest.GetStatus(Gr_recorder)[0]["events"]
    Gr_tms = Gr_spikes["times"]
    Gr_evs = Gr_spikes["senders"]

    # Climbing fibers
    corrected_cf_tms = []
    for idx, (time, sender) in enumerate(zip(cf_tms, cf_evs)):
        corrected_time = time + sender
        corrected_cf_tms.append(corrected_time)

    corrected_cf_tms = [round(t, 1) for t in corrected_cf_tms]
    #print('cf spikes: (', corrected_cf_tms, '), (', cf_evs, ')')

    # MLI
    corrected_MLI_tms = []
    for idx, (time, sender) in enumerate(zip(MLI_tms, MLI_evs)):
        corrected_time = time + (sender - min(MLI_evs) +1)
        corrected_MLI_tms.append(corrected_time)

    corrected_MLI_tms = [round(t, 1) for t in corrected_MLI_tms]
    #print('MLI spikes: (', corrected_MLI_tms, '), (', MLI_evs, ')\n')
    
    MLI_complex_tms = [spike for spike in corrected_MLI_tms if (spike - 1.0) in corrected_cf_tms]
    MLI_complex_evs = [ev for (spike, ev) in zip(corrected_MLI_tms, MLI_evs) if spike in MLI_complex_tms]
    #print('Complex MLI spikes: (', MLI_complex_tms, '), (', MLI_complex_evs, ')\n')

    # Granule cells
    corrected_gr_spikes = [round(t + s, 1) for t, s in zip(Gr_tms, Gr_evs)]
    assert(len(corrected_gr_spikes) == len(Gr_tms))

    #print('Gr spikes: (', corrected_gr_spikes, '), (', Gr_evs, ')\n')

    return corrected_cf_tms, cf_evs, corrected_MLI_tms, MLI_evs, MLI_complex_tms, MLI_complex_evs, corrected_gr_spikes, Gr_evs
