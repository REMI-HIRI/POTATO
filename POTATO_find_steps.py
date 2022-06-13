"""Copyright 2021 Helmholtz-Zentrum für Infektionsforschung GmbH"""

from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import numpy as np
import statistics
from scipy.signal import argrelextrema


# calculates the standard deviation of a dataset
def STD(input_data, column_number):
    dt_STD = statistics.pstdev(input_data[:, column_number])
    return dt_STD


# creates a moving median of the dataset, therefore using slices of len(window_size)
def moving_median(input_data, column_number, window_size):
    # window_half so that the moving median is dependent of values left AND right
    window_half = int(window_size / 2)
    # window vector defines everything from [window_half:-window_half]
    window_vector = range(window_half, len(input_data) - window_half, 1)
    # regions left and right have to be treated differently, everything that is < window_half to the edges is given the same mm value
    window_left = range(len(input_data[:window_half]))
    window_right = range(len(input_data[-window_half:]))
    mov_med = []

    for n in window_left:
        mm_left = np.median(input_data[:window_half + n, column_number])
        mov_med.append(mm_left)

    for n in window_vector:
        mm = np.median(input_data[n - window_half:n + window_half, column_number])
        mov_med.append(mm)

    for n in window_right:
        mm_right = np.median(input_data[-window_half:, column_number])
        mov_med.append(mm_right)

    return mov_med


# sorting the data based on a x times STD threshold (normal distibuted noise vs extreme values from steps)
def cut_off(input_array, column_number, mm, std, z_score):
    # sort values - inside STD region, above STD region and below STD region
    F_values_inside = []
    PD_values_inside = []
    F_dt_inside = []
    PD_dt_inside = []
    inside_indices = []

    F_values_below = []
    PD_values_below = []
    F_dt_below = []
    PD_dt_below = []

    F_values_above = []
    PD_values_above = []
    F_dt_above = []
    PD_dt_above = []

    i = 0
    for n in range(0, len(input_array), 1):
        if input_array[n, column_number] > mm[int(i)] + z_score * std:
            F_dt_above.append(input_array[n, 2])
            F_values_above.append(input_array[n, 0])
            PD_values_above.append(input_array[n, 1])
            PD_dt_above.append(input_array[n, 3])

        elif input_array[n, column_number] < mm[int(i)] - z_score * std:
            F_dt_below.append(input_array[n, 2])
            F_values_below.append(input_array[n, 0])
            PD_values_below.append(input_array[n, 1])
            PD_dt_below.append(input_array[n, 3])

        else:
            F_dt_inside.append(input_array[n, 2])
            F_values_inside.append(input_array[n, 0])
            PD_values_inside.append(input_array[n, 1])
            PD_dt_inside.append(input_array[n, 3])
            inside_indices.append(n)

        i = n * (len(mm) / len(input_array)) - 1

    Above = np.column_stack([F_values_above, PD_values_above, F_dt_above, PD_dt_above])
    Below = np.column_stack([F_values_below, PD_values_below, F_dt_below, PD_dt_below])
    Inside = np.column_stack([F_values_inside, PD_values_inside, F_dt_inside, PD_dt_inside])

    return Above, Inside, Below, inside_indices


# searching for minima in the force derivative to identify unfolding events
def find_steps_F(input_settings, filename_i, Force_Distance, der_arr, orientation):
    global y_vector_F
    global F_mm2_STD2_positive
    global F_mm2_STD2_negative

    results_F = []
    PD_start_F = []

    STD_1 = STD(der_arr, 2)
    F_mm = moving_median(der_arr, 2, input_settings['window_size'])
    Above, Inside, Below, inside_indices_F = cut_off(der_arr, 2, F_mm, STD_1, input_settings['z-score_f'])

    F_mm2_STD2_positive = []
    F_mm2_STD2_negative = []
    n_runs = 1

    while abs(STD_1 - STD(Inside, 2)) / STD_1 > input_settings['STD_diff']:
        F_mm = moving_median(Inside, 2, input_settings['window_size'])
        STD_1 = STD(Inside, 2)
        Above, Inside, Below, inside_indices_F = cut_off(der_arr, 2, F_mm, STD_1, input_settings['z-score_f'])
        n_runs = n_runs + 1

    if STD_1 < 0.05:
        STD_1 = 0.05

    print('STD is', STD_1)

    Above, Inside, Below, inside_indices_F = cut_off(der_arr, 2, F_mm, STD_1, input_settings['z-score_f'])
    F_mm = moving_median(Inside, 2, input_settings['window_size'])

    y_vector_F = []
    last = 0
    for n in range(len(der_arr)):
        if n in inside_indices_F:
            y_vector_F.append(F_mm[n])
            last = n
        else:
            y_vector_F.append(F_mm[last])
            F_mm.insert(n, F_mm[last])

    for i in range(len(F_mm)):
        F_mm2_STD2_positive.append(F_mm[i] + input_settings['z-score_f'] * STD_1)
        F_mm2_STD2_negative.append(F_mm[i] - input_settings['z-score_f'] * STD_1)

    # find the step points
    # for those steps that cross the STD2 threshold -> find the closest 0 values prior/following to the crossing one

    # for local minima

    loc_min = argrelextrema((Below[:, 2]), np.less)

    n_steps = 1

    for k in loc_min[0]:
        F_dt_loc_min = Below[k, 2]
        F_index = np.where(der_arr[:, 2] == F_dt_loc_min)

        # find start and end of the step
        i_start = F_index[0][0]
        i_end = F_index[0][0]
        while der_arr[i_start, 2] < F_mm[int(i_start * len(F_mm) / len(der_arr))] and der_arr[i_start, 2] < der_arr[i_start - 1, 2] and i_start >= 1:
            i_start = i_start - 1
        if i_start == 0:
            i_start = 1
        while der_arr[i_end, 2] < F_mm[int(i_end * len(F_mm) / len(der_arr))] and der_arr[i_end, 2] < der_arr[i_end + 1, 2] and i_end < (len(der_arr) - 2):
            i_end = i_end + 1

        PD_start_F.append(der_arr[i_start, 1])
        dict1 = {
            "filename": filename_i,
            "orientation": orientation,
            "Derivative of": 'Force',
            'step #': n_steps,
            'F1': der_arr[i_start, 0],
            'F2': der_arr[i_end, 0],
            'Fc': (der_arr[i_start, 0] + der_arr[i_end, 0]) / 2,
            'step start': der_arr[i_start, 1],
            'step end': der_arr[i_end, 1],
            'step length': der_arr[i_end, 1] - der_arr[i_start, 1],
        }

        results_F.append(dict1)

        n_steps = n_steps + 1

    return results_F, PD_start_F


# searching for maxima in the distance derivative to identify unfolding events
def find_steps_PD(input_settings, filename_i, Force_Distance, der_arr, orientation):
    global y_vector_PD
    global PD_mm2_STD2_positive
    global PD_mm2_STD2_negative

    results_PD = []
    PD_start_PD = []

    STD_1 = STD(der_arr, 3)
    PD_mm = moving_median(der_arr, 3, input_settings['window_size'])

    Above, Inside, Below, inside_indices_PD = cut_off(der_arr, 3, PD_mm, STD_1, input_settings['z-score_d'])

    PD_mm2_STD2_positive = []
    PD_mm2_STD2_negative = []

    n_runs = 1
    while abs(STD_1 - STD(Inside, 3)) / STD_1 > input_settings['STD_diff']:
        PD_mm = moving_median(Inside, 3, input_settings['window_size'])
        STD_1 = STD(Inside, 3)
        Above, Inside, Below, inside_indices_PD = cut_off(der_arr, 3, PD_mm, STD_1, input_settings['z-score_d'])
        n_runs = n_runs + 1

    if STD_1 < 0.05:
        STD_1 = 0.05

    print('STD is', STD_1)

    Above, Inside, Below, inside_indices_PD = cut_off(der_arr, 3, PD_mm, STD_1, input_settings['z-score_d'])
    PD_mm = moving_median(Inside, 3, input_settings['window_size'])

    y_vector_PD = []
    last = 0
    for n in range(len(der_arr)):
        if n in inside_indices_PD:
            y_vector_PD.append(PD_mm[n])
            last = n
        else:
            y_vector_PD.append(PD_mm[last])
            PD_mm.insert(n, PD_mm[last])

    for i in range(len(PD_mm)):
        PD_mm2_STD2_positive.append(PD_mm[i] + input_settings['z-score_d'] * STD_1)
        PD_mm2_STD2_negative.append(PD_mm[i] - input_settings['z-score_d'] * STD_1)
    # find the step points
    # for those steps that cross the 3*STD2 threshold -> find the closest 0 values prior/following to the crossing one

    loc_max = argrelextrema(Above[:, 3], np.greater)

    n_steps = 1

    for k in loc_max[0]:
        PD_dt_loc_max = Above[k, 3]
        PD_index = np.where(der_arr[:, 3] == PD_dt_loc_max)

        # find start and end of the step
        i_start = PD_index[0][0]
        i_end = PD_index[0][0]

        while der_arr[i_start, 3] > PD_mm[int(i_start * len(PD_mm) / len(der_arr))] and der_arr[i_start - 1, 3] < der_arr[i_start, 3] and i_start >= 1:
            i_start = i_start - 1
        if i_start == 0:
            i_start = 1

        while der_arr[i_end, 3] > PD_mm[int(i_end * len(PD_mm) / len(der_arr))] and der_arr[i_end, 3] > der_arr[i_end + 1, 3] and i_end < (len(der_arr) - 2):
            i_end = i_end + 1

        PD_start_PD.append(der_arr[i_start, 1])

        dict1 = {
            "filename": filename_i,
            "orientation": orientation,
            "Derivative of": 'Distance',
            'step #': n_steps,
            'F1': der_arr[i_start, 0],
            'F2': der_arr[i_end, 0],
            'Fc': (der_arr[i_start, 0] + der_arr[i_end, 0]) / 2,
            'step start': der_arr[i_start, 1],
            'step end': der_arr[i_end, 1],
            'step length': der_arr[i_end, 1] - der_arr[i_start, 1],
        }

        results_PD.append(dict1)

        n_steps = n_steps + 1

    return results_PD, PD_start_PD


# define steps, that were found by Force- and Distance-derivative (used for fitting afterwards)
def find_common_steps(F_steps, PD_steps):
    common_steps = []
    x = 1

    for n in range(0, len(F_steps)):
        F_steps_dict = F_steps[n]
        step_F_middle = (float(F_steps_dict['step start']) + float(F_steps_dict['step end'])) / 2
        for i in range(0, len(PD_steps)):
            PD_steps_dict = PD_steps[i]

            if step_F_middle > PD_steps_dict['step start'] and step_F_middle < PD_steps_dict['step end']:
                new_steps = PD_steps[i]
                new_step_number = {'step #': x}
                new_steps.update(new_step_number)
                common_steps.append(new_steps)
                x += 1

    return common_steps


def calc_integral(area_1, area_2, step_start_d, step_end_d, step_start_f, step_end_f):
    # calculate work from integrals (estimation)
    work_step = area_1 + ((step_end_d - step_start_d) * 0.5 * (step_start_f + step_end_f)) - area_2
    work_in_kT = work_step / 4.114  # 4.114 is the approximate value of kT (Boltzmann constant times temperature) at 298 K

    return work_step, work_in_kT


def save_figure(export_PLOT, timestamp, filename_i, analysis_folder, Force_Distance, derivative_array, F_trimmed, PD_trimmed, PD_start_F, PD_start_PD):
    figure1 = Figure(figsize=(10, 6), dpi=100)
    subplot1 = figure1.add_subplot(221)
    subplot2 = figure1.add_subplot(222)
    subplot3 = figure1.add_subplot(223)
    subplot4 = figure1.add_subplot(224)

    subplot1.set_ylabel("Force (pN)")
    subplot1.set_title("FD-Curve")
    subplot1.scatter(Force_Distance[:, 1], Force_Distance[:, 0], marker='.', s=0.6, linewidths=None, alpha=1)

    legend_elements = [
        Line2D([0], [0], color='red', lw=1),
        Line2D([0], [0], color='green', lw=1)
    ]

    subplot2.set_title("Trimmed FD-Curve - steps marked")
    subplot2.legend(legend_elements, ['Steps found by F-derivative', 'Steps found by D-derivative'])
    subplot2.scatter(PD_trimmed, F_trimmed, marker='.', s=0.6, linewidths=None, alpha=1)

    for i in range(len(PD_start_F)):
        subplot2.axvline(x=PD_start_F[i], ymin=0, ymax=30, color='red', lw=0.5, alpha=0.5)
    for i in range(len(PD_start_PD)):
        subplot2.axvline(x=PD_start_PD[i], ymin=0, ymax=30, color='green', lw=0.5, alpha=0.5)

    subplot3.set_xlabel("Distance (nm)")
    subplot3.set_ylabel("delta Distance (nm/ms)")
    subplot3.set_title("Distance derivative")
    subplot3.plot(derivative_array[:, 1], derivative_array[:, 3])

    subplot3.plot(derivative_array[:, 1], y_vector_PD)
    subplot3.fill_between(derivative_array[:, 1], PD_mm2_STD2_positive, PD_mm2_STD2_negative, color='black', alpha=0.30)

    subplot4.set_xlabel("Distance (nm)")
    subplot4.set_ylabel("delta Force (pN/ms)")
    subplot4.set_title("Force derivative")
    subplot4.plot(derivative_array[:, 1], derivative_array[:, 2])

    subplot4.plot(derivative_array[:, 1], y_vector_F)
    subplot4.fill_between(derivative_array[:, 1], list(F_mm2_STD2_positive), list(F_mm2_STD2_negative), color='black', alpha=0.30)

    if export_PLOT == 1:
        plotname = analysis_folder + "/" + filename_i + "_plot_" + timestamp + ".png"
        figure1.savefig(plotname, dpi=600)
    else:
        pass

    figure1.clf()
