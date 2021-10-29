
from scipy import signal
import numpy as np


def preprocess_RAW(Force, Distance, input_settings):
    # Downsample
    Force_ds = Force[::input_settings['downsample_value']]
    Distance_ds = Distance[::input_settings['downsample_value']]

    # Filter
    b, a = signal.butter(input_settings['filter_degree'], input_settings['filter_cut_off'])
    filteredForce = signal.filtfilt(b, a, Force_ds)
    filteredDistance = signal.filtfilt(b, a, Distance_ds)
    filteredDistance_ready = filteredDistance * 1000

    Force_Distance = np.column_stack((filteredForce, filteredDistance_ready))

    if Force_Distance[0, 1] > Force_Distance[-1, 1]:  # reverse
        Force_Distance = np.flipud(Force_Distance)

    Force_Distance_um = np.column_stack((filteredForce, filteredDistance))

    return Force_Distance, Force_Distance_um


# creates a dataset from min force threshold to max force value
def trim_data(FD_data, F_min):
    global F_trimmed, PD_trimmed, F_low

    F_trimmed = []
    PD_trimmed = []

    F_max = np.where(FD_data[:, 0] == max(FD_data[:, 0]))
    fi = F_max[0][0]

    while FD_data[fi, 0] < FD_data[fi - 10, 0]:
        fi = fi - 1

    while FD_data[fi, 1] < FD_data[fi - 10, 1]:
        fi = fi - 1
    fi0 = fi

    fi = F_max[0][0]
    # print(fi)

    while FD_data[fi, 0] > F_min:
        fi = fi - 1
        # print(Force_Distance[fi, 0])
    F_trimmed = FD_data[fi:fi0, 0]
    PD_trimmed = FD_data[fi:fi0, 1]
    F_low = FD_data[:fi, 0]

    return F_trimmed, PD_trimmed, F_low


# creates derivations for Force and Distance of the trimmed datasets
def create_derivation(input_settings, Frequency_value):
    d_time = 1 / Frequency_value * input_settings['downsample_value'] * input_settings['step_d']

    x = input_settings['step_d']

    derivation_list = []

    while x < len(F_trimmed):
        if PD_trimmed[0] < PD_trimmed[-1]:
            F_value = (F_trimmed[x] + F_trimmed[x - input_settings['step_d']]) / 2
            PD_value = (PD_trimmed[x] + PD_trimmed[x - input_settings['step_d']]) / 2
            delta_PD = PD_trimmed[x] - PD_trimmed[x - input_settings['step_d']]
            delta_F = F_trimmed[x] - F_trimmed[x - input_settings['step_d']]
            F_dt = delta_F / d_time
            PD_dt = delta_PD / d_time
        else:
            F_value = (F_trimmed[x] + F_trimmed[(x - input_settings['step_d'])]) / 2
            PD_value = (PD_trimmed[x] + PD_trimmed[(x - input_settings['step_d'])]) / 2
            delta_PD = PD_trimmed[x] - PD_trimmed[(x - input_settings['step_d'])]
            delta_F = F_trimmed[x] - F_trimmed[(x - input_settings['step_d'])]
            F_dt = (delta_F / d_time) * (-1)
            PD_dt = (delta_PD / d_time) * (-1)

        derivation_list.append([F_value, PD_value, F_dt, PD_dt])
        x = x + input_settings['step_d']

    derivation_array = np.array(derivation_list)

    return derivation_array
