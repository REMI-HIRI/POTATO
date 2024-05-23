"""Copyright 2021 Helmholtz-Zentrum für Infektionsforschung GmbH"""

from scipy import signal
import numpy as np


def preprocess_RAW(Force, Distance, input_settings, input_format):

    if input_format['Augment'] == 1 and input_settings['augment_factor'] != '':
        # Augment data
        new_f = []
        new_d = []
        factor = int(input_settings['augment_factor'])

        for line in range(len(Force) - 1):
            f = Force[line]
            d = Distance[line]
            new_f.append(f)
            new_d.append(d)
            f_next = Force[line + 1]
            d_next = Distance[line + 1]
            delta_f = abs(f_next - f)
            delta_d = abs(d_next - d)

            for point in range(1, factor):
                mu = 0
                sigma_f = (1.5 * delta_f) / 3
                sigma_d = (1.5 * delta_d) / 3
                nf = np.random.normal(mu, sigma_f)
                nd = np.random.normal(mu, sigma_d)
                new_point_f = f + ((point + nf) / factor) * delta_f
                new_point_d = d + ((point + nd) / factor) * delta_d
                new_f.append(new_point_f)
                new_d.append(new_point_d)

        Force = np.array(new_f)
        Distance = np.array(new_d)

    if input_format['preprocess'] == 1:
        # Downsample
        Force_ds = Force[::input_settings['downsample_value']]
        Distance_ds = Distance[::input_settings['downsample_value']]

        # Filter
        b, a = signal.butter(input_settings['filter_degree'], input_settings['filter_cut_off'])
        filteredForce = signal.filtfilt(b, a, Force_ds)
        filteredDistance = signal.filtfilt(b, a, Distance_ds)

        Force_Distance = np.column_stack((filteredForce, filteredDistance * 1000))
        Force_Distance_um = np.column_stack((filteredForce, filteredDistance))

    else:
        Force_Distance = np.column_stack((Force, Distance * 1000))
        Force_Distance_um = np.column_stack((Force, Distance))

    return Force_Distance, Force_Distance_um


# creates a dataset from min force threshold to max force value
def trim_data(FD_data, F_min):
    F_trimmed = np.array([])
    PD_trimmed = np.array([])
    F_low = np.array([])

    F_max = np.where(FD_data[:, 0] == max(FD_data[:, 0]))
    fi = F_max[0][0]

    while FD_data[fi, 0] > F_min and fi > 0:
        fi = fi - 1

    if not fi == F_max[0][0]:
        F_trimmed = FD_data[fi:F_max[0][0], 0]
        PD_trimmed = FD_data[fi:F_max[0][0], 1]
        F_low = FD_data[:fi, 0]
    elif fi == 0:
        print('Could not trim this curve, data below minimum force threshold!')

    return F_trimmed, PD_trimmed, F_low


# creates derivatives for Force and Distance of the trimmed datasets
def create_derivative(input_settings, Frequency_value, F_trimmed, PD_trimmed, F_low):
    d_time = 1 / Frequency_value * input_settings['downsample_value'] * input_settings['step_d']

    x = input_settings['step_d']

    derivative_list = []

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

        derivative_list.append([F_value, PD_value, F_dt, PD_dt])
        x = x + input_settings['step_d']

    derivative_array = np.array(derivative_list)

    return derivative_array
