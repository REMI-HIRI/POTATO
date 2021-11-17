"""Copyright 2021 Helmholtz-Zentrum für Infektionsforschung GmbH"""

from tkinter import filedialog
import matplotlib.pyplot as plt
import h5py
import numpy as np
import time
import pylab
from scipy.optimize import curve_fit
from scipy.integrate import simps
import pandas as pd

from POTATO_preprocessing import preprocess_RAW


# asks for a constant force file, opens and preprocesses the data
def get_constantF(input_settings, input_format, input_constantF):
    global analysis_folder

    # open constant force data in csv or h5 format
    file = filedialog.askopenfilename()
    # get the directory of the selected file and create a path for the analysis folder
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    path_split = file.split('/')
    directory = '/'.join(path_split[:-1])

    # check selected input format
    if input_format['CSV'] == 1:
        with open(file, "r") as f:
            df = pd.read_csv(file)
            # access the raw data
            Force_1x = df.to_numpy()[:, 0]
            Distance_1x = df.to_numpy()[:, 1]
            # accessing the data frequency from user input
            Frequency_value = input_settings['data_frequency']
            Force_Distance, Force_Distance_um = preprocess_RAW(Force_1x, Distance_1x, input_settings)

            filename = path_split[-1][:-4]
            analysis_folder = str(directory + '/Analysis_constantF_' + filename + '_' + timestamp)

    else:
        with h5py.File(file, "r") as f:
            filename = path_split[-1][:-3]
            analysis_folder = str(directory + '/Analysis_constantF_' + filename + '_' + timestamp)

            if input_format['HF'] == 1:
                # opening file and get the raw h5 values
                Force_1x = f.get("Force HF/Force 1x")
                Distance_1x = f.get("Distance/Piezo Distance")
                # accessing the data frequency from the h5 file
                Frequency_value = Force_1x.attrs['Sample rate (Hz)']
                Force_Distance, Force_Distance_um = preprocess_RAW(Force_1x, Distance_1x, input_settings)

            elif input_format['LF'] == 1:
                load_force = f.get("Force LF/Force 1x")
                Force_1x = load_force[:]['Value'][:]
                load_distance = f.get("Distance/Distance 1")[:]
                Distance_1x = load_distance['Value'][:]
                Force_Distance = np.column_stack((Force_1x, Distance_1x))

                # calculating the data frequency based on start- and end-time of the measurement
                size_F_LF = len(Force_1x)
                stop_time_F_LF = load_force.attrs['Stop time (ns)']
                start_time_F_LF = load_force.attrs['Start time (ns)']
                Frequency_value = size_F_LF / ((stop_time_F_LF - start_time_F_LF) / 10**9)

                Force_Distance, Force_Distance_um = preprocess_RAW(Force_1x, Distance_1x, input_settings)

    return Force_Distance, Force_Distance_um, Frequency_value, filename, analysis_folder, timestamp


# function for a gaussian to fit onto the constant force data
def gauss(x, mu, sigma, A):
    return A * pylab.exp(-(x - mu)**2 / 2 / sigma**2)


# depending on how many gaussians should be fitted, this function summs them
def sum_gauss(x, *args):
    sum = 0
    fit_number = int(len(args) / 3)
    for i in range(fit_number):
        sum = sum + gauss(x, args[i], args[i + fit_number], args[i + 2 * fit_number])
    return sum


# for the fit, values have to be guessed so that the fit can be optimized easily
def expected_modal(input_constantF):
    # get the input values for the fit guesses (specified in POTATO_config and the GUI)
    mu = input_constantF['Mean'].split(',')
    mu_map = map(int, mu)
    mu = list(mu_map)
    sigma = input_constantF['STD'].split(',')
    sigma_map = map(int, sigma)
    sigma = list(sigma_map)
    A = input_constantF['Amplitude'].split(',')
    A_map = map(int, A)
    A = list(A_map)

    return mu, sigma, A


# display constant force data in the GUI to perform guesses on the fit
def display_constantF(FD, FD_um, frequency, input_settings, input_constantF):
    global Figure_constantF

    """set constant axes-values"""
    x_min = input_constantF['x min']
    x_max = input_constantF['x max']
    y_min = input_constantF['y min']
    y_max = input_constantF['y max']

    filteredDistance = FD_um[:, 1]
    time_Distance = range(0, len(filteredDistance) * input_settings['downsample_value'], input_settings['downsample_value'])
    time_vector_D = range(0, len(filteredDistance) * input_settings['downsample_value'], input_settings['downsample_value'])

    time_D = []
    for t in time_vector_D:
        time_D.append(t / int(frequency))

    Distance_time = []
    for t in time_Distance:
        Distance_time.append(t / int(frequency))

    # create a Figure
    Figure_constantF, (ax_scatter, ax_histy) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]}, sharey=True)

    ax_scatter.set_xlabel('time, s')
    ax_scatter.set_ylabel('Distance, nm')
    ax_scatter.tick_params(direction='in', top=True, right=True)
    ax_scatter.set_xlim(x_min, x_max)
    ax_scatter.set_ylim(y_min, y_max)

    # the scatter plot:
    filteredDistance_ready = FD[:, 1]
    ax_scatter.plot(time_D, filteredDistance_ready)

    binwidth = 0.1
    bins = np.arange(min(filteredDistance_ready), max(filteredDistance_ready), binwidth)

    ax_histy.tick_params(direction='in', labelleft=False)
    ax_histy.set_xlabel('counts')
    ax_histy.set_ylim(y_min, y_max)
    hist_D = ax_histy.hist(filteredDistance_ready, bins=bins, orientation='horizontal', color='k')

    return Figure_constantF, hist_D, filteredDistance_ready


# open constant force data and try to fit gaussians on the distance distribution
def fit_constantF(hist_D, FD, filteredDistance_ready, frequency, input_settings, input_constantF, filename_i, timestamp):

    hist_der_value = []
    hist_der_position = []

    for i in range(1, len(hist_D[0])):
        der = (hist_D[0][i] - hist_D[0][i - 1]) / (hist_D[1][i] - hist_D[1][i - 1])
        hist_der_value.append(der)
        hist_der_position.append((hist_D[1][i] + hist_D[1][i - 1]) / 2)

    mu, sig, A = expected_modal(input_constantF)
    p0 = np.array([mu, sig, A])
    params, cov = curve_fit(sum_gauss, hist_D[1][0:-1], hist_D[0], p0)
    fit_number = int(len(params) / 3)
    peak_means = []
    sigmas = []
    peak_amplitudes = []
    G = []

    for i in range(0, fit_number):
        peak_means.append(params[i])
        sigmas.append(params[i + fit_number])
        peak_amplitudes.append(params[i + 2 * fit_number])
        G_x = gauss(hist_D[1][0:-1], params[i], params[i + fit_number], params[i + 2 * fit_number])
        G.append(G_x)

    all_param = peak_means + sigmas + peak_amplitudes
    tpl_all_param = tuple(all_param)
    Gsum = sum_gauss(hist_D[1][0:-1], *tpl_all_param)

    figure2, subplot2 = plt.subplots(1, 1)

    binwidth = 0.1
    bins = np.arange(min(filteredDistance_ready), max(filteredDistance_ready), binwidth)
    subplot2.hist(filteredDistance_ready, bins=bins)

    subplot2.plot(hist_D[1][0:-1], Gsum, 'r')
    subplot2.plot(peak_means, peak_amplitudes, color='k', linewidth=0, marker='.', markersize=10)

    Areas = []
    for i in G:
        subplot2.plot(hist_D[1][0:-1], i, 'k', linestyle='dashed')
        Area_gauss = simps(i, hist_D[1][0:-1])
        Areas.append(Area_gauss)

    fractions = []
    for i in range(len(Areas)):
        Area_x = Areas[i]
        # Area_rest = sum(Areas[:i]) + sum(Areas[i + 1:])
        Area_all = sum(Areas[:])
        A_fraction = Area_x / Area_all

        fractions.append(A_fraction)

    # export csv containing all the computed values
    # total_results_export=pd.DataFrame(zip(X_step_data,F_step_data, F, x_WLC, X, x_FJC), columns=['Distance data','Force with step','Force model','WLC distance','WLC + FJC distance', 'FJC distance'])
    smooth_export = pd.DataFrame(zip(FD[:, 0], filteredDistance_ready), columns=['Force [pN]', 'Distance [nm]'])
    smooth_export.to_csv((analysis_folder + "/" + filename_i + "_" + timestamp + '_smooth.csv'), index=False, header=True)

    # histogram & gaussians
    histogram_export = pd.DataFrame(zip(hist_D[1][0:-1], hist_D[0], Gsum, G), columns=['Distance, nm', 'Counts', 'Gaussian sum', 'Gaussians'])
    histogram_export.to_csv((analysis_folder + "/" + filename_i + "_" + timestamp + '_histogram.csv'), index=False, header=True)

    # gaussian parameters + areas
    table_export = pd.DataFrame(zip(peak_means, sigmas, peak_amplitudes, Areas, fractions), columns=['mean', 'sigma', 'amplitude', 'Area', 'fraction'])
    table_export.to_csv((analysis_folder + "/" + filename_i + "_" + timestamp + '_table.csv'), index=False, header=True)
    # plots
    plotname2 = analysis_folder + "/" + filename_i + "_histogram_fitted" + ".png"
    figure2.savefig(plotname2)

    plotname_figure_constant_f = analysis_folder + "/" + filename_i + "_visualization" + ".png"
    Figure_constantF.savefig(plotname_figure_constant_f)

    return figure2
