"""Copyright 2021 Helmholtz-Zentrum für Infektionsforschung GmbH"""

""" POTATO -- 2021-10-14 -- Version 1.1
    Developed by Lukáš Pekárek and Stefan Buck at the Helmholtz Institute for RNA-based Infection Research
    In the research group REMI - Recoding Mechanisms in Infections
    Supervisor - Jun. Prof. Neva Caliskan """
""" This script processes Force-Distance Optical Tweezers data in an automated way, to find unfolding events """
""" The script is developed to handle h5 raw data, produced from the C-Trap OT machine from Lumicks,
    as well as any other FD data prepared in a csv file (2 columns: Force(pN) - Distance(um)) """
""" Furthermore the script can analyse single constant force files """
""" The parameters can be changed in the GUI before each run.
    Alternatively they can be changed permanently in the POTATO_config file"""

import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from tkinter import filedialog
from tkinter import ttk
from PIL import ImageTk, Image
import pandas as pd
import numpy as np
import os
import h5py
import glob
import time
import multiprocessing as mp
import json

# relative imports
from POTATO_ForceRamp import start_subprocess, read_in_data
from POTATO_preprocessing import preprocess_RAW
from POTATO_config import default_values_HF, default_values_LF, default_values_CSV, default_values_FIT, default_values_constantF
from POTATO_constantF import get_constantF, display_constantF, fit_constantF
from POTATO_fitting import fitting_ds


# To avoid blurry GUI - DPI scaling
import ctypes
awareness = ctypes.c_int()
errorCode = ctypes.windll.shcore.GetProcessDpiAwareness(0, ctypes.byref(awareness))
errorCode = ctypes.windll.shcore.SetProcessDpiAwareness(1)


"""define the functions used in the GUI"""


# get settings, get folder directory, create analysis results folder
def start_analysis():
    global p0
    global analysis_folder

    # check user input
    input_settings, input_format, export_data, input_fitting, TOMATO_fitting, input_constantF = check_settings()

    # ask wich directory should be analysed
    folder = filedialog.askdirectory()
    root.title('POTATO -- ' + str(folder))

    # decide which input format was choosen
    if input_format['CSV'] == 1:
        folder_path = str(folder + "/*.csv")
    else:
        folder_path = str(folder + "/*.h5")

    Files = glob.glob(folder_path)

    # print number of files to analyse, if no files found give an error
    print('Files to analyse', len(Files))
    output_window.insert("end", 'Files to analyse: ' + str(len(Files)) + "\n")
    output_window.see("end")
    if not len(Files) == 0:
        output_window.insert("end", 'Analysis in progress. Please do not close the program! \n')

        # print starting time of the analysis
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        print("Timestamp: " + timestamp)
        output_window.insert("end", 'Start of analysis: ' + str(timestamp) + "\n")
        output_window.see("end")

        # create a folder for the analysis results
        analysis_folder = str(folder + '/Analysis_' + timestamp)
        os.mkdir(analysis_folder)

        # export configuration file with used parameters
        export_settings(analysis_folder, timestamp, input_settings, input_fitting)

        # start analysis in a new process
        p0 = mp.Process(target=start_subprocess, name='Process-0', args=(
            analysis_folder,
            timestamp,
            Files,
            input_settings,
            input_format,
            export_data,
            input_fitting,
            output_q,
        ))

        p0.daemon = True
        p0.start()

    else:
        output_window.insert("end", 'No file of the selected data type in the folder! \n')
        output_window.see("end")


# display default values in the GUI
def parameters(frame, default_values, default_fit, default_constantF):
    downsample_value1.delete(0, "end")
    downsample_value1.insert("end", default_values['Downsampling rate'])
    downsample_value2.delete(0, "end")
    downsample_value2.insert("end", default_values['Downsampling rate'])
    Filter_degree1.delete(0, "end")
    Filter_degree1.insert("end", default_values['Butterworth filter degree'])
    Filter_degree2.delete(0, "end")
    Filter_degree2.insert("end", default_values['Butterworth filter degree'])
    Filter_cut_off1.delete(0, "end")
    Filter_cut_off1.insert("end", default_values['Cut-off frequency'])
    Filter_cut_off2.delete(0, "end")
    Filter_cut_off2.insert("end", default_values['Cut-off frequency'])
    Force_Min1.delete(0, "end")
    Force_Min1.insert("end", default_values['Force threshold, pN'])
    Force_Min2.delete(0, "end")
    Force_Min2.insert("end", default_values['Force threshold, pN'])
    Z_score_force1.delete(0, "end")
    Z_score_force1.insert("end", default_values['Z-score force'])
    Z_score_force2.delete(0, "end")
    Z_score_force2.insert("end", default_values['Z-score force'])
    Z_score_distance1.delete(0, "end")
    Z_score_distance1.insert("end", default_values['Z-score distance'])
    Z_score_distance2.delete(0, "end")
    Z_score_distance2.insert("end", default_values['Z-score distance'])
    step_d_value.delete(0, "end")
    step_d_value.insert("end", str(default_values['Step d']))
    window_size_value.delete(0, "end")
    window_size_value.insert("end", str(default_values['Moving median window size']))
    STD_difference_value.delete(0, "end")
    STD_difference_value.insert("end", str(default_values['STD difference threshold']))
    Frequency_value.delete(0, "end")
    Frequency_value.insert("end", str(default_values['Data frequency, Hz']))

    dsLp.delete(0, "end")
    dsLp.insert("end", str(default_fit['Persistance-Length ds, nm']))
    dsLp_up.delete(0, "end")
    dsLp_up.insert("end", str(default_fit['Persistance-Length ds, upper bound, nm']))
    dsLp_low.delete(0, "end")
    dsLp_low.insert("end", str(default_fit['Persistance-Length ds, lower bound, nm']))
    ssLp.delete(0, "end")
    ssLp.insert("end", str(default_fit['Persistance-Length ss, nm']))
    dsLc.delete(0, "end")
    dsLc.insert("end", str(default_fit['Contour-Length ds, nm']))
    ssLc.delete(0, "end")
    ssLc.insert("end", str(default_fit['Persistance-Length ss, nm']))
    stiff_ds.delete(0, "end")
    stiff_ds.insert("end", str(default_fit['Stiffness ds, pN']))
    stiff_ds_up.delete(0, "end")
    stiff_ds_up.insert("end", str(default_fit['Stiffness ds, upper bound, pN']))
    stiff_ds_low.delete(0, "end")
    stiff_ds_low.insert("end", str(default_fit['Stiffness ds, lower bound, pN']))
    stiff_ss.delete(0, "end")
    stiff_ss.insert("end", str(default_fit['Stiffness ss, pN']))
    f_off.delete(0, "end")
    f_off.insert("end", str(default_fit['Force offset, pN']))
    f_off_up.delete(0, "end")
    f_off_up.insert("end", str(default_fit['Force offset, upper bound, pN']))
    f_off_low.delete(0, "end")
    f_off_low.insert("end", str(default_fit['Force offset, lower bound, pN']))
    d_off.delete(0, "end")
    d_off.insert("end", str(default_fit['Distance offset, nm']))
    d_off_up.delete(0, "end")
    d_off_up.insert("end", str(default_fit['Distance offset, upper bound, nm']))
    d_off_low.delete(0, "end")
    d_off_low.insert("end", str(default_fit['Distance offset, lower bound, nm']))

    x_min.delete(0, "end")
    x_min.insert("end", str(default_constantF['x min']))
    x_max.delete(0, "end")
    x_max.insert("end", str(default_constantF['x max']))
    y_min.delete(0, "end")
    y_min.insert("end", str(default_constantF['y min']))
    y_max.delete(0, "end")
    y_max.insert("end", str(default_constantF['y max']))
    number_gauss.delete(0, "end")
    number_gauss.insert("end", str(default_constantF['Number gauss']))
    mean_gauss.delete(0, "end")
    mean_gauss.insert("end", str(default_constantF['Mean']))
    STD_gauss.delete(0, "end")
    STD_gauss.insert("end", str(default_constantF['STD']))
    amplitude_gauss.delete(0, "end")
    amplitude_gauss.insert("end", str(default_constantF['Amplitude']))


# update value that has been changed (needs event <ENTER>)
def user_input(event, param1, param2):
    new_param = param1.get()
    param1.delete(0, "end")
    param2.delete(0, "end")
    param1.insert("end", new_param)
    param2.insert("end", new_param)


# get all settings from the user input before start of the analysis
def check_settings():
    input_settings = {
        'downsample_value': int(downsample_value2.get()),
        'filter_degree': int(Filter_degree2.get()),
        'filter_cut_off': float(Filter_cut_off2.get()),
        'F_min': float(Force_Min2.get()),
        'step_d': int(step_d_value.get()),
        'z-score_f': float(Z_score_force2.get()),
        'z-score_d': float(Z_score_distance2.get()),
        'window_size': int(window_size_value.get()),
        'data_frequency': float(Frequency_value.get()),
        'STD_diff': float(STD_difference_value.get())
    }

    input_format = {
        'HF': check_box_HF.get(),
        'LF': check_box_LF.get(),
        'CSV': check_box_CSV.get(),
        'Trap': check_box_Trap1.get(),
        'length_measure': check_box_um.get(),
        'MultiH5': check_box_multiH5.get(),
        'preprocess': check_box_preprocess.get()
    }

    export_data = {
        'export_SMOOTH': check_box_smooth_data.get(),
        'export_PLOT': check_box_plot.get(),
        'export_STEPS': check_box_steps.get(),
        'export_TOTAL': check_box_total_results.get(),
        'export_FIT': check_box_fitting.get()
    }

    input_fitting = {
        'WLC+WLC': int(check_box_WLC.get()),
        'WLC+FJC': int(check_box_FJC.get()),
        'lp_ds': float(dsLp.get()),
        'lp_ds_up': float(dsLp_up.get()),
        'lp_ds_low': float(dsLp_low.get()),
        'lc_ds': float(dsLc.get()),
        'lp_ss': float(ssLp.get()),
        'lc_ss': float(ssLc.get()),
        'ds_stiff': float(stiff_ds.get()),
        'ds_stiff_up': float(stiff_ds_up.get()),
        'ds_stiff_low': float(stiff_ds_low.get()),
        'ss_stiff': float(stiff_ss.get()),
        'offset_f': float(f_off.get()),
        'offset_f_up': float(f_off_up.get()),
        'offset_f_low': float(f_off_low.get()),
        'offset_d': float(d_off.get()),
        'offset_d_up': float(d_off_up.get()),
        'offset_d_low': float(d_off_low.get())
    }

    TOMATO_fitting = {
        'WLC+WLC': int(check_box_WLC.get()),
        'WLC+FJC': int(check_box_FJC.get()),
        'lp_ds': float(entryText_ds_Lp.get()),
        'lp_ds_up': float(dsLp_up.get()),
        'lp_ds_low': float(dsLp_low.get()),
        'lc_ds': float(entryText_ds_Lc.get()),
        'lp_ss': float(entryText_ss_Lp.get()),
        'lc_ss': float(entryText_ss_Lc.get()),
        'ss_stiff': float(entryText_ss_St.get()),
        'offset_f': float(entryText_shift_F.get()),
        'offset_f_up': float(f_off_up.get()),
        'offset_f_low': float(f_off_low.get()),
        'offset_d': float(entryText_shift_d.get()),
        'offset_d_up': float(d_off_up.get()),
        'offset_d_low': float(d_off_low.get()),
        'ds_stiff': float(entryText_ds_St.get()),
        'ds_stiff_up': float(entryText_ds_St.get()) + float(stiff_ds_up.get()),
        'ds_stiff_low': float(entryText_ds_St.get()) - float(stiff_ds_low.get())
    }

    input_constantF = {
        'x min': int(x_min.get()),
        'x max': int(x_max.get()),
        'y min': int(y_min.get()),
        'y max': int(y_max.get()),
        'Number gauss': int(number_gauss.get()),
        'Mean': mean_gauss.get(),
        'STD': STD_gauss.get(),
        'Amplitude': amplitude_gauss.get()
    }

    return input_settings, input_format, export_data, input_fitting, TOMATO_fitting, input_constantF


def export_settings(analysis_path, timestamp, input_1, input_2):
    with open(str(analysis_path + '/parameters_' + timestamp + '.txt'), 'w') as config_used:
        config_used.write('Data processing:\n')
        config_used.write(json.dumps(input_1, indent=4, sort_keys=False))
        config_used.write('\n\n')
        config_used.write('Fitting parameters:\n')
        config_used.write(json.dumps(input_2, indent=4, sort_keys=False))


# Looks for output of the subprocess
def refresh():
    global new_image
    while output_q.empty() is False:
        output = output_q.get()
        output_window.insert("end", "\n" + output + "\n")
        output_window.see("end")
    try:
        images = str(analysis_folder + "/*plot*.png")
        list_images = glob.glob(images)
        img = Image.open(list_images[-1])
        resized = img.resize((1000, 650), Image.ANTIALIAS)
        new_image = ImageTk.PhotoImage(resized)
        figure_frame.create_image(0, 0, image=new_image, anchor="nw")
    except:
        pass


def readme():
    with open("POTATO_readme.txt", "r") as f:
        help_text = f.read()
        help_window = tk.Toplevel(root)
        help_window.title("Readme")
        text = tk.Text(help_window, height=25, width=200)
        text.grid(row=0, column=0, sticky="nw")
        text.insert("end", help_text)


# display a single h5 file (tab2)
def getRAW_File_h5():
    input_settings, input_format, export_data, input_fitting, TOMATO_fitting, input_constantF = check_settings()
    import_file_path = filedialog.askopenfilename()
    with h5py.File(import_file_path, "r") as raw_data:
        # access the raw data
        if input_format['HF'] == 1:
            if input_format['Trap'] == 1:
                Force = raw_data.get("Force HF/Force 1x")
            elif input_format['Trap'] == 0:
                Force = raw_data.get("Force HF/Force 2x")
            Distance = raw_data.get("Distance/Piezo Distance")

        elif input_format['LF'] == 1:
            if input_format['Trap'] == 1:
                load_force = raw_data.get("Force LF/Force 1x")
                Force = load_force[:]['Value'][:]
                load_distance = raw_data.get("Distance/Distance 1x")[:]
                Distance = load_distance['Value'][:]
            elif input_format['Trap'] == 0:
                load_force = raw_data.get("Force LF/Force 2x")
                Force = load_force[:]['Value'][:]
                load_distance = raw_data.get("Distance/Distance 2x")[:]
                Distance = load_distance['Value'][:]
        if input_format['preprocess'] == 1:
            FD, FD_um = preprocess_RAW(Force, Distance, input_settings)
            display_RAW_FD(FD[:, 0], FD[:, 1], Force[::input_settings['downsample_value']], Distance[::input_settings['downsample_value']] * 1000)
        else:
            Force = np.array(Force)
            Distance = np.array(Distance)
            display_RAW_FD(Force, Distance, Force, Distance)


# display a single csv file (tab2)
def getRAW_File_csv():
    if not check_box_CSV.get() == 1:
        check_box_CSV.set(value=1)
        select_box(check_box_CSV, check_box_HF, check_box_LF)
        parameters(parameter_frame, default_values_CSV, default_values_FIT, default_values_constantF)
        parameters(tab3, default_values_CSV, default_values_FIT, default_values_constantF)
    else:
        pass

    input_settings, input_format, export_data, input_fitting, TOMATO_fitting, input_constantF = check_settings()
    import_file_path = filedialog.askopenfilename()

    df = pd.read_csv(import_file_path)

    Force = df.to_numpy()[:, 0]
    Distance = df.to_numpy()[:, 1]

    if input_format['preprocess'] == 1:
        FD, FD_um = preprocess_RAW(Force, Distance, input_settings)
        display_RAW_FD(FD[:, 0], FD[:, 1], Force[::input_settings['downsample_value']], Distance[::input_settings['downsample_value']] * 1000)
    else:
        display_RAW_FD(Force, Distance, Force, Distance)


# create the plot for tab2
def display_RAW_FD(processed_F, processed_D, raw_F, raw_D):
    single_fd = Figure(figsize=(10, 6), dpi=100)
    subplot1 = single_fd.add_subplot(111)

    legend_elements = [
        Line2D([0], [0], color='C0', lw=4),
        Line2D([0], [0], color='C1', lw=4)
    ]

    subplot1.set_xlabel("Distance (nm)")
    subplot1.set_ylabel("Force (pN)")
    subplot1.plot(raw_D, raw_F, alpha=0.8, color='C0', zorder=0)
    subplot1.scatter(processed_D, processed_F, marker='.', s=0.1, linewidths=None, alpha=1, color='C1', zorder=1)
    subplot1.legend(legend_elements, ['Downsampled FD-Data', 'Filtered FD-Data'])

    figure_raw = FigureCanvasTkAgg(single_fd, figure_frame2)
    figure_raw.get_tk_widget().grid(row=0, column=0, sticky='wens')

    tabControl.select(tab2)


def start_constantF():
    input_settings, input_format, export_data, input_fitting, TOMATO_fitting, input_constantF = check_settings()
    Force_Distance, Force_Distance_um, frequency, filename, analysis_path, timestamp = get_constantF(input_settings, input_format, input_constantF)
    fig_constantF, hist_D, filteredDistance_ready = display_constantF(Force_Distance, Force_Distance_um, frequency, input_settings, input_constantF)
    os.mkdir(analysis_path)
    export_settings(analysis_path, timestamp, input_settings, input_constantF)
    fig_constantF = fit_constantF(hist_D, Force_Distance, filteredDistance_ready, frequency, input_settings, input_constantF, filename, timestamp)
    fig_constantF_tk = FigureCanvasTkAgg(fig_constantF, figure_frame_tab4)
    fig_constantF_tk.get_tk_widget().grid(row=0, column=0, sticky='wens')

    tabControl.select(tab4)


def show_constantF():
    input_settings, input_format, export_data, input_fitting, TOMATO_fitting, input_constantF = check_settings()
    Force_Distance, Force_Distance_um, frequency, filename, analysis_path, timestamp = get_constantF(input_settings, input_format, input_constantF)
    fig_constantF, hist_D, filteredDistance_ready = display_constantF(Force_Distance, Force_Distance_um, frequency, input_settings, input_constantF)
    fig_constantF_tk = FigureCanvasTkAgg(fig_constantF, figure_frame_tab4)
    fig_constantF_tk.get_tk_widget().grid(row=0, column=0, sticky='wens')

    tabControl.select(tab4)


def on_closing():
    # makes sure all python processes/loops are cancelled before exiting
    if tk.messagebox.askokcancel("Quit", "Do you really want to quit?"):
        root.quit()


################ TOMATO ###############################
# from POTATO_TOMATO import open_folder, create_chart, clear_charts, start_click, \
#     end_click, Fitting_WLC_ds_handles, Fitting_WLC_ss_handles, export_table, clear_table_last, \
#     clear_table, reset_parameters, start_work_click, end_work_click, calc_rWork, calc_strWork, next_FD, \
#     previous_FD_key, next_FD_key, save_key, start_click_key, end_click_key, \
#     end_work_click_key, zero_str_work_key, fit_ds_key, fit_ss_key, calc_rWork_key, calc_strWork_key, start_work_click_key, \
#     load_previous_data_key,  previous_FD, write_to_table, export_model, zero_str_work
import lumicks.pylake as lk
from scipy.integrate import simps
from POTATO_TOMATO import plot_TOMATO


############# define the functions for TOMATO ##################
def open_folder():
    global filename_TOMATO
    global Force_Distance_TOMATO
    global import_file_path
    global TOMATO_fig1
    global Files
    global FD_number
    # check user input
    input_settings, input_format, export_data, input_fitting, TOMATO_fitting, input_constantF = check_settings()

    # ask wich directory should be analysed
    folder = filedialog.askdirectory()
    root.title('POTATO -- ' + str(folder))

    # decide which input format was choosen
    if input_format['CSV'] == 1:
        folder_path = str(folder + "/*.csv")
    else:
        folder_path = str(folder + "/*.h5")

    Files = glob.glob(folder_path)

    FD_number = 0
    Force_Distance_TOMATO, Force_Distance_um_TOMATO, Frequency_value, filename_TOMATO = read_in_data(FD_number, Files, input_settings, input_format)
    entryText_filename.set(filename_TOMATO)

    reset_parameters()
    fig = plot_TOMATO(Force_Distance_TOMATO)
    TOMATO_fig1 = FigureCanvasTkAgg(fig, TOMATO_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0, sticky='wens')
    toolbarFrame = tk.Frame(master=TOMATO_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)


def change_FD(direction):
    global TOMATO_fig1
    global filename_TOMATO
    global FD_number
    global Force_Distance_TOMATO
    FD_number = FD_number + direction

    input_settings, input_format, export_data, input_fitting, TOMATO_fitting, input_constantF = check_settings()
    Force_Distance_TOMATO, Force_Distance_um_TOMATO, Frequency_value, filename_TOMATO = read_in_data(FD_number, Files, input_settings, input_format)

    if Force_Distance_TOMATO[0, 1] > Force_Distance_TOMATO[-1, 1]:  # reverse
        Force_Distance_TOMATO = np.flipud(Force_Distance_TOMATO)
        Force_Distance_um_TOMATO = np.flipud(Force_Distance_um_TOMATO)

    entryText_filename.set(filename_TOMATO)

    reset_parameters()
    fig = plot_TOMATO(Force_Distance_TOMATO)
    TOMATO_fig1 = FigureCanvasTkAgg(fig, TOMATO_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0, sticky='wens')
    toolbarFrame = tk.Frame(master=TOMATO_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)


def save_previous_data():
    global TOMATO_dict
    TOMATO_dict = {'shift_d': entryText_shift_d.get(),
                   'shift_F': entryText_shift_F.get(),
                   'dsWork': entryText_dsWork.get(),
                   'ds_St': entryText_ds_St.get(),
                   'ss_St': entryText_ss_St.get(),
                   'ssWork': entryText_ssWork.get(),
                   'ds_Lp': entryText_ds_Lp.get(),
                   'ss_Lp': entryText_ss_Lp.get(),
                   'rWork': entryText_rWork.get(),
                   'ds_Lc': entryText_ds_Lc.get(),
                   'ss_Lc': entryText_ss_Lc.get(),
                   'strWork': entryText_strWork.get(),
                   'fit_end': entryText_end.get(),
                   'fit_start': entryText_start.get(),
                   'start_work_D': entryText_start_work_D.get(),
                   'start_work_F': entryText_start_work_F.get(),
                   'end_work_D': entryText_end_work_D.get(),
                   'end_work_F': entryText_end_work_F.get(),
                   }


    try:
        TOMATO_dict['FD plot'] = Force_Distance_TOMATO
    except:
        pass

    try:
        TOMATO_dict['ds_fitted_region'] = [Force_Distance_TOMATO[:, 1][real_start_2: real_end_2], Force_Distance_TOMATO[:, 0][real_start_2: real_end_2]]
    except:
        pass

    try:
        TOMATO_dict['ss_fitted_region'] = [Force_Distance_TOMATO[:, 1][real_start_3:real_end_3], Force_Distance_TOMATO[:, 0][real_start_3:real_end_3]]
    except:
        pass

    try:
        TOMATO_dict['ds_fit_plot'] = [distance, F_ds_model]
    except:
        pass

    try:
        TOMATO_dict['ss_fit_plot'] = [distance, F_ss_model]
    except:
        pass

    try:
        TOMATO_dict['result_table'] = listBox
    except:
        pass


def load_previous_data():
    entryText_shift_d.set(TOMATO_dict['shift_d'])
    entryText_shift_F.set(TOMATO_dict['shift_F'])
    entryText_dsWork.set(TOMATO_dict['dsWork'])
    entryText_ds_St.set(TOMATO_dict['ds_St'])
    entryText_ss_St.set(TOMATO_dict['ss_St'])
    entryText_ssWork.set(TOMATO_dict['ssWork'])
    entryText_ds_Lp.set(TOMATO_dict['ds_Lp'])
    entryText_ss_Lp.set(TOMATO_dict['ss_Lp'])
    entryText_rWork.set(TOMATO_dict['rWork'])
    entryText_ds_Lc.set(TOMATO_dict['ds_Lc'])
    entryText_ss_Lc.set(TOMATO_dict['ss_Lc'])
    entryText_strWork.set(TOMATO_dict['strWork'])
    entryText_end.set(TOMATO_dict['fit_end'])
    entryText_start.set(TOMATO_dict['fit_start'])
    entryText_start_work_D.set(TOMATO_dict['start_work_D'])
    entryText_start_work_F.set(TOMATO_dict['start_work_F'])
    entryText_end_work_D.set(TOMATO_dict['end_work_D'])
    entryText_end_work_F.set(TOMATO_dict['end_work_F'])

    global TOMATO_fig1
    fig = plot_TOMATO(np.array([TOMATO_dict['FD plot'][1], TOMATO_dict['FD plot'][0]]))

    TOMATO_fig1 = FigureCanvasTkAgg(fig, TOMATO_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0, sticky='wens')

    toolbarFrame = tk.Frame(master=TOMATO_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)

    #  ds fitted region
    try:
        subplot1.plot(TOMATO_dict['ds_fitted_region'][0], TOMATO_dict['ds_fitted_region'][1], color="b")
    except:
        pass

    # ss fitted region
    try:
        subplot1.plot(TOMATO_dict['ss_fitted_region'][0], TOMATO_dict['ss_fitted_region'][1], color="r")
    except:
        pass

    # ds fit plot
    try:
        subplot1.plot(TOMATO_dict['ds_fit_plot'][0], TOMATO_dict['ds_fit_plot'][1], marker=None,linestyle='dashed',linewidth=1,color="black")
    except:
        pass

    # ss fit plot
    try:
        subplot1.plot(TOMATO_dict['ss_fit_plot'][0], TOMATO_dict['ss_fit_plot'][1], marker=None, linestyle='dashed', linewidth=1, color="black")
    except:
        pass


# key binding wrapper functions
def previous_FD_key(event):
    change_FD(-1)


def next_FD_key(event):
    change_FD(+1)


def save_key(event):
    write_to_table()


def start_click_key(event):
    start_click()


def end_click_key(event):
    end_click()


def start_work_click_key(event):
    start_work_click()


def end_work_click_key(event):
    end_work_click()


def zero_str_work_key(event):
    zero_str_work()


def fit_ds_key(event):
    Fitting_WLC_ds_handles()


def fit_ss_key(event):
    Fitting_WLC_ss_handles()


def load_previous_data_key(event):
    load_previous_data()


def calc_rWork_key(event):
    calc_rWork()


def calc_strWork_key(event):
    calc_strWork()


def create_chart():
    global TOMATO_fig1

    fig = plot_TOMATO(Force_Distance_TOMATO)
    TOMATO_fig1 = FigureCanvasTkAgg(fig, TOMATO_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0, sticky='wens')
    toolbarFrame = tk.Frame(master=TOMATO_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)


def start_click():
    global cid
    cid = TOMATO_fig1.mpl_connect('button_press_event', lambda event, arg=1: onclick_start_end(event, arg))


def end_click():
    global cid
    cid = TOMATO_fig1.mpl_connect('button_press_event', lambda event, arg=0: onclick_start_end(event, arg))


def onclick_start_end(event, pos):
    global cid

    PD_position, F_position = float(event.xdata), float(event.ydata)
    print(PD_position, F_position)
    if pos == 1:
        entryText_start.set(round(PD_position, 1))
    elif pos == 0:
        entryText_end.set(round(PD_position, 1))
    TOMATO_fig1.mpl_disconnect(cid)


def start_work_click():
    global cid
    cid = TOMATO_fig1.mpl_connect('button_press_event', lambda event, arg=1: onclick_work(event, arg))


def end_work_click():
    global cid
    cid = TOMATO_fig1.mpl_connect('button_press_event', lambda event, arg=0: onclick_work(event, arg))


def onclick_work(event, pos):
    global cid

    PD_position, F_position = float(event.xdata), float(event.ydata)
    print(PD_position, F_position)
    if pos == 1:
        entryText_start_work_D.set(round(PD_position, 1))
        entryText_start_work_F.set(round(F_position, 3))
    elif pos == 0:
        entryText_end_work_D.set(round(PD_position, 1))
        entryText_end_work_F.set(round(F_position, 3))
    TOMATO_fig1.mpl_disconnect(cid)


def calc_rWork():
    x1 = float(entryText_start_work_D.get())
    x2 = float(entryText_end_work_D.get())

    y1 = float(entryText_start_work_F.get())
    y2 = float(entryText_end_work_F.get())

    rWork = (x2 - x1) * (y1 + y2) / 2 / 4.114

    entryText_rWork.set(rWork)
    entryText_strWork.set("0")


def calc_strWork():
    dsWork = float(entryText_dsWork.get())
    rWork = float(entryText_rWork.get())
    ssWork = float(entryText_ssWork.get())

    strWork = dsWork + rWork - ssWork
    entryText_strWork.set(strWork)


def write_to_table():
    global listBox

    if float(entryText_strWork.get()) == 0:
        work_done = entryText_rWork.get()
    else:
        work_done = entryText_strWork.get()
    listBox.insert("", "end", values=(filename_TOMATO, entryText_start_work_F.get(), entryText_end_work_F.get(), (float(entryText_start_work_F.get()) + float(entryText_end_work_F.get())) / 2, entryText_start_work_D.get(), entryText_end_work_D.get(), float(entryText_end_work_D.get()) - float(entryText_start_work_D.get()), entry_ds_Lc.get(), entry_ds_Lp.get(), entry_ds_St.get(), entry_ss_Lc.get(), entry_ss_Lp.get(), entry_ss_St.get(), entry_shift_F.get(), entry_shift_d.get(), work_done))


def clear_charts():
    TOMATO_fig1.get_tk_widget().grid_forget()


def clear_table():
    global listBox
    list_items = listBox.get_children("")

    for item in list_items:
        listBox.delete(item)


def clear_table_last():
    global listBox
    list_items = listBox.get_children("")

    listBox.delete(list_items[-1])


def reset_parameters():
    entryText_shift_d.set("0")
    entryText_shift_F.set("0")
    entryText_ds_Lp.set("40")
    entryText_ds_Lc.set("1256")
    entryText_ss_Lc.set("0")
    entryText_ss_Lp.set("1")
    entryText_ss_St.set("800")
    entryText_ds_St.set("400")
    entryText_dsWork.set("0")
    entryText_ssWork.set("0")
    entryText_rWork.set("0")
    entryText_strWork.set("0")


def fitting_ss_TOMATO(PD_ss, F_ss, Ds_fit_dict, fix, max_range):

    model_ss = lk.odijk("DNA_2") + lk.odijk("RNA")

    model_ss = model_ss.invert().subtract_independent_offset() + lk.force_offset("DNA")
    fit_ss = lk.FdFit(model_ss)

    fit_ss.add_data("ss_part", F_ss, PD_ss)

    ## ds part parameters

    # Persistance length bounds
    # Lp_ds_range=fit_ds["DNA/Lp"].value/10
    fit_ss["DNA_2/Lp"].value = Ds_fit_dict['Lp_ds']
    fit_ss["DNA_2/Lp"].lower_bound = Ds_fit_dict['Lp_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/Lp"].upper_bound = Ds_fit_dict['Lp_ds'] * (1 + max_range / 100)
    # if fix==1:
    fit_ss["DNA_2/Lp"].fixed = 'True'
    fit_ss["DNA/f_offset"].upper_bound = float(f_off_up.get())
    fit_ss["DNA/f_offset"].lower_bound = float(f_off_low.get())
    fit_ss["DNA/f_offset"].value = Ds_fit_dict['f_offset']
    fit_ss["DNA/f_offset"].fixed = 'True'

    fit_ss["inv(DNA_2_with_RNA)/d_offset"].value = Ds_fit_dict['d_offset']
    fit_ss["inv(DNA_2_with_RNA)/d_offset"].fixed = 'True'

    # contour length
    # Lc_ds_range=Lc_initial_guess/100 # nm
    fit_ss["DNA_2/Lc"].upper_bound = Ds_fit_dict['Lc_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lc"].lower_bound = Ds_fit_dict['Lc_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/Lc"].value = Ds_fit_dict['Lc_ds']
    fit_ss["DNA_2/Lc"].unit = 'nm'
    # if fix==1:
    fit_ss["DNA_2/Lc"].fixed = 'True'

    # stifness

    fit_ss["DNA_2/St"].upper_bound = Ds_fit_dict['St_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/St"].lower_bound = Ds_fit_dict['St_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/St"].value = Ds_fit_dict['St_ds']
    if fix == 1:
        fit_ss["DNA_2/St"].fixed = 'True'

    ## ss part parameters
    # Persistance length bounds

    fit_ss["RNA/Lp"].value = float(entryText_ss_Lp.get())
    fit_ss["RNA/Lp"].lower_bound = 0.8
    fit_ss["RNA/Lp"].upper_bound = 2
    if fix == 1:
        fit_ss["RNA/Lp"].fixed = 'True'

    # stiffnes
    fit_ss["RNA/St"].value = float(entryText_ss_St.get())
    fit_ss["RNA/St"].lower_bound = 300
    fit_ss["RNA/St"].upper_bound = 1500
    # contour length

    fit_ss["RNA/Lc"].upper_bound = float(entryText_ss_Lc.get()) + 100
    fit_ss["RNA/Lc"].lower_bound = 0
    fit_ss["RNA/Lc"].value = float(entryText_ss_Lc.get())
    fit_ss["RNA/Lc"].unit = 'nm'

    fit_ss.fit()

    Fit_dict = {'model': model_ss, 'fit_model': fit_ss, 'Lc_ds': fit_ss["DNA_2/Lc"].value, 'Lp_ds': fit_ss["DNA_2/Lp"].value, 'St_ds': fit_ss["DNA_2/St"].value, 'Lc_ss': fit_ss["RNA/Lc"].value, 'Lp_ss': fit_ss["RNA/Lp"].value, 'St_ss': fit_ss["RNA/St"].value, 'f_offset': fit_ss["DNA/f_offset"].value, 'd_offset': fit_ss["inv(DNA_2_with_RNA)/d_offset"].value}
    return Fit_dict


def Fitting_WLC_ds_handles():
    # create a sublist of the ROI PD_nm
    global ds_fit_dict_TOMATO
    global F_region
    global F_ds_model
    global distance
    global real_start, real_end
    global real_start_2, real_end_2
    # find match with PD
    save_previous_data()
    input_settings, input_format, export_data, input_fitting, TOMATO_fitting, input_constantF = check_settings()

    real_PD = []
    start_PD = float(entry_start.get())
    end_PD = float(entry_end.get())

    for i in [start_PD, end_PD]:
        absolute_difference_function = lambda cPD: abs(cPD - i)
        real_PD.append(min(Force_Distance_TOMATO[:, 1], key=absolute_difference_function))

    PD_nm_list = list(Force_Distance_TOMATO[:, 1])

    real_start = PD_nm_list.index(real_PD[0])
    real_end = PD_nm_list.index(real_PD[1])

    PD_region = []
    F_region = []
    if real_start < real_end:
        for i in range(real_start, real_end, 10):
            PD_region.append(Force_Distance_TOMATO[:, 1][i])
            F_region.append(Force_Distance_TOMATO[:, 0][i])

    else:
        for i in range(real_end, real_start, 10):
            PD_region.append(Force_Distance_TOMATO[:, 1][i])
            F_region.append(Force_Distance_TOMATO[:, 0][i])

    Force_Distance_ds_fit = np.array([F_region, PD_region])
    ds_fit_dict_TOMATO, area_ds_TOMATO = fitting_ds(filename_TOMATO, input_settings, export_data, TOMATO_fitting, real_end, Force_Distance_ds_fit, None, None, 1)

    entryText_ds_Lp.set(ds_fit_dict_TOMATO['Lp_ds'])
    entryText_shift_F.set(ds_fit_dict_TOMATO['f_offset'])
    entryText_shift_d.set(ds_fit_dict_TOMATO["d_offset"])
    entryText_ds_Lc.set(ds_fit_dict_TOMATO['Lc_ds'])
    entryText_ds_St.set(ds_fit_dict_TOMATO['St_ds'])

    # plot the marked region and fitted WLC
    global TOMATO_fig1
    global figure1
    global subplot1
    # model data
    distance = np.arange(min(Force_Distance_TOMATO[:, 1]), max(Force_Distance_TOMATO[:, 1]) + 50, 2)
    F_ds_model = ds_fit_dict_TOMATO['model_ds'](distance, ds_fit_dict_TOMATO['fit_model'])

    figure1 = plot_TOMATO(Force_Distance_TOMATO)

    if real_start < real_end:
        real_start_2 = real_start
        real_end_2 = real_end
    else:
        real_start_2 = real_end
        real_end_2 = real_start

    subplot1 = figure1.add_subplot(111)
    subplot1.plot(Force_Distance_TOMATO[:, 1][real_start_2: real_end_2], Force_Distance_TOMATO[:, 0][real_start_2:real_end_2], color="b")
    subplot1.plot(distance, F_ds_model, marker=None, linestyle='dashed', linewidth=1, color="black")
    subplot1.set_ylim([min(Force_Distance_TOMATO[:, 0]), max(Force_Distance_TOMATO[:, 0])])
    subplot1.set_xlim([min(Force_Distance_TOMATO[:, 1]) - 10, max(Force_Distance_TOMATO[:, 1]) + 10])
    subplot1.tick_params('both', direction='in')

    TOMATO_fig1 = FigureCanvasTkAgg(figure1, TOMATO_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0)

    toolbarFrame = tk.Frame(master=TOMATO_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)

    entryText_dsWork.set(area_ds_TOMATO)
    print("area_ds = " + str(area_ds_TOMATO))
    # add the parameters to table


## fitting the ss RNA part combined with ds handles part
def Fitting_WLC_ss_handles():
    # create a sublist of the ROI PD_nm
    global F_region
    global F_ss_model
    global distance
    global real_start_3, real_end_3
    # find match with PD

    save_previous_data()

    real_PD = []
    start_PD = float(entry_start.get())
    end_PD = float(entry_end.get())

    for i in [start_PD, end_PD]:
        absolute_difference_function = lambda cPD: abs(cPD - i)
        real_PD.append(min(Force_Distance_TOMATO[:, 1], key=absolute_difference_function))
    # print(real_PD)

    PD_nm_list = list(Force_Distance_TOMATO[:, 1])
    real_start = PD_nm_list.index(real_PD[0])
    real_end = PD_nm_list.index(real_PD[1])

    PD_region = []
    F_region = []
    if abs(real_start - real_end) > 1000:
        if real_start < real_end:
            for i in range(real_start, real_end, 1000):
                PD_region.append(Force_Distance_TOMATO[:, 1][i])
                F_region.append(Force_Distance_TOMATO[:, 0][i])

        else:
            for i in range(real_end, real_start, 1000):
                PD_region.append(Force_Distance_TOMATO[:, 1][i])
                F_region.append(Force_Distance_TOMATO[:, 0][i])
    else:
        if real_start < real_end:
            for i in range(real_start, real_end, 100):
                PD_region.append(Force_Distance_TOMATO[:, 1][i])
                F_region.append(Force_Distance_TOMATO[:, 0][i])

        else:
            for i in range(real_end, real_start, 100):
                PD_region.append(Force_Distance_TOMATO[:, 1][i])
                F_region.append(Force_Distance_TOMATO[:, 0][i])

    #fitting itself
    Fit_ss = fitting_ss_TOMATO(PD_region, F_region, ds_fit_dict_TOMATO, 1, 1)
    entryText_ds_Lp.set(Fit_ss['Lp_ds'])
    entryText_shift_F.set(Fit_ss['f_offset'])
    entryText_shift_d.set(Fit_ss["d_offset"])
    entryText_ds_Lc.set(Fit_ss['Lc_ds'])
    entryText_ss_Lc.set(Fit_ss['Lc_ss'])
    entryText_ss_Lp.set(Fit_ss['Lp_ss'])
    entryText_ss_St.set(Fit_ss['St_ss'])
    entryText_ds_St.set(Fit_ss['St_ds'])

    # model data
    # distance = np.arange(min(PD_nm), max(PD_nm), 1)
    F_ss_model = Fit_ss['model'](distance, Fit_ss['fit_model'])
    # plot the marked region and fitted WLC
    global TOMATO_fig1
    global figure1
    global subplot1

    if real_start < real_end:
        real_start_3 = real_start
        real_end_3 = real_end
    else:
        real_start_3 = real_end
        real_end_3 = real_start

    subplot1.plot(Force_Distance_TOMATO[:, 1][real_start_3:real_end_3], Force_Distance_TOMATO[:, 0][real_start_3:real_end_3], color="r")
    subplot1.plot(distance, F_ss_model, marker=None, linewidth=1, linestyle='dashed', color="black")

    subplot1.set_ylim([min(Force_Distance_TOMATO[:, 0]), max(Force_Distance_TOMATO[:, 0])])
    subplot1.set_xlim([min(Force_Distance_TOMATO[:, 1]) - 10, max(Force_Distance_TOMATO[:, 1]) + 10])

    TOMATO_fig1 = FigureCanvasTkAgg(figure1, TOMATO_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0)

    toolbarFrame = tk.Frame(master=TOMATO_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)

    distance_integral = np.arange(float(entryText_start_work_D.get()), float(entryText_end_work_D.get()))
    ss_integral = Fit_ss['model'](distance_integral, Fit_ss['fit_model'])
    area_ss = simps(ss_integral) / 4.114
    entryText_ssWork.set(area_ss)
    # add the parameters to table


def export_table():
    global listBox
    global name
    global Fit_results
    ''' exporting the table results '''
    results = []
    for child in listBox.get_children():
        results.append(listBox.item(child)['values'])

    Fit_results = pd.DataFrame(results,
                            columns=[
                                'Filename',
                                'F1',
                                'F2',
                                'F1/2',
                                'Step start',
                                'Step end',
                                'Step length',
                                'ds Contour length',
                                'ds Persistance Length',
                                'ds St',
                                'ss Contour Length',
                                'ss Persistance Length',
                                'ss St',
                                'Shift F',
                                'Shift x',
                                'Work'
                            ]
    )

    name = filedialog.asksaveasfile(mode='w', defaultextension=".csv")
    print(name)
    Fit_results.to_csv(name.name, index=False, header=True)


def export_model():
    global listBox
    global name
    global Fit_results
    ''' exporting ds and ss model '''
    try:
        F_ss_model
        model_data = pd.DataFrame(list(zip(distance, F_ds_model, F_ss_model)), columns=['Distance [nm]', 'Force WLC data [pN]', 'Force WLC+FJC data [pN]'])
    except NameError:
        model_data = pd.DataFrame(list(zip(distance, F_ds_model)), columns=['Distance [nm]', 'Force WLC data [pN]'])

    name = filedialog.asksaveasfile(mode='w', defaultextension=".csv")
    name_model = name.name[:-4] + '_model_data.csv'
    model_data.to_csv(name_model, index=False, header=True)

    ''' exporting figure '''
    plotname = name.name[:-4] + '_graph.png'
    figure1.savefig(plotname, dpi=600)


def zero_str_work():
    entryText_strWork.set("0")
############## TOMATO functions end ###################


""" start the main process and Tkinter application """
if __name__ == '__main__':
    mp.freeze_support()
    root = tk.Tk()
    root.iconbitmap('POTATO.ico')
    root.title("POTATO -- Practical Optical Tweezers Analysis TOol")

    output_q = mp.Queue()

    # create a drop down menu
    drop_down_menu = tk.Menu(root)
    root.config(menu=drop_down_menu)

    # first drop down possibility: File
    file_menu = tk.Menu(drop_down_menu, tearoff=0)
    drop_down_menu.add_cascade(label='File', menu=file_menu)
    file_menu.add_command(label='Analyse folder (FD curves)', command=start_analysis)
    file_menu.add_command(label='Display single FD curve (h5)', command=getRAW_File_h5)
    file_menu.add_command(label='Display single FD curve (csv)', command=getRAW_File_csv)
    file_menu.add_separator()
    file_menu.add_command(label='Display constant force', command=show_constantF)
    file_menu.add_command(label='Fit constant force', command=start_constantF)

    # second drop down possibility: Settings
    settings_menu = tk.Menu(drop_down_menu, tearoff=0)
    drop_down_menu.add_cascade(label='Settings', menu=settings_menu)
    settings_menu.add_command(label='Set advanced settings', command=lambda: tabControl.select(tab3))

    # third drop down possibility: Help
    help_menu = tk.Menu(drop_down_menu, tearoff=0)
    drop_down_menu.add_cascade(label='Help', menu=help_menu)
    help_menu.add_command(label='Readme', command=readme)

    # Create different GUI tabs
    tabControl = ttk.Notebook(root)
    tabControl.grid(row=0, column=0, padx=2, pady=2)

    tab1 = ttk.Frame(tabControl, width=800, height=600)
    tab2 = ttk.Frame(tabControl, width=800, height=600)
    tab3 = ttk.Frame(tabControl, width=800, height=600)
    tab4 = ttk.Frame(tabControl, width=800, height=600)
    tab5 = ttk.Frame(tabControl, width=800, height=600)

    tab1.grid(row=0, column=0, padx=2, pady=2)
    tab2.grid(row=0, column=0, padx=2, pady=2)
    tab3.grid(row=0, column=0, padx=2, pady=2)
    tab4.grid(row=0, column=0, padx=2, pady=2)
    tab5.grid(row=0, column=0, padx=2, pady=2)

    # ATTENTION - tab3 and tab4 are displayed the other way round in the GUI
    tabControl.add(tab1, text="Folder Analysis")
    tabControl.add(tab2, text="Show Single File")
    tabControl.add(tab4, text="Constant Force Analysis")
    tabControl.add(tab3, text="Advanced Settings")
    tabControl.add(tab5, text="Manual Analysis - TOMATO")

    """ divide the tab1 into frames """
    # output window
    output_frame = tk.Frame(tab1, height=50)
    output_frame.grid(row=0, column=0)
    output_window = tk.Text(output_frame, height=6, width=115)
    output_window.grid(row=0, column=0)
    output_window.insert(
        "end",
        "Welcome to POTATO! \n"
        "Please make sure to select the right datatype -----------------------------------------------------------------> \n"
        "Parameters should be adjusted and validated with <ENTER>.\n"
        "Folders with multiple files can be analysed at once.\n"
    )

    refresh_button = tk.Button(
        output_frame,
        text='Refresh',
        command=refresh,
        bg='#df4c4c',
        activebackground='#eaa90d',
        font='Helvetica 7 bold',
        height=3,
        width=6,
        cursor="exchange"
    )

    refresh_button.grid(row=0, column=1, padx=5)

    # check boxes
    check_box = tk.Frame(tab1)
    check_box.grid(row=0, column=1)

    def select_box(*check_box):
        for i in check_box:
            if i.get() == 1:
                for n in check_box:
                    if not n == i:
                        n.set(value=0)
        boxes = [check_box[x].get() for x in range(len(check_box))]
        if all(boxes) == 0:
            check_box[0].set(value=1)

    check_box_HF = tk.IntVar(value=1)
    check_box_LF = tk.IntVar()
    check_box_CSV = tk.IntVar()
    check_box_Trap1 = tk.IntVar()
    check_box_Trap2 = tk.IntVar(value=1)
    check_box_um = tk.IntVar(value=1)
    check_box_nm = tk.IntVar()
    check_box_multiH5 = tk.IntVar()
    check_box_preprocess = tk.IntVar(value=1)

    check_HF = tk.Checkbutton(
        check_box,
        text="High Frequency (Piezo Distance)",
        variable=check_box_HF,
        command=lambda: [select_box(check_box_HF, check_box_LF, check_box_CSV), parameters(parameter_frame, default_values_HF, default_values_FIT, default_values_constantF)]
    ).grid(row=0, column=0, sticky='W')

    check_LF = tk.Checkbutton(
        check_box,
        text="Low Frequency",
        variable=check_box_LF,
        command=lambda: [select_box(check_box_LF, check_box_HF, check_box_CSV), parameters(parameter_frame, default_values_LF, default_values_FIT, default_values_constantF)]
    ).grid(row=1, column=0, sticky='W')

    check_CSV = tk.Checkbutton(
        check_box,
        text="CSV (F(pN) | d)",
        variable=check_box_CSV,
        command=lambda: [select_box(check_box_CSV, check_box_HF, check_box_LF), parameters(parameter_frame, default_values_CSV, default_values_FIT, default_values_constantF)]
    ).grid(row=2, column=0, sticky='W')

    check_Trap1 = tk.Checkbutton(
        check_box,
        text="Trap 1x",
        variable=check_box_Trap1,
        command=lambda: select_box(check_box_Trap1, check_box_Trap2)
    ).grid(row=0, column=1, padx=8, sticky='W')

    check_Trap2 = tk.Checkbutton(
        check_box,
        text="Trap 2x",
        variable=check_box_Trap2,
        command=lambda: select_box(check_box_Trap2, check_box_Trap1)
    ).grid(row=1, column=1, padx=8, sticky='W')

    check_um = tk.Checkbutton(
        check_box,
        text="µm input",
        variable=check_box_um,
        command=lambda: select_box(check_box_um, check_box_nm)
    ).grid(row=2, column=1, padx=8, sticky='W')

    check_nm = tk.Checkbutton(
        check_box,
        text="nm input",
        variable=check_box_nm,
        command=lambda: select_box(check_box_nm, check_box_um)
    ).grid(row=3, column=1, padx=8, sticky='W')

    check_Multi = tk.Checkbutton(
        check_box,
        text="MultiH5",
        variable=check_box_multiH5
    ).grid(row=4, column=0, sticky='W')

    figure_frame = tk.Canvas(tab1, height=650, width=1000, borderwidth=1, relief='ridge')
    figure_frame.grid(row=1, column=0)

    parameter_frame = tk.Frame(tab1)
    parameter_frame.grid(row=1, column=1, sticky='NE')

    """ parameter frame """
    Cluster_preprocessing = tk.Label(parameter_frame, text='PREPROCESSING', font='Helvetica 9 bold')
    check_preprocess = tk.Checkbutton(
        parameter_frame,
        variable=check_box_preprocess
    ).grid(row=0, column=1, pady=(20, 2), sticky='W')
    Label_downsample = tk.Label(parameter_frame, text='Downsampling rate')
    Label_Filter1 = tk.Label(parameter_frame, text='Butterworth filter degree')
    Label_Filter2 = tk.Label(parameter_frame, text='Cut-off frequency')
    Label_ForceMin = tk.Label(parameter_frame, text='Force threshold, pN')
    Cluster_statistics = tk.Label(parameter_frame, text='STATISTICS', font='Helvetica 9 bold')
    Label_Zscore_F = tk.Label(parameter_frame, text='Z-score force')
    Label_Zscore_D = tk.Label(parameter_frame, text='Z-score distance')

    downsample_value1 = tk.Entry(parameter_frame)
    downsample_value1.bind("<Return>", lambda event: user_input(event, downsample_value1, downsample_value2))

    Filter_degree1 = tk.Entry(parameter_frame)
    Filter_degree1.bind("<Return>", lambda event: user_input(event, Filter_degree1, Filter_degree2))

    Filter_cut_off1 = tk.Entry(parameter_frame)
    Filter_cut_off1.bind("<Return>", lambda event: user_input(event, Filter_cut_off1, Filter_cut_off2))

    Force_Min1 = tk.Entry(parameter_frame)
    Force_Min1.bind("<Return>", lambda event: user_input(event, Force_Min1, Force_Min2))

    Z_score_force1 = tk.Entry(parameter_frame)
    Z_score_force1.bind("<Return>", lambda event: user_input(event, Z_score_force1, Z_score_force2))

    Z_score_distance1 = tk.Entry(parameter_frame)
    Z_score_distance1.bind("<Return>", lambda event: user_input(event, Z_score_distance1, Z_score_distance2))

    Cluster_preprocessing.grid(row=0, column=0, padx=2, pady=(20, 2))
    Label_downsample.grid(row=1, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    downsample_value1.grid(row=1, column=1, padx=2, pady=2)

    Label_Filter1.grid(row=2, column=0, padx=2, pady=2)
    Filter_degree1.grid(row=2, column=1, padx=2, pady=2)

    Label_Filter2.grid(row=3, column=0, padx=2, pady=2)
    Filter_cut_off1.grid(row=3, column=1, padx=2, pady=2)

    Label_ForceMin.grid(row=4, column=0, padx=2, pady=2)
    Force_Min1.grid(row=4, column=1, padx=2, pady=2)

    Cluster_statistics.grid(row=5, column=0, padx=2, pady=(20, 2))
    Label_Zscore_F.grid(row=6, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    Z_score_force1.grid(row=6, column=1, padx=2, pady=2)

    Label_Zscore_D.grid(row=7, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    Z_score_distance1.grid(row=7, column=1, padx=2, pady=2)

    BUTTON1 = tk.Button(
        parameter_frame,
        text='Select Folder to Analyse!',
        command=start_analysis,
        bg='#df4c4c',
        activebackground='#eaa90d',
        font='Helvetica 12 bold',
        height=2,
        width=20
    )

    BUTTON1.grid(row=9, column=0, columnspan=2, pady=125)

    """organize tab2"""
    figure_frame2 = tk.Canvas(tab2, height=650, width=650, borderwidth=1, relief='ridge')
    figure_frame2.grid(row=0, column=0)

    parameter_frame2 = tk.Frame(tab2)
    parameter_frame2.grid(row=0, column=1, sticky='NE')

    BUTTON2 = tk.Button(
        parameter_frame2,
        text='Open h5 file',
        command=getRAW_File_h5,
        bg='#df4c4c',
        activebackground='#eaa90d',
        font='Helvetica 11 bold',
        height=1,
        width=15
    )

    BUTTON2.grid(row=0, column=0, pady=20, sticky='E')

    BUTTON3 = tk.Button(
        parameter_frame2,
        text='Open csv file',
        command=getRAW_File_csv,
        bg='#df4c4c',
        activebackground='#eaa90d',
        font='Helvetica 11 bold',
        height=1,
        width=15
    )

    BUTTON3.grid(row=1, column=0, pady=20, sticky='E')

    """organize tab3 """
    frame1 = tk.Frame(tab3, borderwidth=1, relief='ridge')
    frame1.grid(row=0, column=0)
    frame2 = tk.Frame(tab3, borderwidth=1, relief='ridge')
    frame2.grid(row=0, column=1, sticky='N', padx=(50, 20))
    frame3 = tk.Frame(tab3, borderwidth=1, relief='ridge')
    frame3.grid(row=0, column=2, sticky='N', padx=(50, 20))

    """ parameters in advanced settings """
    Cluster_preprocessing = tk.Label(frame1, text='PREPROCESSING', font='Helvetica 9 bold')
    Label_downsample = tk.Label(frame1, text='Downsampling rate')
    Label_Filter1 = tk.Label(frame1, text='Butterworth filter degree')
    Label_Filter2 = tk.Label(frame1, text='Cut-off frequency')
    Label_ForceMin = tk.Label(frame1, text='Force threshold, pN')
    Cluster_derivative = tk.Label(frame1, text="DERIVATIVE", font='Helvetica 9 bold')
    Label_step_d = tk.Label(frame1, text='Step d')
    Label_Frequency = tk.Label(frame1, text='Data frequency, Hz')
    Cluster_statistics = tk.Label(frame1, text='STATISTICS', font='Helvetica 9 bold')
    Label_Zscore_F = tk.Label(frame1, text='Z-score force')
    Label_Zscore_D = tk.Label(frame1, text='Z-score distance')
    Label_window_size = tk.Label(frame1, text='Moving median window size')
    Label_STD_difference = tk.Label(frame1, text='SD difference threshold')

    downsample_value2 = tk.Entry(frame1)
    downsample_value2.bind("<Return>", lambda event: user_input(event, downsample_value2, downsample_value1))

    Filter_degree2 = tk.Entry(frame1)
    Filter_degree2.bind("<Return>", lambda event: user_input(event, Filter_degree2, Filter_degree1))

    Filter_cut_off2 = tk.Entry(frame1)
    Filter_cut_off2.bind("<Return>", lambda event: user_input(event, Filter_cut_off2, Filter_cut_off1))

    Force_Min2 = tk.Entry(frame1)
    Force_Min2.bind("<Return>", lambda event: user_input(event, Force_Min2, Force_Min1))

    Z_score_force2 = tk.Entry(frame1)
    Z_score_force2.bind("<Return>", lambda event: user_input(event, Z_score_force2, Z_score_force1))

    Z_score_distance2 = tk.Entry(frame1)
    Z_score_distance2.bind("<Return>", lambda event: user_input(event, Z_score_distance2, Z_score_distance1))

    step_d_value = tk.Entry(frame1)
    window_size_value = tk.Entry(frame1)
    STD_difference_value = tk.Entry(frame1)
    Frequency_value = tk.Entry(frame1)

    Cluster_preprocessing.grid(row=0, column=0, padx=2, pady=(20, 2))
    Label_downsample.grid(row=1, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    downsample_value2.grid(row=1, column=1, padx=(0, 20), pady=2)

    Label_Filter1.grid(row=2, column=0, padx=2, pady=2)
    Filter_degree2.grid(row=2, column=1, padx=(0, 20), pady=2)

    Label_Filter2.grid(row=3, column=0, padx=2, pady=2)
    Filter_cut_off2.grid(row=3, column=1, padx=(0, 20), pady=2)

    Label_ForceMin.grid(row=4, column=0, padx=2, pady=2)
    Force_Min2.grid(row=4, column=1, padx=(0, 20), pady=2)

    Cluster_derivative.grid(row=5, column=0, padx=2, pady=(20, 2))
    Label_step_d.grid(row=6, column=0, padx=2, pady=2)
    step_d_value.grid(row=6, column=1, padx=(0, 20), pady=2)

    Label_Frequency.grid(row=7, column=0, padx=2, pady=2)
    Frequency_value.grid(row=7, column=1, padx=(0, 20), pady=2)

    Cluster_statistics.grid(row=8, column=0, padx=2, pady=(20, 2))
    Label_Zscore_F.grid(row=9, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    Z_score_force2.grid(row=9, column=1, padx=(0, 20), pady=2)

    Label_Zscore_D.grid(row=10, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    Z_score_distance2.grid(row=10, column=1, padx=(0, 20), pady=2)

    Label_window_size.grid(row=12, column=0, padx=2, pady=2)
    window_size_value.grid(row=12, column=1, padx=(0, 20), pady=2)

    Label_STD_difference.grid(row=13, column=0, padx=2, pady=2)
    STD_difference_value.grid(row=13, column=1, padx=(0, 20), pady=2)

    """ Output settings """
    check_box_smooth_data = tk.IntVar(value=1)
    check_box_plot = tk.IntVar(value=1)
    check_box_steps = tk.IntVar(value=1)
    check_box_total_results = tk.IntVar(value=1)
    check_box_fitting = tk.IntVar(value=1)

    Label_export = tk.Label(frame2, text="Select exported data", font='Helvetica 9 bold').grid(row=0, column=0, padx=20, pady=20)

    check_1 = tk.Checkbutton(
        frame2,
        text="Processed FD data",
        variable=check_box_smooth_data,
    ).grid(row=1, column=0, sticky='W')

    check_2 = tk.Checkbutton(
        frame2,
        text="Plot",
        variable=check_box_plot,
    ).grid(row=2, column=0, sticky='W')

    check_3 = tk.Checkbutton(
        frame2,
        text="Steps found",
        variable=check_box_steps,
    ).grid(row=3, column=0, sticky='W')

    check_4 = tk.Checkbutton(
        frame2,
        text="Total results (All steps from all files)",
        variable=check_box_total_results,
    ).grid(row=4, column=0, sticky='W')

    check_5 = tk.Checkbutton(
        frame2,
        text="Fitting",
        variable=check_box_fitting,
    ).grid(row=5, column=0, sticky='W')

    """ Fitting parameters """
    Cluster_fitting = tk.Label(frame3, text='FITTING', font='Helvetica 9 bold')
    check_box_WLC = tk.IntVar(value=1)
    check_box_FJC = tk.IntVar(value=0)
    Label_dsLp = tk.Label(frame3, text='dsLp, nm')
    Label_dsLp_up = tk.Label(frame3, text='dsLp, upper bound, nm')
    Label_dsLp_low = tk.Label(frame3, text='dsLp, lower bound, nm')
    Label_dsLc = tk.Label(frame3, text='dsLc, nm')
    Label_ssLp = tk.Label(frame3, text='ssLp, nm')
    Label_ssLc = tk.Label(frame3, text='ssLc, nm')
    Label_stiffness_ds = tk.Label(frame3, text='dsK0, pN')
    Label_stiffness_ds_up = tk.Label(frame3, text='dsK0, upper bound, pN')
    Label_stiffness_ds_low = tk.Label(frame3, text='dsK0, lower bound, pN')
    Label_stiffness_ss = tk.Label(frame3, text='ssK0, pN')
    Label_f_offset = tk.Label(frame3, text='Force offset, pN')
    Label_f_offset_up = tk.Label(frame3, text='Force offset, upper bound, pN')
    Label_f_offset_low = tk.Label(frame3, text='Force offset, lower bound, pN')
    Label_d_offset = tk.Label(frame3, text='Distance offset, nm')
    Label_d_offset_up = tk.Label(frame3, text='Distance offset, upper bound, nm')
    Label_d_offset_low = tk.Label(frame3, text='Distance offset, lower bound, nm')

    dsLp = tk.Entry(frame3)
    dsLp_up = tk.Entry(frame3)
    dsLp_low = tk.Entry(frame3)
    dsLc = tk.Entry(frame3)
    ssLp = tk.Entry(frame3)
    ssLc = tk.Entry(frame3)
    stiff_ds = tk.Entry(frame3)
    stiff_ds_up = tk.Entry(frame3)
    stiff_ds_low = tk.Entry(frame3)
    stiff_ss = tk.Entry(frame3)
    f_off = tk.Entry(frame3)
    f_off_up = tk.Entry(frame3)
    f_off_low = tk.Entry(frame3)
    d_off = tk.Entry(frame3)
    d_off_up = tk.Entry(frame3)
    d_off_low = tk.Entry(frame3)

    Cluster_fitting.grid(row=0, column=0, padx=20, pady=20)

    check_WLC = tk.Checkbutton(
        frame3,
        text="WLC+WLC",
        variable=check_box_WLC,
        command=lambda: [check_box_WLC.set(value=1), check_box_FJC.set(value=0)]
    ).grid(row=1, column=0, sticky='W', pady=20)

    check_FJC = tk.Checkbutton(
        frame3,
        text="WLC+FJC",
        variable=check_box_FJC,
        command=lambda: [check_box_WLC.set(value=0), check_box_FJC.set(value=1)]
    ).grid(row=1, column=1, sticky='W', pady=20)

    Label_dsLp.grid(row=2, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    dsLp.grid(row=2, column=1, padx=(0, 20), pady=2)

    Label_dsLp_up.grid(row=3, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    dsLp_up.grid(row=3, column=1, padx=(0, 20), pady=2)

    Label_dsLp_low.grid(row=4, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    dsLp_low.grid(row=4, column=1, padx=(0, 20), pady=2)

    Label_dsLc.grid(row=5, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    dsLc.grid(row=5, column=1, padx=(0, 20), pady=2)

    Label_ssLp.grid(row=2, column=2, sticky=tk.E + tk.W, padx=2, pady=2)
    ssLp.grid(row=2, column=3, padx=(0, 20), pady=2)

    Label_ssLc.grid(row=3, column=2, sticky=tk.E + tk.W, padx=2, pady=2)
    ssLc.grid(row=3, column=3, padx=(0, 20), pady=2)

    Label_stiffness_ds.grid(row=6, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    stiff_ds.grid(row=6, column=1, padx=(0, 20), pady=2)

    Label_stiffness_ds_up.grid(row=7, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    stiff_ds_up.grid(row=7, column=1, padx=(0, 20), pady=2)

    Label_stiffness_ds_low.grid(row=8, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    stiff_ds_low.grid(row=8, column=1, padx=(0, 20), pady=2)

    Label_stiffness_ss.grid(row=6, column=2, sticky=tk.E + tk.W, padx=2, pady=2)
    stiff_ss.grid(row=6, column=3, padx=(0, 20), pady=2)

    Label_f_offset.grid(row=9, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    f_off.grid(row=9, column=1, padx=(0, 20), pady=2)

    Label_f_offset_up.grid(row=10, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    f_off_up.grid(row=10, column=1, padx=(0, 20), pady=2)

    Label_f_offset_low.grid(row=11, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    f_off_low.grid(row=11, column=1, padx=(0, 20), pady=2)

    Label_d_offset.grid(row=12, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    d_off.grid(row=12, column=1, padx=(0, 20), pady=2)

    Label_d_offset_up.grid(row=13, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    d_off_up.grid(row=13, column=1, padx=(0, 20), pady=2)

    Label_d_offset_low.grid(row=14, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    d_off_low.grid(row=14, column=1, padx=(0, 20), pady=2)

    """organize tab4"""
    # split tab into 2 frames, one for the figure to be displayed and one for the parameters
    figure_frame_tab4 = tk.Canvas(tab4, height=650, width=650, borderwidth=1, relief='ridge')
    figure_frame_tab4.grid(row=0, column=0)

    parameter_frame_tab4 = tk.Frame(tab4)
    parameter_frame_tab4.grid(row=0, column=1, sticky='NE')

    BUTTON1_tab4 = tk.Button(
        parameter_frame_tab4,
        text='Fit Constant Force Data',
        command=start_constantF,
        bg='#df4c4c',
        activebackground='#eaa90d',
        font='Helvetica 11 bold',
        height=1,
        width=25
    )

    BUTTON2_tab4 = tk.Button(
        parameter_frame_tab4,
        text='Display Constant Force Data',
        command=show_constantF,
        bg='#df4c4c',
        activebackground='#eaa90d',
        font='Helvetica 11 bold',
        height=1,
        width=25
    )

    BUTTON1_tab4.grid(row=0, column=0, columnspan=2, padx=20, pady=20, sticky='E')
    BUTTON2_tab4.grid(row=0, column=2, columnspan=2, padx=20, pady=20, sticky='E')

    # organize settings
    Cluster_axes = tk.Label(parameter_frame_tab4, text='SET AXES', font='Helvetica 9 bold')

    Label_x_min = tk.Label(parameter_frame_tab4, text='x min')
    x_min = tk.Entry(parameter_frame_tab4)

    Label_x_max = tk.Label(parameter_frame_tab4, text='x max')
    x_max = tk.Entry(parameter_frame_tab4)

    Label_y_min = tk.Label(parameter_frame_tab4, text='y min')
    y_min = tk.Entry(parameter_frame_tab4)

    Label_y_max = tk.Label(parameter_frame_tab4, text='y max')
    y_max = tk.Entry(parameter_frame_tab4)

    Cluster_expected_fit = tk.Label(parameter_frame_tab4, text='EXPECTED VALUES', font='Helvetica 9 bold')

    Label_number_gauss = tk.Label(parameter_frame_tab4, text='Number of expected gaussians')
    number_gauss = tk.Entry(parameter_frame_tab4)

    Label_mean_gauss = tk.Label(parameter_frame_tab4, text='Expected mean of each gaussian')
    mean_gauss = tk.Entry(parameter_frame_tab4)

    Label_STD_gauss = tk.Label(parameter_frame_tab4, text='Expected SD of each gaussian')
    STD_gauss = tk.Entry(parameter_frame_tab4)

    Label_amplitude_gauss = tk.Label(parameter_frame_tab4, text='Expected amplitude of each gaussian')
    amplitude_gauss = tk.Entry(parameter_frame_tab4)

    Cluster_axes.grid(row=2, column=0, padx=2, pady=(20, 2))
    Label_x_min.grid(row=3, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    x_min.grid(row=3, column=1, padx=(0, 20), pady=2)

    Label_x_max.grid(row=4, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    x_max.grid(row=4, column=1, padx=(0, 20), pady=2)

    Label_y_min.grid(row=3, column=2, sticky=tk.E + tk.W, padx=2, pady=2)
    y_min.grid(row=3, column=3, padx=(0, 20), pady=2)

    Label_y_max.grid(row=4, column=2, sticky=tk.E + tk.W, padx=2, pady=2)
    y_max.grid(row=4, column=3, padx=(0, 20), pady=2)

    Cluster_expected_fit.grid(row=5, column=0, padx=2, pady=(20, 2))
    Label_number_gauss.grid(row=6, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    number_gauss.grid(row=6, column=1, sticky=tk.E + tk.W, padx=2, pady=2)

    Label_mean_gauss.grid(row=7, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    mean_gauss.grid(row=7, column=1, sticky=tk.E + tk.W, padx=2, pady=2)

    Label_STD_gauss.grid(row=8, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    STD_gauss.grid(row=8, column=1, sticky=tk.E + tk.W, padx=2, pady=2)

    Label_amplitude_gauss.grid(row=9, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    amplitude_gauss.grid(row=9, column=1, sticky=tk.E + tk.W, padx=2, pady=2)

    """organize tab5 ---- TOMATO"""
    tab5.columnconfigure([0, 1], weight=1, minsize=75)
    tab5.rowconfigure(0, weight=1, minsize=50)

    canvas1 = tk.Canvas(tab5, width=650, height=700)
    canvas1.grid(row=0, column=1)

    canvas2 = tk.Canvas(tab5, width=650, height=100)
    canvas2.grid(row=1, column=1)

    TOMATO_frame = tk.Frame(tab5, width=500, height=700)
    TOMATO_frame.grid(row=0, column=0)

    canvas_name = tk.Canvas(tab5, width=400, height=50)
    canvas_name.grid(row=1, column=0)

    frame_table = tk.Frame(tab5, width=400, height=100)
    frame_table.grid(row=2, column=0)

    label_shift_d = tk.Label(tab5, text='Shift x [nm]')
    label_shift_d.config(font=('Arial', 10))
    canvas1.create_window(100, 160, window=label_shift_d)

    entryText_shift_d = tk.StringVar()
    entry_shift_d = tk.Entry(tab5, textvariable=entryText_shift_d)
    canvas1.create_window(100, 180, window=entry_shift_d)
    entryText_shift_d.set("0")

    # K0 for both
    # ds
    label_ds_St = tk.Label(tab5, text='K0 ds (St)')
    label_ds_St.config(font=('Arial', 10))
    canvas1.create_window(100, 200, window=label_ds_St)

    entryText_ds_St = tk.StringVar()
    entry_ds_St = tk.Entry(tab5, textvariable=entryText_ds_St)
    canvas1.create_window(100, 220, window=entry_ds_St)
    entryText_ds_St.set("450")

    # ss
    label_ss_St = tk.Label(tab5, text='K0 ss (St)')
    label_ss_St.config(font=('Arial', 10))
    canvas1.create_window(300, 200, window=label_ss_St)

    entryText_ss_St = tk.StringVar()
    entry_ss_St = tk.Entry(tab5, textvariable=entryText_ss_St)
    canvas1.create_window(300, 220, window=entry_ss_St)
    entryText_ss_St.set("800")

    # shift in F
    label_shift_F = tk.Label(tab5, text='shift F [pN]')
    label_shift_F.config(font=('Arial', 10))
    canvas1.create_window(300, 160, window=label_shift_F)

    entryText_shift_F = tk.StringVar()
    entry_shift_F = tk.Entry(tab5, textvariable=entryText_shift_F)
    canvas1.create_window(300, 180, window=entry_shift_F)
    entryText_shift_F.set("0")

    ## ds handle part
    # ds handle persistance length
    label_ds_Lp = tk.Label(tab5, text='dsHandles Lp [nm]')
    label_ds_Lp.config(font=('Arial', 10))
    canvas1.create_window(100, 240, window=label_ds_Lp)

    entryText_ds_Lp = tk.StringVar()
    entry_ds_Lp = tk.Entry(tab5, textvariable=entryText_ds_Lp)
    canvas1.create_window(100, 260, window=entry_ds_Lp)
    entryText_ds_Lp.set("40")

    # ds handle  contour length

    label_ds_Lc = tk.Label(tab5, text='dsHandles Lc [nm]')
    label_ds_Lc.config(font=('Arial', 10))
    canvas1.create_window(100, 290, window=label_ds_Lc)

    entryText_ds_Lc = tk.StringVar()
    entry_ds_Lc = tk.Entry(tab5, textvariable=entryText_ds_Lc)
    canvas1.create_window(100, 310, window=entry_ds_Lc)
    entryText_ds_Lc.set("1256")

    ## ss RNA part
    # ss RNA persistance length
    label_ss_Lp = tk.Label(tab5, text=' ssRNA Lp [nm]')
    label_ss_Lp.config(font=('Arial', 10))
    canvas1.create_window(300, 240, window=label_ss_Lp)

    entryText_ss_Lp = tk.StringVar()
    entry_ss_Lp = tk.Entry(tab5, textvariable=entryText_ss_Lp)
    canvas1.create_window(300, 260, window=entry_ss_Lp)
    entryText_ss_Lp.set("1")

    # ss RNA contour length
    label_ss_Lc = tk.Label(tab5, text=' ssRNA Lc [nm]')
    label_ss_Lc.config(font=('Arial', 10))
    canvas1.create_window(300, 290, window=label_ss_Lc)

    entryText_ss_Lc = tk.StringVar()
    entry_ss_Lc = tk.Entry(tab5, textvariable=entryText_ss_Lc)
    canvas1.create_window(300, 310, window=entry_ss_Lc)
    entryText_ss_Lc.set("0")

    # work done
    # ds work done
    label_dsWork = tk.Label(tab5, text='ds Work done [KbT]')
    label_dsWork.config(font=('Arial', 10))
    canvas1.create_window(520, 160, window=label_dsWork)

    entryText_dsWork = tk.StringVar()
    entry_dsWork = tk.Entry(tab5, textvariable=entryText_dsWork)
    canvas1.create_window(520, 180, window=entry_dsWork)
    entryText_dsWork.set("0")

    # ss work done
    label_ssWork = tk.Label(tab5, text='ss Work done [KbT]')
    label_ssWork.config(font=('Arial', 10))
    canvas1.create_window(520, 200, window=label_ssWork)

    entryText_ssWork = tk.StringVar()
    entry_ssWork = tk.Entry(tab5, textvariable=entryText_ssWork)
    canvas1.create_window(520, 220, window=entry_ssWork)
    entryText_ssWork.set("0")

    # rectangle work done
    label_rWork = tk.Label(tab5, text='rect. Work done [KbT]')
    label_rWork.config(font=('Arial', 10))
    canvas1.create_window(520, 240, window=label_rWork)

    entryText_rWork = tk.StringVar()
    entry_rWork = tk.Entry(tab5, textvariable=entryText_rWork)
    canvas1.create_window(520, 260, window=entry_rWork)
    entryText_rWork.set("0")

    #  work done by structure

    label_strWork = tk.Label(tab5, text='str. Work done [KbT]')
    label_strWork.config(font=('Arial', 10))
    canvas1.create_window(520, 290, window=label_strWork)

    entryText_strWork = tk.StringVar()
    entry_strWork = tk.Entry(tab5, textvariable=entryText_strWork)
    canvas1.create_window(520, 310, window=entry_strWork)
    entryText_strWork.set("0")

    # start position
    entryText_start = tk.StringVar()
    entry_start = tk.Entry(tab5, textvariable=entryText_start)
    canvas1.create_window(250, 350, window=entry_start, width=50, height=20)

    # end position
    entryText_end = tk.StringVar()
    entry_end = tk.Entry(tab5, textvariable=entryText_end)
    canvas1.create_window(250, 390, window=entry_end, width=50, height=20)

    ## Work start
    # start position
    label_start_work_D = tk.Label(tab5, text='D')
    label_start_work_D.config(font=('Arial', 10))
    canvas1.create_window(500, 330, window=label_start_work_D)

    entryText_start_work_D = tk.StringVar()
    entry_start_work_D = tk.Entry(tab5, textvariable=entryText_start_work_D)
    canvas1.create_window(500, 350, window=entry_start_work_D, width=50, height=20)
    entryText_start_work_D.set('0')

    label_start_work_F = tk.Label(tab5, text='F')
    label_start_work_F.config(font=('Arial', 10))
    canvas1.create_window(550, 330, window=label_start_work_F)

    entryText_start_work_F = tk.StringVar()
    entry_start_work_F = tk.Entry(tab5, textvariable=entryText_start_work_F)
    canvas1.create_window(550, 350, window=entry_start_work_F, width=50, height=20)
    entryText_start_work_F.set('0')

    ## work end
    # end positio
    label_start_work_D = tk.Label(tab5, text='D')
    label_start_work_D.config(font=('Arial', 10))
    canvas1.create_window(500, 390, window=label_start_work_D)

    entryText_end_work_D = tk.StringVar()
    entry_end_work_D = tk.Entry(tab5, textvariable=entryText_end_work_D)
    canvas1.create_window(500, 390, window=entry_end_work_D, width=50, height=20)
    entryText_end_work_D.set('0')

    label_start_work_F = tk.Label(tab5, text='F')
    label_start_work_F.config(font=('Arial', 10))
    canvas1.create_window(550, 390, window=label_start_work_F)

    entryText_end_work_F = tk.StringVar()
    entry_end_work_F = tk.Entry(tab5, textvariable=entryText_end_work_F)
    canvas1.create_window(550, 390, window=entry_end_work_F, width=50, height=20)
    entryText_end_work_F.set('0')

    entryText_filename = tk.StringVar()
    entry_filename = tk.Entry(tab5, textvariable=entryText_filename)
    canvas_name.create_window(200, 30, window=entry_filename, width=800, height=30)

    ## shortcut descriptions
    label_shortcut = tk.Label(tab5, text='Keyboard shortcuts:')
    label_shortcut.config(font=('Arial', 10))
    canvas1.create_window(100, 650, window=label_shortcut)

    label_shortcut_arrows = tk.Label(tab5, text='Left/right/a/d arrow - switching between curves')
    label_shortcut_arrows.config(font=('Arial', 10))
    canvas1.create_window(180, 680, window=label_shortcut_arrows)

    label_shortcut_save = tk.Label(tab5, text='Enter or Shift_L - Save')
    label_shortcut_save.config(font=('Arial', 10))
    canvas1.create_window(105, 710, window=label_shortcut_save)

    label_shortcut_start = tk.Label(tab5, text='s - fit start')
    label_shortcut_start.config(font=('Arial', 10))
    canvas2.create_window(60, 30, window=label_shortcut_start)

    label_shortcut_end = tk.Label(tab5, text='e - fit end')
    label_shortcut_end.config(font=('Arial', 10))
    canvas2.create_window(60, 55, window=label_shortcut_end)

    label_shortcut_work_start = tk.Label(tab5, text='space+s - work start')
    label_shortcut_work_start.config(font=('Arial', 10))
    canvas2.create_window(290, 30, window=label_shortcut_work_start)

    label_shortcut_work_end = tk.Label(tab5, text='space+e - work end')
    label_shortcut_work_end.config(font=('Arial', 10))
    canvas2.create_window(290, 55, window=label_shortcut_work_end)

    label_shortcut_zero = tk.Label(tab5, text='o or ;/° - zero str work')
    label_shortcut_zero.config(font=('Arial', 10))
    canvas2.create_window(105, 80, window=label_shortcut_zero)

    label_shortcut_fit_ds = tk.Label(tab5, text='space+f - fit ds')
    label_shortcut_fit_ds.config(font=('Arial', 10))
    canvas2.create_window(80, 105, window=label_shortcut_fit_ds)

    label_shortcut_fit_ss = tk.Label(tab5, text='space+g - fit ss')
    label_shortcut_fit_ss.config(font=('Arial', 10))
    canvas2.create_window(270, 105, window=label_shortcut_fit_ss)

    label_shortcut_rWork = tk.Label(tab5, text='r - rWork')
    label_shortcut_rWork.config(font=('Arial', 10))
    canvas2.create_window(450, 30, window=label_shortcut_rWork)

    label_shortcut_start = tk.Label(tab5, text='t - strWork')
    label_shortcut_start.config(font=('Arial', 10))
    canvas2.create_window(455, 55, window=label_shortcut_start)

    ## create button widgets that use the defined functions
    browseButton_CSV = tk.Button(tab5, text="      Choose folder     ", command=open_folder, bg='green', fg='white', font=('Arial', 11, 'bold'))
    canvas1.create_window(200, 50, window=browseButton_CSV)

    button_create = tk.Button(tab5, text=' Create Charts ', command=create_chart, bg='palegreen2', font=('Arial', 11, 'bold'))
    canvas1.create_window(200, 90, window=button_create)

    button_clear = tk.Button(tab5, text='  Clear Charts  ', command=clear_charts, bg='lightskyblue2', font=('Arial', 11, 'bold'))
    canvas1.create_window(200, 130, window=button_clear)

    button_start = tk.Button(tab5, text='Set start', command=start_click, bg='lightsteelblue2', font=('Arial', 11, 'bold'))
    canvas1.create_window(150, 350, window=button_start)

    button_end = tk.Button(tab5, text='Set end', command=end_click, bg='lightsteelblue2', font=('Arial', 11, 'bold'))
    canvas1.create_window(150, 390, window=button_end)

    button_fit_Lp_shift_x = tk.Button(tab5, text='Fit ds Lp & shift_x', command=Fitting_WLC_ds_handles, bg='PeachPuff', font=('Arial', 10, 'bold'))
    canvas1.create_window(100, 450, window=button_fit_Lp_shift_x, width=150)

    button_fit_Lc = tk.Button(tab5, text='Fit ss Lc', command=Fitting_WLC_ss_handles, bg='PeachPuff', font=('Arial', 10, 'bold'))
    canvas1.create_window(300, 450, window=button_fit_Lc, width=90)

    button_export = tk.Button(tab5, text='Export', command=export_table, bg='palegreen2', font=('Arial', 14, 'bold'))
    canvas1.create_window(100, 550, window=button_export)

    button_clear_last = tk.Button(tab5, text='Delete last', command=clear_table_last, bg='red', font=('Arial', 10, 'bold'))
    canvas1.create_window(250, 520, window=button_clear_last)

    button_clear_table = tk.Button(tab5, text='Delete all', command=clear_table, bg='red', font=('Arial', 10, 'bold'))
    canvas1.create_window(250, 600, window=button_clear_table)

    button_reset_parameters = tk.Button(tab5, text='Reset prmtrs', command=reset_parameters, bg='PeachPuff', font=('Arial', 10, 'bold'))
    canvas1.create_window(400, 100, window=button_reset_parameters)

    button_start_work = tk.Button(tab5, text='Set W start', command=start_work_click, bg='lightsteelblue2', font=('Arial', 11, 'bold'))
    canvas1.create_window(400, 350, window=button_start_work)

    button_end_work = tk.Button(tab5, text='Set W end', command=end_work_click, bg='lightsteelblue2', font=('Arial', 11, 'bold'))
    canvas1.create_window(400, 390, window=button_end_work)

    button_rWork = tk.Button(tab5, text='rWork', command=calc_rWork, bg='PeachPuff', font=('Arial', 10, 'bold'))
    canvas1.create_window(500, 450, window=button_rWork, width=90)

    button_rWork = tk.Button(tab5, text='rWork', command=calc_rWork, bg='PeachPuff', font=('Arial', 10, 'bold'))
    canvas1.create_window(500, 450, window=button_rWork, width=90)

    button_strWork = tk.Button(tab5, text='strWork', command=calc_strWork, bg='PeachPuff', font=('Arial', 10, 'bold'))
    canvas1.create_window(500, 500, window=button_strWork, width=90)

    button_next_FD = tk.Button(tab5, text='Next FD >', command=lambda: change_FD(1), bg='lightsteelblue2', font=('Arial', 10, 'bold'))
    canvas1.create_window(500, 50, window=button_next_FD, width=90)

    root.bind("<Right>", next_FD_key)
    root.bind("<d>", next_FD_key)
    root.bind("<Left>", previous_FD_key)
    root.bind("<a>", previous_FD_key)
    root.bind("<Return>", save_key)
    root.bind("<Shift_L>", save_key)
    root.bind("<s>", start_click_key)
    root.bind("<e>", end_click_key)
    root.bind("<space><s>", start_work_click_key)
    root.bind("<space><e>", end_work_click_key)
    root.bind("<o>", zero_str_work_key)
    root.bind("<;>", zero_str_work_key)
    root.bind("<space><f>", fit_ds_key)
    root.bind("<space><g>", fit_ss_key)
    root.bind("<r>", calc_rWork_key)
    root.bind("<t>", calc_strWork_key)
    root.bind("<Control-z>", load_previous_data_key)

    button_previous_FD = tk.Button(tab5, text='< Prev. FD', command=lambda: change_FD(-1), bg='lightsteelblue2', font=('Arial', 10, 'bold'))
    canvas1.create_window(400, 50, window=button_previous_FD, width=90)

    button_save = tk.Button(tab5, text='Save', command=write_to_table, bg='palegreen2', font=('Arial', 14, 'bold'))
    canvas1.create_window(400, 500, window=button_save, width=90)

    button_export = tk.Button(tab5, text='Export model', command=export_model, bg='palegreen2', font=('Arial', 10, 'bold'))
    canvas1.create_window(100, 600, window=button_export)

    button_zero = tk.Button(tab5, text='0', command=zero_str_work, bg='palegreen2', font=('Arial', 14, 'bold'))
    canvas1.create_window(620, 300, window=button_zero, width=30)

    ## show the fitting parameters in a table
    # create Treeview with 3 columns
    cols = ('Filename', 'F1', 'F2', 'F1/2', 'Step start', 'Step end', 'Step length', 'ds Lc', 'ds Lp', 'ds St', 'ss Lc', 'ss Lp', 'ss St', 'Shift F', 'Shift x', 'Work')
    listBox = ttk.Treeview(frame_table, columns=cols, show='headings', height=5)
    # set column headings
    for col in cols:
        listBox.heading(col, text=col)
        listBox.column(col, minwidth=0, width=65)
    listBox.grid(row=1, column=0, columnspan=1, padx=5, pady=5)
    ######### TOMATO end ############

    ############ POTATO last part ###############
    # put default values into the widgets
    parameters(parameter_frame, default_values_HF, default_values_FIT, default_values_constantF)

    # loop ensuring the GUI is running until closed
    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()
