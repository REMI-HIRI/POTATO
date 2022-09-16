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
from POTATO_preprocessing import preprocess_RAW, create_derivative
from POTATO_config import default_values_HF, default_values_LF, default_values_CSV, default_values_FIT, default_values_constantF
from POTATO_constantF import get_constantF, display_constantF, fit_constantF
from POTATO_fitting import fitting_ds, fitting_ss
from POTATO_find_steps import calc_integral

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
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()

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
def parameters(default_values, default_fit, default_constantF):
    if not default_values == 0:
        downsample_value.set(default_values['Downsampling rate'])
        Filter_degree.set(default_values['Butterworth filter degree'])
        Filter_cut_off.set(default_values['Cut-off frequency'])
        Force_Min.set(default_values['Force threshold, pN'])
        Z_score_force.set(default_values['Z-score force'])
        Z_score_distance.set(default_values['Z-score distance'])

        step_d_variable.set(str(default_values['Step d']))
        window_size_variable.set(str(default_values['Moving median window size']))
        STD_difference_variable.set(str(default_values['STD difference threshold']))
        Frequency_variable.set(str(default_values['Data frequency, Hz']))

    dsLp_variable.set(str(default_fit['Persistance-Length ds, nm']))
    dsLp_up_variable.set(str(default_fit['Persistance-Length ds, upper bound, nm']))
    dsLp_low_variable.set(str(default_fit['Persistance-Length ds, lower bound, nm']))
    ssLp_variable.set(str(default_fit['Persistance-Length ss, nm']))
    dsLc_variable.set(str(default_fit['Contour-Length ds, nm']))
    ssLc_variable.set(str(default_fit['Persistance-Length ss, nm']))
    stiff_ds_variable.set(str(default_fit['Stiffness ds, pN']))
    stiff_ds_up_variable.set(str(default_fit['Stiffness ds, upper bound, pN']))
    stiff_ds_low_variable.set(str(default_fit['Stiffness ds, lower bound, pN']))
    stiff_ss_variable.set(str(default_fit['Stiffness ss, pN']))
    f_off_variable.set(str(default_fit['Force offset, pN']))
    f_off_up_variable.set(str(default_fit['Force offset, upper bound, pN']))
    f_off_low_variable.set(str(default_fit['Force offset, lower bound, pN']))
    d_off_variable.set(str(default_fit['Distance offset, nm']))
    d_off_up_variable.set(str(default_fit['Distance offset, upper bound, nm']))
    d_off_low_variable.set(str(default_fit['Distance offset, lower bound, nm']))

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

    return input_settings, input_format, export_data, input_fitting, input_constantF


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
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
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
        parameters(default_values_CSV, default_values_FIT, default_values_constantF)
    else:
        pass

    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
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
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
    Force_Distance, Force_Distance_um, frequency, filename, analysis_path, timestamp = get_constantF(input_settings, input_format, input_constantF)
    fig_constantF, hist_D, filteredDistance_ready = display_constantF(Force_Distance, Force_Distance_um, frequency, input_settings, input_constantF)
    os.mkdir(analysis_path)
    export_settings(analysis_path, timestamp, input_settings, input_constantF)
    fig_constantF = fit_constantF(hist_D, Force_Distance, filteredDistance_ready, frequency, input_settings, input_constantF, filename, timestamp)
    fig_constantF_tk = FigureCanvasTkAgg(fig_constantF, figure_frame_tab4)
    fig_constantF_tk.get_tk_widget().grid(row=0, column=0, sticky='wens')

    tabControl.select(tab4)


def show_constantF():
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
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
from POTATO_TOMATO import plot_TOMATO


############# define the functions for TOMATO ##################
def open_folder():
    global filename_TOMATO
    global Force_Distance_TOMATO
    global der_arr_TOMATO
    global TOMATO_fig1
    global Files
    global FD_number
    # check user input
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()

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
    der_arr_TOMATO = create_derivative(input_settings, Frequency_value, Force_Distance_TOMATO[:, 0], Force_Distance_TOMATO[:, 1], 0)

    entryText_filename.set(filename_TOMATO)

    fig = plot_TOMATO(Force_Distance_TOMATO)
    TOMATO_fig1 = FigureCanvasTkAgg(fig, TOMATO_figure_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0, sticky='wens')
    toolbarFrame = tk.Frame(master=TOMATO_figure_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)


def change_FD(direction):
    global TOMATO_fig1
    global filename_TOMATO
    global FD_number
    global Force_Distance_TOMATO
    global orientation

    FD_number = FD_number + direction
    if FD_number == len(Files):
        FD_number = FD_number - len(Files)
    if FD_number < 0:
        FD_number = FD_number + len(Files)

    delete_all_steps()
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
    Force_Distance_TOMATO, Force_Distance_um_TOMATO, Frequency_value, filename_TOMATO = read_in_data(FD_number, Files, input_settings, input_format)

    orientation = 'forward'
    if Force_Distance_TOMATO[0, 1] > Force_Distance_TOMATO[-1, 1]:  # reverse
        Force_Distance_TOMATO = np.flipud(Force_Distance_TOMATO)
        Force_Distance_um_TOMATO = np.flipud(Force_Distance_um_TOMATO)
        orientation = 'reverse'

    entryText_filename.set(filename_TOMATO)

    parameters(0, default_values_FIT, default_values_constantF)

    TOMATO_fig1.get_tk_widget().destroy()

    fig = plot_TOMATO(Force_Distance_TOMATO)
    TOMATO_fig1 = FigureCanvasTkAgg(fig, TOMATO_figure_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0, sticky='wens')
    toolbarFrame = tk.Frame(master=TOMATO_figure_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)


# key binding wrapper functions
def previous_FD_key(event):
    change_FD(-1)


def next_FD_key(event):
    change_FD(+1)


def save_step_key(event):
    save_step()


def start_analysis_key(event):
    analyze_steps()


def start_click_key(event):
    start_click()


def end_click_key(event):
    end_click()


def start_click():
    global cid
    cid = TOMATO_fig1.mpl_connect('button_press_event', lambda event, arg=1: onclick_start_end(event, arg))


def end_click():
    global cid
    cid = TOMATO_fig1.mpl_connect('button_press_event', lambda event, arg=0: onclick_start_end(event, arg))


def onclick_start_end(event, pos):
    global cid

    PD_position, F_position = float(event.xdata), float(event.ydata)

    if pos == 1:
        entryText_startF.set(round(F_position, 1))
        entryText_startD.set(round(PD_position, 1))
    elif pos == 0:
        entryText_endF.set(round(F_position, 1))
        entryText_endD.set(round(PD_position, 1))
    TOMATO_fig1.mpl_disconnect(cid)


def save_step():
    global step_number
    try:
        tree_steps.insert('', 'end', values=(step_number, entryText_startF.get(), entryText_endF.get(), entryText_startD.get(), entryText_endD.get()))
        step_number += 1
    except:
        print('Please make sure step start and step end are selected!')


def analyze_steps():
    global TOMATO_fig1
    global subplot1
    global tree_results

    # get input settings
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
    timestamp = time.strftime("%Y%m%d-%H%M%S")

    # write step list into pandas dataframe
    row_list = []
    columns = ('Step number', 'F start', 'F end', 'Step start', 'Step end')
    for row in tree_steps.get_children():
        row_list.append(tree_steps.item(row)["values"])
    treeview_df = pd.DataFrame(row_list, columns=columns)

    # iterate through dataframe and fit each part of the curve
    TOMATO_fig1.get_tk_widget().destroy()

    figure1 = plot_TOMATO(Force_Distance_TOMATO)
    diff_colors = ['b', 'r', 'c', 'g', 'y', 'm', 'b', 'r', 'c', 'g', 'y', 'm']
    subplot1 = figure1.add_subplot(111)
    subplot1.plot(Force_Distance_TOMATO[:, 1], Force_Distance_TOMATO[:, 0], color='gray')
    distance = np.arange(min(Force_Distance_TOMATO[:, 1]), max(Force_Distance_TOMATO[:, 1]) + 50, 2)

    export_fit = []
    fit = []
    start_force_ss = []
    start_distance_ss = []
    integral_ss_fit_start = []
    integral_ss_fit_end = []

    for i in treeview_df.index:
        # part before first step is fitted with a single WLC model (ds part)
        if treeview_df['Step number'][i] == 1:
            j = i
            ds_fit_dict_TOMATO, TOMATO_area_ds, real_step_start = fitting_ds(filename_TOMATO, input_settings, export_data, input_fitting, float(treeview_df['Step start'][i]), Force_Distance_TOMATO, der_arr_TOMATO, [], 1)
            ds_fit_region_end = real_step_start

            dsLp_variable.set(ds_fit_dict_TOMATO['Lp_ds'])
            f_off_variable.set(ds_fit_dict_TOMATO['f_offset'])
            d_off_variable.set(ds_fit_dict_TOMATO["d_offset"])
            dsLc_variable.set(ds_fit_dict_TOMATO['Lc_ds'])
            stiff_ds_variable.set(ds_fit_dict_TOMATO['St_ds'])

            tree_results.insert("", "end", iid='{}no step'.format(timestamp), values=(
                filename_TOMATO,
                i,
                '',
                '',
                '',
                '',
                '',
                '',
                dsLc_variable.get(),
                dsLp_variable.get(),
                stiff_ds_variable.get(),
                '',
                '',
                '',
                f_off_variable.get(),
                d_off_variable.get(),
                '',
                ''
            )
            )

            export_fit.append(ds_fit_dict_TOMATO)

            F_ds_model = ds_fit_dict_TOMATO['model_ds'](distance, ds_fit_dict_TOMATO['fit_model'].params)
            # plot the marked ds region and fits
            subplot1.plot(Force_Distance_TOMATO[:, 1][:real_step_start], Force_Distance_TOMATO[:, 0][:real_step_start], color=diff_colors[i])
            subplot1.plot(distance, F_ds_model, marker=None, linestyle='dashed', linewidth=1, color="black")

        # fit the other ss parts
        elif treeview_df['Step number'][i] > 1:
            j = i
            fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict_TOMATO, area_ss_fit_start, area_ss_fit_end = fitting_ss(filename_TOMATO, input_settings, export_data, input_fitting, float(treeview_df['Step end'][i - 1]), float(treeview_df['Step start'][i]), Force_Distance_TOMATO, 1, 1, der_arr_TOMATO, [], 1)

            real_step_start = np.where(Force_Distance_TOMATO[:, 0] == f_fitting_region_ss[-1])
            real_step_start = real_step_start[0][0]

            fit.append(fit_ss)
            start_force_ss.append(f_fitting_region_ss)
            start_distance_ss.append(d_fitting_region_ss)
            export_fit.append(ss_fit_dict_TOMATO)
            integral_ss_fit_start.append(area_ss_fit_start)
            integral_ss_fit_end.append(area_ss_fit_end)

            ssLp_variable.set(ss_fit_dict_TOMATO['Lp_ss'])
            f_off_variable.set(ss_fit_dict_TOMATO['f_offset'])
            d_off_variable.set(ss_fit_dict_TOMATO["d_offset"])
            ssLc_variable.set(ss_fit_dict_TOMATO['Lc_ss'])
            stiff_ss_variable.set(ss_fit_dict_TOMATO['St_ss'])

            tree_results.insert("", "end", iid='{}step{}'.format(timestamp, j), values=(
                filename_TOMATO,
                i,
                Force_Distance_TOMATO[real_step_start, 0],
                f_fitting_region_ss[0],
                f_fitting_region_ss[0] - Force_Distance_TOMATO[real_step_start, 0],
                Force_Distance_TOMATO[real_step_start, 1],
                d_fitting_region_ss[0],
                d_fitting_region_ss[0] - Force_Distance_TOMATO[real_step_start, 1],
                '',
                '',
                '',
                ssLc_variable.get(),
                ssLp_variable.get(),
                stiff_ss_variable.get(),
                f_off_variable.get(),
                d_off_variable.get(),
                '',
                ''
            )
            )

            # plot the marked regions and fits
            # model data
            F_ss_model = ss_fit_dict_TOMATO['model_ss'](distance, fit_ss.params)

            # plot the marked ss region and fits
            subplot1.plot(d_fitting_region_ss[:], f_fitting_region_ss, color=diff_colors[i])
            subplot1.plot(distance, F_ss_model, marker=None, linewidth=1, linestyle='dashed', color="black")

    # fit the last part of the curve
    fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict_TOMATO, area_ss_fit_start, area_ss_fit_end = fitting_ss(
        filename_TOMATO,
        input_settings,
        export_data,
        input_fitting,
        float(treeview_df['Step end'][len(treeview_df) - 1]),
        max(Force_Distance_TOMATO[:, 1]),
        Force_Distance_TOMATO,
        1,
        1,
        der_arr_TOMATO,
        [],
        1
    )

    fit.append(fit_ss)
    start_force_ss.append(f_fitting_region_ss)
    start_distance_ss.append(d_fitting_region_ss)
    export_fit.append(ss_fit_dict_TOMATO)
    integral_ss_fit_start.append(area_ss_fit_start)
    integral_ss_fit_end.append(area_ss_fit_end)

    ssLp_variable.set(ss_fit_dict_TOMATO['Lp_ss'])
    f_off_variable.set(ss_fit_dict_TOMATO['f_offset'])
    d_off_variable.set(ss_fit_dict_TOMATO["d_offset"])
    ssLc_variable.set(ss_fit_dict_TOMATO['Lc_ss'])
    stiff_ss_variable.set(ss_fit_dict_TOMATO['St_ss'])

    tree_results.insert("", "end", iid='{}step{}'.format(timestamp, j + 1), values=(
        filename_TOMATO,
        j + 1,
        Force_Distance_TOMATO[:, 0][real_step_start],
        f_fitting_region_ss[0],
        f_fitting_region_ss[0] - Force_Distance_TOMATO[:, 0][real_step_start],
        Force_Distance_TOMATO[:, 1][real_step_start],
        d_fitting_region_ss[0],
        d_fitting_region_ss[0] - Force_Distance_TOMATO[:, 1][real_step_start],
        '',
        '',
        '',
        ssLc_variable.get(),
        ssLp_variable.get(),
        stiff_ss_variable.get(),
        f_off_variable.get(),
        d_off_variable.get(),
        '',
        ''
    )
    )

    work_first_step, kT_1 = calc_integral(
        TOMATO_area_ds,
        integral_ss_fit_start[0],
        Force_Distance_TOMATO[ds_fit_region_end, 1],
        start_distance_ss[0][0],
        Force_Distance_TOMATO[ds_fit_region_end, 0],
        start_force_ss[0][0]
    )

    tree_results.set('{}step1'.format(timestamp), column='Work [pN*nm]', value=work_first_step)
    tree_results.set('{}step1'.format(timestamp), column='Work [kT]', value=kT_1)

    if j > 1:
        for n in range(1, j):
            work_step_n, kT_n = calc_integral(
                integral_ss_fit_end[n],
                integral_ss_fit_start[n - 1],
                start_distance_ss[n - 1][-1],
                start_distance_ss[n][0],
                start_force_ss[n - 1][-1],
                start_force_ss[n][0]
            )

            tree_results.set('{}step{}'.format(timestamp, n), column='Work [pN*nm]', value=work_step_n)
            tree_results.set('{}step{}'.format(timestamp, n), column='Work [kT]', value=kT_n)

    # plot the marked regions and fits
    # model data
    F_ss_model = ss_fit_dict_TOMATO['model_ss'](distance, fit_ss.params)

    # plot the marked ss region and fits
    subplot1.plot(d_fitting_region_ss[:], f_fitting_region_ss, color=diff_colors[j + 1])
    subplot1.plot(distance, F_ss_model, marker=None, linewidth=1, linestyle='dashed', color="black")

    subplot1.set_ylim([min(Force_Distance_TOMATO[:, 0]), max(Force_Distance_TOMATO[:, 0])])
    subplot1.set_xlim([min(Force_Distance_TOMATO[:, 1]) - 10, max(Force_Distance_TOMATO[:, 1]) + 10])
    subplot1.tick_params('both', direction='in')

    TOMATO_fig1 = FigureCanvasTkAgg(figure1, TOMATO_figure_frame)
    TOMATO_fig1.get_tk_widget().grid(row=0, column=0)

    toolbarFrame = tk.Frame(master=TOMATO_figure_frame)
    toolbarFrame.grid(row=2, column=0)
    toolbar = NavigationToolbar2Tk(TOMATO_fig1, toolbarFrame)


def delete_step():
    global step_number
    list_items = tree_steps.get_children("")
    tree_steps.delete(list_items[-1])
    step_number -= 1


def delete_all_steps():
    global step_number

    tree_steps.delete(*tree_steps.get_children())
    step_number = 1


def delete_result(event):
    selected_items = tree_results.selection()
    for selected_item in selected_items:
        tree_results.delete(selected_item)


def clear_table():
    global tree_results
    list_items = tree_results.get_children("")

    for item in list_items:
        tree_results.delete(item)


def clear_table_last():
    global tree_results
    list_items = tree_results.get_children("")

    tree_results.delete(list_items[-1])


def export_table():
    global tree_results
    global name
    global Fit_results
    ''' exporting the table results '''
    results = []
    for child in tree_results.get_children():
        results.append(tree_results.item(child)['values'])

    Fit_results = pd.DataFrame(results,
                            columns=[
                                'Filename',
                                'step number',
                                'Force step start [pN]',
                                'Force step end [pN]',
                                'delta force [pN]',
                                'extension step start [nm]',
                                'extension step end [nm]',
                                'Step length [nm]',
                                'ds contour length',
                                'ds persistance Length',
                                'ds stiffness (K0) [pN]',
                                'ss contour Length',
                                'ss persistance Length',
                                'ss stiffness (K0) [pN]',
                                'Force offset',
                                'Distance offset',
                                'Work [pN*nm]',
                                'Work [kT]'
                            ]
    )

    name = filedialog.asksaveasfile(mode='w', defaultextension=".csv")
    Fit_results.to_csv(name.name, index=False, header=True)


def tab_bind(event=None):
    if tabControl.index(tabControl.select()) == 4:
        root.bind("<Right>", next_FD_key)
        root.bind("<d>", next_FD_key)
        root.bind("<Left>", previous_FD_key)
        root.bind("<a>", previous_FD_key)
        root.bind("<s>", start_click_key)
        root.bind("<e>", end_click_key)
        root.bind("<Control-s>", save_step_key)
        root.bind("<f>", start_analysis_key)
        root.bind("<Delete>", delete_result)
    else:
        root.unbind("<Right>")
        root.unbind("<d>")
        root.unbind("<Left>")
        root.unbind("<a>")
        root.unbind("<s>")
        root.unbind("<e>")
        root.unbind("<Control-s>")
        root.unbind("<f>")
        root.unbind("<Delete>")

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
    tabControl.pack(expand=4, fill='both')
    root.bind('<<NotebookTabChanged>>', tab_bind)

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
        "Parameters should be adjusted prior to analysis.\n"
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
        command=lambda: [select_box(check_box_HF, check_box_LF, check_box_CSV), parameters(default_values_HF, default_values_FIT, default_values_constantF)]
    ).grid(row=0, column=0, sticky='W')

    check_LF = tk.Checkbutton(
        check_box,
        text="Low Frequency",
        variable=check_box_LF,
        command=lambda: [select_box(check_box_LF, check_box_HF, check_box_CSV), parameters(default_values_LF, default_values_FIT, default_values_constantF)]
    ).grid(row=1, column=0, sticky='W')

    check_CSV = tk.Checkbutton(
        check_box,
        text="CSV (F(pN) | d)",
        variable=check_box_CSV,
        command=lambda: [select_box(check_box_CSV, check_box_HF, check_box_LF), parameters(default_values_CSV, default_values_FIT, default_values_constantF)]
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

    downsample_value = tk.StringVar()
    downsample_value1 = tk.Entry(parameter_frame, textvariable=downsample_value)

    Filter_degree = tk.StringVar()
    Filter_degree1 = tk.Entry(parameter_frame, textvariable=Filter_degree)

    Filter_cut_off = tk.StringVar()
    Filter_cut_off1 = tk.Entry(parameter_frame, textvariable=Filter_cut_off)

    Force_Min = tk.StringVar()
    Force_Min1 = tk.Entry(parameter_frame, textvariable=Force_Min)

    Z_score_force = tk.StringVar()
    Z_score_force1 = tk.Entry(parameter_frame, textvariable=Z_score_force)

    Z_score_distance = tk.StringVar()
    Z_score_distance1 = tk.Entry(parameter_frame, textvariable=Z_score_distance)

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
    frame1.grid(row=0, column=0, sticky='N')
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

    # parameters that occur double (tab1 and tab4)
    downsample_value2 = tk.Entry(frame1, textvariable=downsample_value)
    Filter_degree2 = tk.Entry(frame1, textvariable=Filter_degree)
    Filter_cut_off2 = tk.Entry(frame1, textvariable=Filter_cut_off)
    Force_Min2 = tk.Entry(frame1, textvariable=Force_Min)
    Z_score_force2 = tk.Entry(frame1, textvariable=Z_score_force)
    Z_score_distance2 = tk.Entry(frame1, textvariable=Z_score_distance)

    # parameters only in advanced settings
    step_d_variable = tk.StringVar()
    step_d_value = tk.Entry(frame1, textvariable=step_d_variable)

    window_size_variable = tk.StringVar()
    window_size_value = tk.Entry(frame1, textvariable=window_size_variable)

    STD_difference_variable = tk.StringVar()
    STD_difference_value = tk.Entry(frame1, textvariable=STD_difference_variable)

    Frequency_variable = tk.StringVar()
    Frequency_value = tk.Entry(frame1, textvariable=Frequency_variable)

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
    # Labels
    Cluster_fitting = tk.Label(frame3, text='FITTING', font='Helvetica 9 bold')
    check_box_WLC = tk.IntVar(value=1)
    check_box_FJC = tk.IntVar(value=0)
    Label_dsLp = tk.Label(frame3, text='dsLp [nm]')
    Label_dsLp_up = tk.Label(frame3, text='dsLp upper bound [nm]')
    Label_dsLp_low = tk.Label(frame3, text='dsLp lower bound [nm]')
    Label_dsLc = tk.Label(frame3, text='dsLc [nm]')
    Label_ssLp = tk.Label(frame3, text='ssLp [nm]')
    Label_ssLc = tk.Label(frame3, text='ssLc [nm]')
    Label_stiffness_ds = tk.Label(frame3, text='dsK0 [pN]')
    Label_stiffness_ds_up = tk.Label(frame3, text='dsK0 upper bound [pN]')
    Label_stiffness_ds_low = tk.Label(frame3, text='dsK0 lower bound [pN]')
    Label_stiffness_ss = tk.Label(frame3, text='ssK0 [pN]')
    Label_f_offset = tk.Label(frame3, text='Force offset [pN]')
    Label_f_offset_up = tk.Label(frame3, text='Force offset upper bound [pN]')
    Label_f_offset_low = tk.Label(frame3, text='Force offset lower bound [pN]')
    Label_d_offset = tk.Label(frame3, text='Distance offset [nm]')
    Label_d_offset_up = tk.Label(frame3, text='Distance offset upper bound [nm]')
    Label_d_offset_low = tk.Label(frame3, text='Distance offset lower bound [nm]')

    # Entry widgets
    dsLp_variable = tk.StringVar()
    dsLp = tk.Entry(frame3, textvariable=dsLp_variable)
    dsLp_up_variable = tk.StringVar()
    dsLp_up = tk.Entry(frame3, textvariable=dsLp_up_variable)
    dsLp_low_variable = tk.StringVar()
    dsLp_low = tk.Entry(frame3, textvariable=dsLp_low_variable)
    dsLc_variable = tk.StringVar()
    dsLc = tk.Entry(frame3, textvariable=dsLc_variable)
    ssLp_variable = tk.StringVar()
    ssLp = tk.Entry(frame3, textvariable=ssLp_variable)
    ssLc_variable = tk.StringVar()
    ssLc = tk.Entry(frame3, textvariable=ssLc_variable)
    stiff_ds_variable = tk.StringVar()
    stiff_ds = tk.Entry(frame3, textvariable=stiff_ds_variable)
    stiff_ds_up_variable = tk.StringVar()
    stiff_ds_up = tk.Entry(frame3, textvariable=stiff_ds_up_variable)
    stiff_ds_low_variable = tk.StringVar()
    stiff_ds_low = tk.Entry(frame3, textvariable=stiff_ds_low_variable)
    stiff_ss_variable = tk.StringVar()
    stiff_ss = tk.Entry(frame3, textvariable=stiff_ss_variable)
    f_off_variable = tk.StringVar()
    f_off = tk.Entry(frame3, textvariable=f_off_variable)
    f_off_up_variable = tk.StringVar()
    f_off_up = tk.Entry(frame3, textvariable=f_off_up_variable)
    f_off_low_variable = tk.StringVar()
    f_off_low = tk.Entry(frame3, textvariable=f_off_low_variable)
    d_off_variable = tk.StringVar()
    d_off = tk.Entry(frame3, textvariable=d_off_variable)
    d_off_up_variable = tk.StringVar()
    d_off_up = tk.Entry(frame3, textvariable=d_off_up_variable)
    d_off_low_variable = tk.StringVar()
    d_off_low = tk.Entry(frame3, textvariable=d_off_low_variable)

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
    TOMATO_figure_frame = tk.Canvas(tab5, height=650, width=1000, borderwidth=1, relief='ridge')
    TOMATO_figure_frame.grid(row=0, column=0, rowspan=2)

    TOMATO_button_frame = tk.Frame(tab5)
    TOMATO_button_frame.grid(row=0, column=1, sticky=tk.N + tk.W)

    TOMATO_parameter_frame = tk.Frame(tab5, borderwidth=1, relief='ridge')
    TOMATO_parameter_frame.grid(row=1, column=1, sticky=tk.N + tk.W)

    TOMATO_frame_table = tk.Frame(tab5, width=400, height=100)
    TOMATO_frame_table.grid(row=2, column=0)

    ### create entry widgets ###
    # shift in distance
    label_shift_d = tk.Label(TOMATO_parameter_frame, text='Distance offset [nm]')
    label_shift_d.grid(row=0, column=0, sticky=tk.E + tk.W, padx=(4, 2), pady=2)
    entry_shift_d = tk.Entry(TOMATO_parameter_frame, textvariable=d_off_variable).grid(row=0, column=1)

    # shift in F
    label_shift_F = tk.Label(TOMATO_parameter_frame, text='Force offset [pN]')
    label_shift_F.grid(row=0, column=2, sticky=tk.E + tk.W, padx=(4, 2), pady=2)
    entry_shift_F = tk.Entry(TOMATO_parameter_frame, textvariable=f_off_variable).grid(row=0, column=3)

    # K0 for both
    # ds
    label_ds_St = tk.Label(TOMATO_parameter_frame, text='dsK0 [pN]')
    label_ds_St.grid(row=1, column=0, sticky=tk.E + tk.W, padx=(4, 2), pady=2)
    entry_ds_St = tk.Entry(TOMATO_parameter_frame, textvariable=stiff_ds_variable).grid(row=1, column=1)

    # ss
    label_ss_St = tk.Label(TOMATO_parameter_frame, text='ssK0 [pN]')
    label_ss_St.grid(row=1, column=2, sticky=tk.E + tk.W, padx=(4, 2), pady=2)
    entry_ss_St = tk.Entry(TOMATO_parameter_frame, textvariable=stiff_ss_variable).grid(row=1, column=3)

    ## ds handle part
    # ds handle persistance length
    label_ds_Lp = tk.Label(TOMATO_parameter_frame, text='dsLp [nm]')
    label_ds_Lp.grid(row=2, column=0, sticky=tk.E + tk.W, padx=(4, 2), pady=2)
    entry_ds_Lp = tk.Entry(TOMATO_parameter_frame, textvariable=dsLp_variable).grid(row=2, column=1)

    # ds handle  contour length
    label_ds_Lc = tk.Label(TOMATO_parameter_frame, text='dsLc [nm]')
    label_ds_Lc.grid(row=3, column=0, sticky=tk.E + tk.W, padx=(4, 2), pady=2)
    entry_ds_Lc = tk.Entry(TOMATO_parameter_frame, textvariable=dsLc_variable).grid(row=3, column=1)

    ## ss RNA part
    # ss RNA persistance length
    label_ss_Lp = tk.Label(TOMATO_parameter_frame, text=' ssLp [nm]')
    label_ss_Lp.grid(row=2, column=2, sticky=tk.E + tk.W, padx=(4, 2), pady=2)
    entry_ss_Lp = tk.Entry(TOMATO_parameter_frame, textvariable=ssLp_variable).grid(row=2, column=3)

    # ss RNA contour length
    label_ss_Lc = tk.Label(TOMATO_parameter_frame, text=' ssLc [nm]')
    label_ss_Lc.grid(row=3, column=2, sticky=tk.E + tk.W, padx=(4, 2), pady=2)
    entry_ss_Lc = tk.Entry(TOMATO_parameter_frame, textvariable=ssLc_variable).grid(row=3, column=3)

    # start position
    label_strF = tk.Label(TOMATO_parameter_frame, text='Force [pN]')
    label_strF.grid(row=4, column=1, sticky=tk.E + tk.W, pady=(25, 0))

    label_strD = tk.Label(TOMATO_parameter_frame, text='Distance [nm]')
    label_strD.grid(row=4, column=2, sticky=tk.E + tk.W, pady=(25, 0))

    entryText_startF = tk.StringVar()
    entry_startF = tk.Entry(TOMATO_parameter_frame, textvariable=entryText_startF).grid(row=5, column=1, pady=2)

    entryText_startD = tk.StringVar()
    entry_startD = tk.Entry(TOMATO_parameter_frame, textvariable=entryText_startD).grid(row=5, column=2, pady=2)

    # end position

    entryText_endF = tk.StringVar()
    entry_endF = tk.Entry(TOMATO_parameter_frame, textvariable=entryText_endF).grid(row=6, column=1, pady=2)

    entryText_endD = tk.StringVar()
    entry_endD = tk.Entry(TOMATO_parameter_frame, textvariable=entryText_endD).grid(row=6, column=2, pady=2)

    entryText_filename = tk.StringVar()
    entry_filename = tk.Entry(TOMATO_frame_table, textvariable=entryText_filename, width=100)
    entry_filename.grid(row=0, column=0)

    # create button widgets that use the defined functions
    browseButton_folder = tk.Button(TOMATO_button_frame, text="      Choose folder     ", command=open_folder, bg='green', fg='white', font=('Arial', 11, 'bold'))
    browseButton_folder.grid(row=0, column=0, padx=4, pady=4)

    button_save = tk.Button(TOMATO_button_frame, text='Save results table', command=export_table, bg='palegreen2', font=('Arial', 11, 'bold'))
    button_save.grid(row=1, column=0, padx=4, pady=4)

    button_reset_parameters = tk.Button(TOMATO_button_frame, text='Reset parameters', command=lambda: parameters(0, default_values_FIT, default_values_constantF) if check_box_HF == 1 else (parameters(default_values_HF, default_values_FIT, default_values_constantF) if check_box_LF == 1 else parameters(default_values_CSV, default_values_FIT, default_values_constantF)), bg='palegreen2', font=('Arial', 11, 'bold'))
    button_reset_parameters.grid(row=1, column=1, padx=4, pady=4)

    label_info = tk.Label(TOMATO_button_frame, text='TOMATO is unresponsive during fitting, please be patient\n\nmark step start <s>\n mark step end <e>\n save marked step <Ctrl+s>\n start analysis <f>\n delete results line <mark+del>\n next curve <Right arrow> or <d>\n previous curve <Left arrow> or <a>')
    label_info.grid(row=2, column=0, padx=4, pady=4)

    button_start = tk.Button(TOMATO_parameter_frame, text='Set start', command=start_click, bg='lightsteelblue2', font=('Arial', 10, 'bold'))
    button_start.grid(row=5, column=0, pady=2)

    button_end = tk.Button(TOMATO_parameter_frame, text='Set end', command=end_click, bg='lightsteelblue2', font=('Arial', 10, 'bold'))
    button_end.grid(row=6, column=0, pady=2)

    step_number = 1
    button_save_step = tk.Button(TOMATO_parameter_frame, text='Save step', command=save_step, bg='lightsteelblue2', font=('Arial', 10, 'bold'))
    button_save_step.grid(row=7, column=0, pady=2)

    button_delete_step = tk.Button(TOMATO_parameter_frame, text='Delete step', command=delete_step, bg='lightsteelblue2', font=('Arial', 10, 'bold'))
    button_delete_step.grid(row=8, column=0, pady=2)

    button_start_analysis = tk.Button(TOMATO_parameter_frame, text='Analyze curve', command=analyze_steps, bg='#df4c4c', font=('Arial', 10, 'bold'))
    button_start_analysis.grid(row=7, column=1, rowspan=2, pady=2)

    ## show the fitting parameters in a table
    # create Treeview for results table
    cols = (
        'Filename',
        'step number',
        'Force step start [pN]',
        'Force step end [pN]',
        'delta force [pN]',
        'extension step start [nm]',
        'extension step end [nm]',
        'Step length [nm]',
        'ds contour length',
        'ds persistance Length',
        'ds stiffness (K0) [pN]',
        'ss contour Length',
        'ss persistance Length',
        'ss stiffness (K0) [pN]',
        'Force offset',
        'Distance offset',
        'Work [pN*nm]',
        'Work [kT]'
    )

    tree_results = ttk.Treeview(TOMATO_frame_table, columns=cols, show='headings', height=5)
    # set column headings
    for col in cols:
        tree_results.heading(col, text=col)
        tree_results.column(col, minwidth=25, width=65)
    tree_results.grid(row=1, column=0, padx=5, pady=5)

    # create Treeview for the steps to analyze
    cols_steps = ('Step number', 'F start', 'F end', 'Step start', 'Step end')
    tree_steps = ttk.Treeview(TOMATO_parameter_frame, columns=cols_steps, show='headings', height=5)
    # set column headings
    for col in cols_steps:
        tree_steps.heading(col, text=col)
        tree_steps.column(col, minwidth=25, width=65)
    tree_steps.grid(row=9, column=0, columnspan=2, pady=5)
    ######### TOMATO end ############

    ############ POTATO last part ###############
    # put default values into the widgets
    parameters(default_values_HF, default_values_FIT, default_values_constantF)

    # loop ensuring the GUI is running until closed
    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()
