""" POTATO -- 2021-10-14 -- Version 0.1
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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
from PIL import ImageTk, Image
import pandas as pd
import os
import h5py
import glob
import time
import multiprocessing as mp
import json

# relative imports
from POTATO_ForceRamp import start_subprocess
from POTATO_preprocessing import preprocess_RAW
from POTATO_config import default_values_HF, default_values_LF, default_values_CSV, default_values_FIT, default_values_constantF
from POTATO_constantF import get_constantF, display_constantF, fit_constantF

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
    global image_number
    # check user input
    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
    image_number = []
    # ask wich directory should be analysed
    folder = filedialog.askdirectory()
    root.title('POTATO -- ' + str(folder))

    # decide which input format was choosen
    if input_format['CSV'] == 1:
        folder_path = str(folder + "/*.csv")
    else:
        folder_path = str(folder + "/*.h5")

    Files = glob.glob(folder_path)
    print('files to analyse', len(Files))

    # print number of files to analyse, if no files found give an error
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
    STD1_threshold1.delete(0, "end")
    STD1_threshold1.insert("end", default_values['Z-score'])
    STD1_threshold2.delete(0, "end")
    STD1_threshold2.insert("end", default_values['Z-score'])
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
        'z-score': float(STD1_threshold2.get()),
        'window_size': int(window_size_value.get()),
        'data_frequency': float(Frequency_value.get()),
        'STD_diff': float(STD_difference_value.get())
    }

    input_format = {
        'HF': check_box_HF.get(),
        'LF': check_box_LF.get(),
        'CSV': check_box_CSV.get()
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
            Force_1x = raw_data.get("Force HF/Force 1x")
            Distance_1x = raw_data.get("Distance/Piezo Distance")

        elif input_format['LF'] == 1:
            load_force = raw_data.get("Force LF/Force 1x")
            Force_1x = load_force[:]['Value'][:]
            load_distance = raw_data.get("Distance/Distance 1")[:]
            Distance_1x = load_distance['Value'][:]
            print(Force_1x, Distance_1x)
        FD, FD_um = preprocess_RAW(Force_1x, Distance_1x, input_settings)
        display_RAW_FD(FD[:, 0], FD[:, 1], Force_1x[::input_settings['downsample_value']], Distance_1x[::input_settings['downsample_value']] * 1000)


# display a single csv file (tab2)
def getRAW_File_csv():
    if not check_box_CSV.get() == 1:
        check_box_CSV.set(value=1)
        select_box(check_box_CSV, check_box_HF, check_box_LF)
        parameters(parameter_frame, default_values_CSV, default_values_FIT, default_values_constantF)
        parameters(tab3, default_values_CSV, default_values_FIT, default_values_constantF)
    else:
        pass

    input_settings, input_format, export_data, input_fitting, input_constantF = check_settings()
    import_file_path = filedialog.askopenfilename()

    df = pd.read_csv(import_file_path)

    Force_HF_1x = df.to_numpy()[:, 0]
    Distance_HF_1x = df.to_numpy()[:, 1]

    FD, FD_um = preprocess_RAW(Force_HF_1x, Distance_HF_1x, input_settings)
    display_RAW_FD(FD[:, 0], FD[:, 1], Force_HF_1x[::input_settings['downsample_value']], Distance_HF_1x[::input_settings['downsample_value']] * 1000)


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
    if messagebox.askokcancel("Quit", "Do you really want to quit?"):
        root.quit()


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

    tab1.grid(row=0, column=0, padx=2, pady=2)
    tab2.grid(row=0, column=0, padx=2, pady=2)
    tab3.grid(row=0, column=0, padx=2, pady=2)
    tab4.grid(row=0, column=0, padx=2, pady=2)

    # ATTENTION - tab3 and tab4 are displayed the other way round in the GUI
    tabControl.add(tab1, text="Folder Analysis")
    tabControl.add(tab2, text="Show Single File")
    tabControl.add(tab4, text="Constant Force Analysis")
    tabControl.add(tab3, text="Advanced Settings")

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

    def select_box(check_box_1, check_box_2, check_box_3):
        if check_box_1.get() == 1:
            check_box_2.set(value=0)
            check_box_3.set(value=0)
        elif check_box_1.get() == 0 and check_box_2.get() == 0 and check_box_3.get() == 0:
            check_box_1.set(value=1)

    check_box_HF = tk.IntVar(value=1)
    check_box_LF = tk.IntVar()
    check_box_CSV = tk.IntVar()

    check1 = tk.Checkbutton(
        check_box,
        text="High Frequency (Piezo Distance)",
        variable=check_box_HF,
        command=lambda: [select_box(check_box_HF, check_box_LF, check_box_CSV), parameters(parameter_frame, default_values_HF, default_values_FIT, default_values_constantF)]
    ).grid(row=0, column=0, sticky='W')

    check2 = tk.Checkbutton(
        check_box,
        text="Low Frequency",
        variable=check_box_LF,
        command=lambda: [select_box(check_box_LF, check_box_HF, check_box_CSV), parameters(parameter_frame, default_values_LF, default_values_FIT, default_values_constantF)]
    ).grid(row=1, column=0, sticky='W')

    check3 = tk.Checkbutton(
        check_box,
        text="CSV (F/D)",
        variable=check_box_CSV,
        command=lambda: [select_box(check_box_CSV, check_box_HF, check_box_LF), parameters(parameter_frame, default_values_CSV, default_values_FIT, default_values_constantF)]
    ).grid(row=2, column=0, sticky='W')

    figure_frame = tk.Canvas(tab1, height=650, width=1000, borderwidth=1, relief='ridge')
    figure_frame.grid(row=1, column=0)

    parameter_frame = tk.Frame(tab1)
    parameter_frame.grid(row=1, column=1, sticky='NE')

    """ parameter frame """
    Cluster_preprocessing = tk.Label(parameter_frame, text='PREPROCESSING', font='Helvetica 9 bold')
    Label_downsample = tk.Label(parameter_frame, text='Downsampling rate')
    Label_Filter1 = tk.Label(parameter_frame, text='Butterworth filter degree')
    Label_Filter2 = tk.Label(parameter_frame, text='Cut-off frequency')
    Label_ForceMin = tk.Label(parameter_frame, text='Force threshold, pN')
    Cluster_statistics = tk.Label(parameter_frame, text='STATISTICS', font='Helvetica 9 bold')
    Label_STD_1 = tk.Label(parameter_frame, text='Z-score')

    downsample_value1 = tk.Entry(parameter_frame)
    downsample_value1.bind("<Return>", lambda event: user_input(event, downsample_value1, downsample_value2))

    Filter_degree1 = tk.Entry(parameter_frame)
    Filter_degree1.bind("<Return>", lambda event: user_input(event, Filter_degree1, Filter_degree2))

    Filter_cut_off1 = tk.Entry(parameter_frame)
    Filter_cut_off1.bind("<Return>", lambda event: user_input(event, Filter_cut_off1, Filter_cut_off2))

    Force_Min1 = tk.Entry(parameter_frame)
    Force_Min1.bind("<Return>", lambda event: user_input(event, Force_Min1, Force_Min2))

    STD1_threshold1 = tk.Entry(parameter_frame)
    STD1_threshold1.bind("<Return>", lambda event: user_input(event, STD1_threshold1, STD1_threshold2))

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
    Label_STD_1.grid(row=6, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    STD1_threshold1.grid(row=6, column=1, padx=2, pady=2)

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
    Cluster_derivation = tk.Label(frame1, text="DERIVATION", font='Helvetica 9 bold')
    Label_step_d = tk.Label(frame1, text='Step d')
    Label_Frequency = tk.Label(frame1, text='Data frequency, Hz')
    Cluster_statistics = tk.Label(frame1, text='STATISTICS', font='Helvetica 9 bold')
    Label_STD_1 = tk.Label(frame1, text='Z-score')
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

    STD1_threshold2 = tk.Entry(frame1)
    STD1_threshold2.bind("<Return>", lambda event: user_input(event, STD1_threshold2, STD1_threshold1))

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

    Cluster_derivation.grid(row=5, column=0, padx=2, pady=(20, 2))
    Label_step_d.grid(row=6, column=0, padx=2, pady=2)
    step_d_value.grid(row=6, column=1, padx=(0, 20), pady=2)

    Label_Frequency.grid(row=7, column=0, padx=2, pady=2)
    Frequency_value.grid(row=7, column=1, padx=(0, 20), pady=2)

    Cluster_statistics.grid(row=8, column=0, padx=2, pady=(20, 2))
    Label_STD_1.grid(row=9, column=0, sticky=tk.E + tk.W, padx=2, pady=2)
    STD1_threshold2.grid(row=9, column=1, padx=(0, 20), pady=2)

    Label_window_size.grid(row=11, column=0, padx=2, pady=2)
    window_size_value.grid(row=11, column=1, padx=(0, 20), pady=2)

    Label_STD_difference.grid(row=12, column=0, padx=2, pady=2)
    STD_difference_value.grid(row=12, column=1, padx=(0, 20), pady=2)

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

    # put default values into the widgets
    parameters(parameter_frame, default_values_HF, default_values_FIT, default_values_constantF)

    # loop ensuring the GUI is running until closed
    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.mainloop()
