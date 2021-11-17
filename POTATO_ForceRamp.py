"""Copyright 2021 Helmholtz-Zentrum für Infektionsforschung GmbH"""

import pandas as pd
import h5py
import numpy as np
from pathlib import Path

# relative imports
from POTATO_fitting import fitting_ds, fitting_ss, plot_fit
from POTATO_preprocessing import preprocess_RAW, trim_data, create_derivation
from POTATO_find_steps import find_steps_F, find_steps_PD, find_common_steps, calc_integral, save_figure
"""define the functions of the subprocess processing the data"""


# open a folder containing raw data and lead through the analysis process
def start_subprocess(analysis_folder, timestamp, Files, input_settings, input_format, export_data, input_fitting, output_q):
    # empty dataframe to store all step results of all curves in the folder
    total_results_steps = pd.DataFrame()

    # create dataframe to store all fitting parameters of all curves in the folder
    header_fit = [
        "filename",
        "model",
        "log_likelihood",
        "Lc_ds",
        "Lp_ds",
        "Lp_ds_stderr",
        "St_ds",
        "Lc_ss",
        "Lc_ss_stderr",
        "Lp_ss",
        "St_ss",
        "f_offset",
        "d_offset",
        "Work_(pN*nm)",
        "Work_(kB*T)"
    ]

    total_results_fit = pd.DataFrame(columns=header_fit)

    # iterate through the files in the selected folder
    i = 0
    # proceed differently with h5 and csv files
    while i < len(Files):
        if input_format['CSV'] == 1:
            df = pd.read_csv(Files[i])
            directory_i = Path(Files[i])
            filename_i = directory_i.name[:-4]
            # access the raw data
            Force_1x = df.to_numpy()[:, 0]
            Distance_1x = df.to_numpy()[:, 1]
            # accessing the data frequency from user input
            Frequency_value = input_settings['data_frequency']
            Force_Distance, Force_Distance_um = preprocess_RAW(Force_1x, Distance_1x, input_settings)
        else:
            with h5py.File(Files[i], "r") as f:
                directory_i = Path(Files[i])
                filename_i = directory_i.name[:-3]

                # access the raw data
                if input_format['HF'] == 1:
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
                    timestamp_F_LF = load_force.attrs['Start time (ns)']

                    Frequency_value = size_F_LF / ((stop_time_F_LF - timestamp_F_LF) / 10**9)
                    Force_Distance, Force_Distance_um = preprocess_RAW(Force_1x, Distance_1x, input_settings)

        # Export down sampled and smoothened FD values
        if export_data['export_SMOOTH'] == 1:
            filename = analysis_folder + "/" + filename_i + "_smooth_" + timestamp + ".csv"
            np.savetxt(filename, Force_Distance_um, delimiter=",")
        else:
            pass

        # trim data below specified force thresholds
        F_trimmed, PD_trimmed, F_low = trim_data(Force_Distance, input_settings['F_min'])

        # create force and distance derivation of the pre-processed data to be able to identify steps
        derivation_array = create_derivation(input_settings, Frequency_value)

        """find steps based on force derivation"""
        filename_results = analysis_folder + "/" + filename_i + "_results_" + timestamp + ".csv"

        try:
            results_F, PD_start_F = find_steps_F(
                input_settings,
                filename_i,
                Force_Distance,
                derivation_array
            )

            results_F_list = list(results_F)

            if export_data['export_STEPS'] == 1:
                steps_results_F = pd.DataFrame(results_F_list)
                with open(filename_results, 'a+') as f:
                    f.write('\nSteps found by force derivation:\n')
                steps_results_F.to_csv(filename_results, mode='a', index=False, header=True)
            else:
                pass

        except:
            results_F = []
            PD_start_F = []
            print("Error in finding steps for file " + str(filename_i) + '\n' 'There was an error in finding Force steps')
            pass

        """find steps based on distance derivation"""

        try:
            results_PD, PD_start_PD = find_steps_PD(
                input_settings,
                filename_i,
                Force_Distance,
                derivation_array
            )

            results_PD_list = list(results_PD)

            if export_data['export_STEPS'] == 1:
                steps_results_PD = pd.DataFrame(results_PD_list)
                with open(filename_results, 'a+') as f:
                    f.write('\nSteps found by distance derivation:\n')
                steps_results_PD.to_csv(filename_results, mode='a', index=False, header=True)

        except:
            results_PD = []
            PD_start_PD = []
            err_PD = str("Error in finding steps for file " + str(filename_i) + '\n' 'There was an error in finding Distance steps')
            print(err_PD)
            pass

        # save plot with FD-curve, derivations and found steps
        save_figure(
            export_data['export_PLOT'],
            timestamp,
            filename_i,
            analysis_folder,
            Force_Distance,
            derivation_array,
            F_trimmed,
            PD_trimmed,
            PD_start_F,
            PD_start_PD
        )

        # when steps are found by force AND distance derivation, they are considered common steps
        try:
            common_steps = find_common_steps(results_F_list, results_PD_list)
            # to match with the fitting rows (always one more than steps) put a 'step 0' as first line
            common_steps_results = [{'filename': filename_i, 'Derivation of': '', 'step #': 0, 'F1': '', 'F2': '', 'Fc': '', 'step start': '', 'step end': '', 'step length': ''}]
        except:
            err_FCS = str("Error in finding common steps" + str(filename_i) + '\n' 'There was an error in finding common steps')
            output_q.put(err_FCS)
            pass

        # append common steps to the 'step 0'
        if common_steps:
            for x in range(len(common_steps)):
                common_steps_results.append(common_steps[x])

            # convert common steps to dataframe for export
            common_steps_results = pd.DataFrame(common_steps_results)

            # export the steps into the results for ONLY this file
            with open(filename_results, 'a+') as f:
                f.write('\nCommon steps:\n')
            common_steps_results.to_csv(filename_results, mode='a', index=False, header=True)

            # put common steps into a total_results dataframe so all steps from all files of the analysed folder can be exported together
            total_results_steps = total_results_steps.append(common_steps_results, ignore_index=True, sort=False)

        else:
            common_steps_results = [{'filename': filename_i, 'Derivation of': '', 'step #': 'no common steps', 'F1': '', 'F2': '', 'Fc': '', 'step start': '', 'step end': '', 'step length': ''}]
            total_results_steps = total_results_steps.append(common_steps_results, ignore_index=True, sort=False)

        '''if common steps were found, try to fit FD-Curve'''
        empty = {
            'filename': filename_i,
            'model': 'None',
            'log_likelihood': 'None',
            'Lc_ds': 'None',
            'Lp_ds': 'None',
            'Lp_ds_stderr': 'None',
            'St_ds': 'None',
            'f_offset': 'None',
            'd_offset': 'None'
        }

        if export_data['export_FIT'] == 1:
            try:
                export_fit = []
                fit = []
                start_force_ss = []
                start_distance_ss = []
                integral_ss_fit_start = []
                integral_ss_fit_end = []

                # try to fit all parts of curve based on the common steps
                try:
                    # fit part between start of the FD-cure up to the first common step
                    export_fit_ds, area_ds = fitting_ds(
                        filename_i,
                        input_settings,
                        export_data,
                        input_fitting,
                        float(common_steps[0]['step start']),
                        Force_Distance,
                        derivation_array,
                        F_low
                    )

                    export_fit.append(export_fit_ds)

                    # fit parts after steps, when more than one common step was found, there are multiple parts to fit
                    if len(common_steps) > 1:
                        for n in range(0, len(common_steps) - 1):
                            # try to fit each part of the curve, if one of the parts can not be fitted, still try to fit the others
                            try:
                                fit_ss, f_fitting_region_ss, d_fitting_region_ss, export_fit_ss, area_ss_fit_start, area_ss_fit_end = fitting_ss(
                                    filename_i,
                                    input_settings,
                                    export_data,
                                    input_fitting,
                                    float(common_steps[n]['step end']),
                                    float(common_steps[n + 1]['step start']),
                                    Force_Distance, 1, 1,
                                    derivation_array,
                                    F_low
                                )

                                fit.append(fit_ss)
                                start_force_ss.append(f_fitting_region_ss)
                                start_distance_ss.append(d_fitting_region_ss)
                                export_fit.append(export_fit_ss)
                                integral_ss_fit_start.append(area_ss_fit_start)
                                integral_ss_fit_end.append(area_ss_fit_end)

                            except:
                                export_fit.append(empty)
                                pass

                    # fit the last part of the curve
                    try:
                        fit_ss, f_fitting_region_ss, d_fitting_region_ss, export_fit_ss, area_ss_fit_start, area_ss_fit_end = fitting_ss(
                            filename_i,
                            input_settings,
                            export_data,
                            input_fitting,
                            float(common_steps[len(common_steps) - 1]['step end']),
                            max(derivation_array[:, 1]),
                            Force_Distance, 1, 1,
                            derivation_array,
                            F_low
                        )

                        fit.append(fit_ss)
                        start_force_ss.append(f_fitting_region_ss)
                        start_distance_ss.append(d_fitting_region_ss)
                        export_fit.append(export_fit_ss)
                        integral_ss_fit_start.append(area_ss_fit_start)
                        integral_ss_fit_end.append(area_ss_fit_end)

                    except:
                        export_fit.append(empty)
                        pass

                    '''from the fits, work put into the system is calculated'''
                    if common_steps:
                        work_per_step = [0]  # in pN*nm
                        kT_per_step = [0]    # in kT

                        work_first_step, kT_1 = calc_integral(
                            area_ds,
                            integral_ss_fit_start[0],
                            common_steps[0]['step start'],
                            common_steps[0]['step end'],
                            common_steps[0]['F1'],
                            common_steps[0]['F2']
                        )

                        print("Work of first step: " + str(work_first_step))
                        work_per_step.append(work_first_step)
                        kT_per_step.append(kT_1)

                        if len(common_steps) > 1:
                            for n in range(0, len(common_steps) - 1):
                                work_step_n, kT_n = calc_integral(
                                    integral_ss_fit_end[n],
                                    integral_ss_fit_start[n + 1],
                                    common_steps[n + 1]['step start'],
                                    common_steps[n + 1]['step end'],
                                    common_steps[n + 1]['F1'],
                                    common_steps[n + 1]['F2']
                                )

                                work_per_step.append(work_step_n)
                                kT_per_step.append(kT_n)

                        j = 0
                        for dict in export_fit:
                            dict["Work_(pN*nm)"] = work_per_step[j]
                            dict["Work_(kB*T)"] = kT_per_step[j]
                            j += 1

                # if no step was found, the common step index 0 is not available and will raise an IndexError.
                # So in this case the fit will be performed for the whole curve from beginning to end.
                except IndexError:
                    if not common_steps:
                        export_fit_ds, area_ds = fitting_ds(
                            filename_i,
                            input_settings,
                            export_data,
                            input_fitting,
                            derivation_array[-1, 1],
                            Force_Distance,
                            derivation_array,
                            F_low
                        )

                        export_fit.append(export_fit_ds)

                total_results_fit = total_results_fit.append(export_fit, ignore_index=True, sort=False)

                # create a plot for the fitted curve
                plot_fit(fit, start_force_ss, start_distance_ss, Force_Distance, analysis_folder, filename_i, timestamp)

            except:
                print('Something went wrong with fitting')
                pass

        if i == 0:
            print('\nHard work ahead!\n')
            output_q.put('Hard work ahead!')
        elif i == int(len(Files) / 2):
            print('\nHalf way there!\n')
            output_q.put('Half way there!')
            print()
        elif i == len(Files) - 1:
            print('\nAlmost there!\n')
            output_q.put('Almost there!')
        elif i == len(Files):
            print('Analysis finished! \nProgram can be closed.')
            output_q.put('Analysis finished! \nProgram can be closed.')

        i = i + 1
        print('done', i, 'from', len(Files))
        out_progress = str('Done ' + str(i) + ' from ' + str(len(Files)))
        output_q.put(out_progress)

        print(filename_i)
        output_q.put(filename_i)

    '''after folder analysis is done, export total results (all steps + fit parameters) in one file'''
    if export_data['export_TOTAL'] == 1:
        filename_total_results = analysis_folder + '/total_results_' + timestamp + '.csv'

        with open(filename_total_results, 'w') as f:
            f.write('Common steps from all curves of the folder:\n')

        results_total_total = pd.concat([total_results_steps, total_results_fit], axis=1)
        results_total_total.to_csv((filename_total_results), mode='a', index=False)
    else:
        pass
