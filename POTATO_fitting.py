import matplotlib.pyplot as plt
import lumicks.pylake as lk
import numpy as np
from scipy.integrate import simps
from matplotlib.lines import Line2D


"""define the functions used for fitting"""


def fitting_ds(filename_i, input_settings, export_data, input_fitting, i_start, Force_Distance, derivation_array, F_low):
    global model_ds, fit_ds
    global ds_fit_dict
    global f_fitting_region_ds, d_fitting_region_ds
    global export_fit_ds
    global fitting_model

    start_step1 = np.where(derivation_array[:, 1] == i_start)
    start_step1 = start_step1[0][0]

    f_fitting_region_ds = Force_Distance[0:start_step1 * input_settings['step_d'] + len(F_low), 0]
    d_fitting_region_ds = Force_Distance[0:start_step1 * input_settings['step_d'] + len(F_low), 1]

    model_ds = lk.inverted_odijk("ds_part").subtract_independent_offset() + lk.force_offset("ds_part")

    fit_ds = lk.FdFit(model_ds)

    fit_ds.add_data("Double stranded", f_fitting_region_ds, d_fitting_region_ds)
    # Persistance length bounds
    fit_ds["ds_part/Lp"].value = input_fitting['lp_ds']
    fit_ds["ds_part/Lp"].upper_bound = input_fitting['lp_ds_up']
    fit_ds["ds_part/Lp"].lower_bound = input_fitting['lp_ds_low']

    # Force shift bounds
    fit_ds["ds_part/f_offset"].value = input_fitting['offset_f']
    fit_ds["ds_part/f_offset"].upper_bound = input_fitting['offset_f_up']
    fit_ds["ds_part/f_offset"].lower_bound = input_fitting['offset_f_low']

    # distance shift bounds
    fit_ds["ds_part/d_offset"].value = input_fitting['offset_d']
    fit_ds["ds_part/d_offset"].upper_bound = input_fitting['offset_d_up']
    fit_ds["ds_part/d_offset"].lower_bound = input_fitting['offset_d_low']

    # stiffnes
    fit_ds["ds_part/St"].value = input_fitting['ds_stiff']
    fit_ds["ds_part/St"].upper_bound = input_fitting['ds_stiff_up']
    fit_ds["ds_part/St"].lower_bound = input_fitting['ds_stiff_low']

    # contour length
    Lc_initial_guess = input_fitting['lc_ds']  # nm
    Lc_range = 5
    fit_ds["ds_part/Lc"].upper_bound = Lc_initial_guess + Lc_range
    fit_ds["ds_part/Lc"].lower_bound = Lc_initial_guess - Lc_range
    fit_ds["ds_part/Lc"].value = Lc_initial_guess
    fit_ds["ds_part/Lc"].unit = 'nm'

    fit_ds.fit()
    fit_qual = fit_ds.log_likelihood()
    print(fit_ds.params)

    # calculate the integral until the first unfolding step
    # used to calculate the work done by the machine
    distance_integral = np.arange(min(Force_Distance[:, 1]), i_start)
    ds_integral = model_ds(distance_integral, fit_ds)
    area_ds = simps(ds_integral)
    print("area_ds = " + str(area_ds))

    # export the fitting parameters
    ds_fit_dict = {
        'filename': filename_i,
        'model': 'WLC',
        'log_likelihood': fit_qual,
        'Lc_ds': fit_ds["ds_part/Lc"].value,
        'Lp_ds': fit_ds["ds_part/Lp"].value,
        'Lp_ds_stderr': fit_ds["ds_part/Lp"].stderr,
        'St_ds': fit_ds["ds_part/St"].value,
        'f_offset': fit_ds["ds_part/f_offset"].value,
        'd_offset': fit_ds["ds_part/d_offset"].value
    }

    return ds_fit_dict, area_ds


def fitting_ss(filename_i, input_settings, export_data, input_fitting, i_start, i_end, Force_Distance, fix, max_range, derivation_array, F_low):
    global model_ss
    global ss_fit_dict

    start_fitting_region = np.where(derivation_array[:, 1] == i_start)
    end_fitting_region = np.where(derivation_array[:, 1] == i_end)
    start_fitting_region = start_fitting_region[0][0]
    end_fitting_region = end_fitting_region[0][0]

    raw_f_fitting_region = Force_Distance[start_fitting_region * input_settings['step_d'] + len(F_low):end_fitting_region * input_settings['step_d'] + len(F_low), 0]
    raw_d_fitting_region = Force_Distance[start_fitting_region * input_settings['step_d'] + len(F_low):end_fitting_region * input_settings['step_d'] + len(F_low), 1]

    # downsample the data used for fitting to 200 datapoints
    if len(raw_f_fitting_region) > 200:
        f_fitting_region_ss = raw_f_fitting_region[::int(len(raw_f_fitting_region) / 200)]
        d_fitting_region_ss = raw_d_fitting_region[::int(len(raw_f_fitting_region) / 200)]
    else:
        f_fitting_region_ss = raw_f_fitting_region
        d_fitting_region_ss = raw_d_fitting_region

    if input_fitting['WLC+FJC'] == 1:
        model_ss = lk.odijk("DNA_2") + lk.freely_jointed_chain("RNA")
    elif input_fitting['WLC+WLC'] == 1:
        model_ss = lk.odijk("DNA_2") + lk.odijk("RNA")

    model_ss = model_ss.invert().subtract_independent_offset() + lk.force_offset("DNA")
    fit_ss = lk.FdFit(model_ss)

    fit_ss.add_data("ss_part", f_fitting_region_ss, d_fitting_region_ss)

    # ds part parameters
    # Persistance length bounds

    # Lp_ds_range=fit_ds["DNA/Lp"].value/10
    fit_ss["DNA_2/Lp"].value = ds_fit_dict['Lp_ds']
    fit_ss["DNA_2/Lp"].upper_bound = ds_fit_dict['Lp_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lp"].lower_bound = ds_fit_dict['Lp_ds'] * (1 - max_range / 100)
    # if fix==1:
    fit_ss["DNA_2/Lp"].fixed = 'True'

    fit_ss["DNA/f_offset"].upper_bound = 5
    fit_ss["DNA/f_offset"].lower_bound = -5
    fit_ss["DNA/f_offset"].value = ds_fit_dict['f_offset']
    fit_ss["DNA/f_offset"].fixed = 'True'

    fit_ss["inv(DNA_2_with_RNA)/d_offset"].value = ds_fit_dict['d_offset']
    fit_ss["inv(DNA_2_with_RNA)/d_offset"].fixed = 'True'

    # contour length
    # Lc_ds_range=Lc_initial_guess/100 # nm
    fit_ss["DNA_2/Lc"].upper_bound = ds_fit_dict['Lc_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/Lc"].lower_bound = ds_fit_dict['Lc_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/Lc"].value = ds_fit_dict['Lc_ds']
    fit_ss["DNA_2/Lc"].unit = 'nm'
    # if fix==1:
    fit_ss["DNA_2/Lc"].fixed = 'True'

    # stifness

    fit_ss["DNA_2/St"].upper_bound = ds_fit_dict['St_ds'] * (1 + max_range / 100)
    fit_ss["DNA_2/St"].lower_bound = ds_fit_dict['St_ds'] * (1 - max_range / 100)
    fit_ss["DNA_2/St"].value = ds_fit_dict['St_ds']
    if fix == 1:
        fit_ss["DNA_2/St"].fixed = 'True'

    # ss part parameters

    # Persistance length bounds
    fit_ss["RNA/Lp"].value = input_fitting['lp_ss']
    fit_ss["RNA/Lp"].lower_bound = 0.8
    fit_ss["RNA/Lp"].upper_bound = 2
    if fix == 1:
        fit_ss["RNA/Lp"].fixed = 'True'

    # stiffnes
    fit_ss["RNA/St"].value = input_fitting['ss_stiff']
    fit_ss["RNA/St"].lower_bound = 300
    fit_ss["RNA/St"].upper_bound = 1500

    # contour length
    fit_ss["RNA/Lc"].value = input_fitting['lc_ss']
    fit_ss["RNA/Lc"].lower_bound = 0
    fit_ss["RNA/Lc"].upper_bound = input_fitting['lc_ss'] + 100

    fit_ss["RNA/Lc"].unit = 'nm'

    fit_ss.fit()
    print(fit_ss.params)

    # calculate the integrals of the fitted functions
    distance_integral_fit_start = np.arange(min(Force_Distance[:, 1]), i_start)
    ss_integral_start = model_ss(distance_integral_fit_start, fit_ss)
    area_ss_fit_start = simps(ss_integral_start)
    print("area_ss_start = " + str(area_ss_fit_start))

    distance_integral_fit_end = np.arange(min(Force_Distance[:, 1]), i_end)
    ss_integral_end = model_ss(distance_integral_fit_end, fit_ss)
    area_ss_fit_end = simps(ss_integral_end)
    print("area_ss_end = " + str(area_ss_fit_end))

    fit_qual = fit_ss.log_likelihood()

    if input_fitting["WLC+WLC"] == 1:
        fitting_model = "WLC+WLC"
    elif input_fitting["WLC+FJC"] == 1:
        fitting_model = "WLC+FJC"

    ss_fit_dict = {
        'filename': filename_i,
        'model': fitting_model,
        'log_likelihood': fit_qual,
        'Lc_ds': fit_ss["DNA_2/Lc"].value,
        'Lp_ds': fit_ss["DNA_2/Lp"].value,
        'St_ds': fit_ss["DNA_2/St"].value,
        'Lc_ss': fit_ss["RNA/Lc"].value,
        'Lc_ss_stderr': fit_ss["RNA/Lc"].stderr,
        'Lp_ss': fit_ss["RNA/Lp"].value,
        'St_ss': fit_ss["RNA/St"].value,
        'f_offset': fit_ss["DNA/f_offset"].value,
        'd_offset': fit_ss["inv(DNA_2_with_RNA)/d_offset"].value
    }

    return fit_ss, f_fitting_region_ss, d_fitting_region_ss, ss_fit_dict, area_ss_fit_start, area_ss_fit_end


def plot_fit(fit, start_force_ss, start_distance_ss, Force_Distance, save_folder, filename_i, start_time):
    distance = np.arange(min(Force_Distance[:, 1]), max(Force_Distance[:, 1]) + 50, 2)
    F_ds_model = model_ds(distance, fit_ds)

    legend_elements = [
        Line2D([0], [0], color='k', lw=1, alpha=0.85),
        Line2D([0], [0], color='r', lw=1),
        Line2D([0], [0], color='gray', linestyle='dashed', lw=1)
    ]

    plt.plot(Force_Distance[:, 1], Force_Distance[:, 0], 'k', alpha=0.85)
    plt.scatter(d_fitting_region_ds, f_fitting_region_ds, color='r', s=4)
    plt.plot(distance, F_ds_model, linestyle='dashed', color='gray')
    plt.ylabel("Force [pN]")
    plt.xlabel("Distance [nm]")
    plt.legend(legend_elements, ['FD-Curve', 'Part used for fitting', 'Fitted WLC model'])

    for i in range(0, len(fit)):
        F_ss_model = model_ss(distance, fit[i])
        plt.scatter(start_distance_ss[i], start_force_ss[i], s=4)
        plt.plot(distance, F_ss_model, linestyle='dashed', color='gray')

    plotname = save_folder + "/" + filename_i + "_fit_" + start_time + ".png"

    plt.savefig(plotname, dpi=600)
    plt.clf()
