Copyright 2021 Helmholtz-Zentrum für Infektionsforschung GmbH

POTATO -- 2022-06-10

Developed by Lukáš Pekárek and Stefan Buck at the Helmholtz Institute for RNA-based Infection Research

in the research group REMI - Recoding Mechanisms in Infections.

Supervisor - Jun. Prof. Neva Caliskan


Description

This script processes Force-Distance Optical Tweezers data in an automated way, to find unfolding events,

perform fitting of polymer models and calculate the work applied onto the structure.

The script is developed to handle h5 raw data, produced by the C-Trap OT instrument from LUMICKS,

as well as any other FD data prepared in a CSV file (2 columns: Force(pN) - Distance(um or nm))

Furthermore the script can analyse single constant force markers.

The parameters can be changed in the GUI before each run.

Alternatively they can be changed permanently in the POTATO_config file.

For more detailed description please refer to: doi: https://doi.org/10.1101/2021.11.11.468103


Dependencies

Python3 Packages: csv, glob, h5py, json, lumicks.pylake, matplotlib, multiprocessing, numpy, os, pandas, pathlib, PIL, scipy, statistics, time, tkinter

There is also a standalone executable POTATO version for Windows available.


Navigation

The POTATO GUI is structured in five tabs: "Folder Analysis", "Show Single File", "Constant Force Analysis", "Advanced Settings"and "TOMATO - Manual analysis". For each analysis step, buttons are displayed in the different tabs and some of them can also be found in the drop-down menu.


Folder Analysis

In this tab, all force ramp experimental files in the specified folder will be analyzed automatically.

Input

POTATO supports three different input formats for folder analysis. The appropriate dataset has to be selected prior to the analysis. ATTENTION: When changing between input formats, all parameters are reset to the default values!

High frequency (Piezo Distance): Analyses the high frequency data (Piezo tracking) of all h5 files in a given directory. The data gathering frequency is derived directly out of the h5 files.
Low Frequency: Analyses the low frequency data of all h5 files in a given directory. The data gathering frequency is calculated directly out of the h5 files.
CSV (F/D): Analyses all csv files in a given directory. The architecture of these files need to consist out of two columns (Force and Distance) without headers. Force needs to be in pN. Distance can be either be in µm or nm. The script can process forward (unfolding) as well as reverse (folding) curves. If the script should distinguish, the reverse curves need to start with the highest distance and decrease from there. The data gathering frequency and all other parameters are derived from the user input in the GUI.
Parameters

Parameters can be changed directly in the Graphical User Interface (GUI). For more setting options refer to the "Advanced Settings" tab, which includes all adjustable parameters. When the parameters are optimized, default parameters can be changed in the POTATO_config file, so they will be loaded when the GUI is started. The parameters are read in once the analysis starts and for the force-ramp analysis the used parameters are exported in json format.

Output POTATO creates an "Analysis" folder with timestamp in the analysed directory. The "Refresh" button loads the last saved image and displays the progress in the GUI. In the "Advanced Settings" tab, several export settings can be set.

Processed FD data: Exports the down-sampled and filtered Force-Distance-Array in CSV format. The exported filename consists of the original filename and additional suffix "_smooth".
Plot: Exports a figure (PNG) containing - the visualized processed data with and without marked unfolding events - the corresponding force- and distance derivatives
Steps found: Exports the identified steps for each curve into a separate CSV file. The step length does NOT correspond to the contour length change, but is only the distance difference between step start and step end.
Fitting: Exports a plot with the fitted models and a table of the fitting parameters for each section in CSV format. When 'Fitting' is not selected, the script skips all fitting steps and therefore the analysis is much faster.
Total results: Exports all found steps for all analysed files to a single CSV file.


Show Single File

Input

Single FD-curves of all three input formats (h5-HF, h5-LF, CSV) can be displayed.

Output

A FD-curve of the original input values, as well as the down sampled values are plotted in the GUI. This may help identify potential causes of errors.


Constant Force Analysis

Input

First the constant force marker should be displayed in a time-distance plot with adjacent histogram (Button -> Display Constant Force Data). All three input formats (h5-HF, h5-LF, CSV) are supported. After providing an initial guess of the fitting parameters, the histogram of the same file can be fitted with gaussians (Button -> Fit Constant Force Data).

Parameters

The down sample and filtering parameters are applied also onto the constant force data. This can have a huge impact on the displayed plot. Furthermore, the axes have to be set manually and the fitting parameters have to be guessed (dependent on the # of expected gaussians). To guess the fitting parameters, they are entered from lowest to highest values, separated with a comma.

Output

POTATO creates an "Analysis_constantF" folder with timestamp in the analysed directory. Both, the time-distance plot and the fitted histogram are exported as separate figures in PNG format. In addition, the smoothened values as well as the histogram distribution are exported into CSV files. Last, a table with values of the fitted parameters is exported in CSV format.


Advanced Settings

This tab contains all the adjustable parameters in the POTATO. The parameters are divided into several groups based on the part of analysis, in which they are used.

Preprocessing

Downsample rate - Only every nth value is taken for analysis; speeds up subsequent processing.

Butterworth Filter degree - Defines the stringency of the filter.

Cut off frequency - Signals with a frequency above this threshold are suppressed.

Force threshold, pN - Values with lower force are excluded from the analysis.

Derivative

Step d - Characterizes the interval between two values used for numerical derivative calculation.

Data frequency, Hz - The frequency at which data is recorded.

Statistics

z-score force/distance - The number of standard deviations used to determine whether a given value is part of a normal distribution.

moving median window size - The number of values considered for each median calculation. Can be set independently for force and distance

SD difference threshold - Statistical analysis and data sorting are iterated until the difference between two consecutive SDs is below this value.

Fitting

"WLC+WLC" or "WLC+FJC" tick option - determines whether the unfolded regions of FD curves will be fitted with model combining two WLC models or WLC and FJC models, repectively.

dsLp, nm - Persistence length of the double-stranded (folded) part of the tethered construct.

dsLc, nm - Contour length of double-stranded (folded) part of the tethered construct.

dsK0, pN - Stretch modulus of double-stranded (folded) part of the tethered construct.

Force_offset, pN - Force offset of a given dataset; compensates for a shift in the dataset.

Distance_offset, nm - Distance offset of a given dataset; compensates for a shift in the dataset.

ssLp, nm - Persistence length of the single-stranded (unfolded) part of the tethered construct.

ssLc, nm - Contour length of single-stranded (unfolded) part of the tethered construct.

ssK0, pN - Stretch modulus of single-stranded (unfolded) part of the tethered construct.

Export

Consists of tick options for marking the data files to be exported (ticked) or not (unticked) during the analysis.

Processed FD data - Exports the down-sampled and filtered Force-Distance-Array in CSV format.The exported filename consists of the original filename and additional suffix "_smooth".

Plot - Exports a figure (PNG) containing - the visualized processed data with and without marked unfolding events - the corresponding force- and distance derivatives

Steps found - Exports the identified steps for each curve into a separate CSV file. The step length does NOT correspond to the contour length change, but is only the distance difference between step start and step end.

Total results - Exports all steps from both derivatives in CSV format.

Fitting - Exports a plot with the fitted models and a table of the fitting parameters for each section in CSV format. When 'Fitting' is not selected, the script skips all fitting steps and therefore the analysis is much faster.


TOMATO - manual analysis

This tab allows user to perform manual analysis of the FD curves. TOMATO uses the same input parameters as POTATO.

Choose the folder containing the curves to analyse with the "Choose folder" button. Attention the file format and preprocessing values are set in the "Analysis folder" tab.

Manually mark the start and end of each (un)folding step:

press "s" or click on "Set start" button --> click on the position of start of a step
press "e" or click on "Set end" button --> click on the position of end of the same step
press "ctrl+s" or click on "Save step" button to save marked step coordinates repeat until all steps are marked. Attention, the steps have to be marked in the correct order from the "lowest" to the "highest" step.
if a mistake occurs during the marking, the last saved value can be removed by the "Delete step" button.
Click on "Analyze curve" TOMATO proceeds with step coordinates input and performes the rest analysis as done by POTATO (fitting + work calculations). This may take some time during which the GUI is unresponsive. Optimize fitting parameters until a good fit is obtained. Entries of bad fits can be deleted in the "results table" by selecting them and pressing "Delete" The relevant parameters are loaded into the table below the figure. The column width can be adjusted to see the parameter labels.

Switch to the next curve of the folder with the arrow keys (left, right) or press "a" or "d". The "step table" will be reseted. Analysis is performed as explained above and the parameters are added to the "results table".

When all curves are manually analysed, the results table can be saved as csv file with the "Save results table" button.
