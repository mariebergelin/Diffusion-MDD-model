# He3-Diffusion-MDD-model

The following codes in this folder are provided for the sole purpose of evaluation of an accompanying scientific paper preprint entitled ’Diffusion kinetics of 3He in pyroxene and plagioclase and applications to cosmogenic exposure dating and paleothermometry in mafic rocks’ by Marie Bergelin and three other authors. It is not licensed for any other reuse or redistribution.

This folder contains the Matlab codes needed for (i) plotting diffusion experiments, (ii) optimization of the best-fit forward modeled MDD diffusion kinetics, and (iii) plotting of the MDD results, and includes the following files:

- optimizer_MDD_Domain.m
- unpack_Ddata.m
- initDiff_GUI.m
- objective_MDD_model.d
- nonlcon_MDD.m
- get_MDD_misfit.m
- plot_Ddata.m
- plot_MDD_forward_model.m
- make_MDD_fig_3_4.m
- make_MDD_fig_6_7.m

In addition, we include a 'Ddata' folder containing the diffusion experiment txt-files and two subfolders with the MDD results m-files described in the scientific paper. The 'MDDdata' folder contains the best-fit forward-modeled MDD results, whereas the 'MDDdom_results' folder contains m-files of the full workspace for the best-fit MDD results. 


- To run the forward MDD model optimization, run 'optimizer_MDD_Domain.m'

This script runs the forward MDD optimization and requires a txt-file containing the diffusion experiment results (see script for details). The forward model and statistical misfit is rooted in the 'optimizer_MDD_domain.m' script. Therefore, you only need to run this script to obtain the modeled MDD results.

When running: 
1. a file selection window will pop up, allowing you to select the txt-file containing the diffusion experiment results that you would like to fit the forward MDD model to.
2. the 'InitDiff_GUI' will appear. This gui allows you to input the time and temperature a sample has experienced during irradiation and/or storage. Up to three initial gas loss stages can be added. If values are left as zero, they will not be included in the forward MDD model. Similar, if you select 'cancel', the initial gas loss will be excluded from the Forward MDD model. 
3. At the end of the run, a window will appear and ask if you would like to save the data. If 'yes', two folders will be created: In one folder, the full workspace will be saved. In the other, the best-fit MDDdata results will be saved as a matlab structure, which later can be loaded in the 'explore_MDD.m' script for plotting of the forward MDD model results. If 'no', the data will not be saved. If files with identical samplename already exist, the window will ask if you would like to override the existing data. 

- To plot a diffusion experiment data, run 'plot_Ddata.m'

This function plots the diffusion experiment data (Ddata). When running, this script asks you to select a single diffusion experiment txt-file for plotting.

- To explore a given set of diffusion kinetics, run 'explore_MDD.m'

This script compares diffusion experiment results with a given set of diffusion kinetics, either defined within the code or selected as a m-file. In addition, this script allows you to (i) plot multiple diffusion experiments and (ii) explore the effect initial irradiation and/or storage loss conditions on the apparent diffusivity.

when running:
1. a file selection window will pop-up, allowing you to select the txt-file containing the diffusion experiment results that you would like to fit the forward MDD model to. If multiple txt-files are selected for plotting, the forward MDD model is calculated only for the first selected data.
2. If the diffusion kinetics parameters are defined within the script, then the 'InitDiff_GUI' will appear. This gui allows you to input the time and temperature of initil gas loss conditions to be include in the forward model. selceting 'cancel' or inputting 'zero' values will exclude initial gas loss from the forward MDD model.
3. If no diffusion kinetic parameters are defined (hence, commented out), then a window will pop-up to selected m-files containing MDD diffusion kinetic results. This m-file is the MDDdata optimization results saved when running the 'optimizer_MDD_Domain.m' script.











