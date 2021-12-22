
The forward model predicts the accumulation of cosmogenic nuclide 10Be, 26Al, and 21Ne in the ice core and the above laying supraglacial debris during a single event of exposure, where sublimation, surface erosion and accumulation of supraglacial debris has occurred. The predicted nuclide concentrations are solved numerically using the default integration algorithm in MATLAB. The predicted nuclide concentrations are then fitted to the measured nuclide concentrations in the ice core. A best fit is found using the minimum constrained nonlinear multivariable optimizing function (fmincon) in MATLAB and optimize for the free parameters; (i) the age the ice was emplaced, (ii) sublimation rate of the ice since emplacement, (iii) surface erosion rate of the accumulating supraglacial debris, and (iv) the inherited nuclide concentrations in the englacial debris at the time of ice emplacement.

The forward model and statistical error obtained from a Monte Carlo simulation is rooted in the MC_optimizer.m script. Therefore, one must only run the ‘MC_optimizer.m’ script to obtain the modeled results. Note, before running the script, all variable in the 'Set/define all model parameters and dataset' section should be adjusted to fit the data set and desired model constraints.

For a successful run, the following MATLAB files are required

- get_burial_age.m
- get_misfit.m
- nonlcon.m
- P_of_z.m
- plot_forward_model.m
- plot_forward_model_sample_depth.m
- Precalculate_P_of_z.m
- objective_sublimation_model.m
- save_data_core1.m
- sublimation_model_params.m
- z_of_t.m
- plot_forward_model.m
- plot_forward_model_sample_depth.m

In addition, two m-files originally published as part of the online exposure age calculators described by Balco et al. (2008) are needed:

- antatm.m
- stone2000

And one file described in Balco (2017):

- P_mu_total_alpha1.m

