%% Execute this to perform the experiments shown in 
%% Bradley, Hülse, LaRowe, Arndt (2022, NatComms)


%% plot boundary conditions
plot_BC_MS

%% Run global analysis and create summary output folders
% default values for porosity and sedimentation rate
SA_value = 1.0;
[res, BE_depth, Flux_depth, Total_burial_depth, Total_burial_age, BE_age, Flux_age, Flux_TOC_swi, dxdy, depth_levels, age_levels, str_outdir] = benthic_test.test_benthic_BE_global(SA_value);
% 0.9 x porosity and 1.1 x sedimentation rate
SA_value = 0.9;
[res, BE_depth, Flux_depth, Total_burial_depth, Total_burial_age, BE_age, Flux_age, Flux_TOC_swi, dxdy, depth_levels, age_levels, str_outdir] = benthic_test.test_benthic_BE_global(SA_value);
% 1.1 x porosity and 0.9 x sedimentation rate
SA_value = 1.1;
[res, BE_depth, Flux_depth, Total_burial_depth, Total_burial_age, BE_age, Flux_age, Flux_TOC_swi, dxdy, depth_levels, age_levels, str_outdir] = benthic_test.test_benthic_BE_global(SA_value);

%% plot figures

plot_total_burial_SA