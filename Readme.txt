Readme for the OMEN-SED Burial Efficiency experiments for the paper
Bradley, Hülse, LaRowe, Arndt (2022, NatComms)

To excute and plot all experiments presented in the paper just run
> fun_execute_OMEN

Not that you need to install M_Map for plotting the maps (https://www.eoas.ubc.ca/~rich/map.html)

 
%%%%%%%%%%%%%%%%%%%  For SI figure

1) Plot 5 sediment profiles with observations
Execute
> benthic_test.run_OMEN_RCM_Obs(pObs)
with:
pObs = 1 for Self: OMEXDIA_2809_108m
pObs = 2 for Margin: Sanata Barbara 585m
pObs = 3 for Margin: OMEXDIA_2809_2213m
pObs = 4 for Abyss: Estes et al. (2019, NatGeosc) NA 11 & NA 12


2) To calculate and plot the EETs, change into directory ./safe_R1.1 and execute
> workflow_eet_OMEN_BurEff

