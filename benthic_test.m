% Code for Bradley, Hülse, LaRowe, Arndt (2022, Nature Comms) 
% OMEN-SED 1.0 BENTHIC-MODEL Stand-alone matlab code
% Hülse et al (2017) GMD paper, incl. the RCM approximation of Pika et al. (2021, GMD)
% 

% benthic_test.m
% functions to run OMEN-SED and plot the results

% Command to run the model as used in this paper:
% SA_value = 1.0;
% [res, BE_depth, Flux_depth, Total_burial_depth, Total_burial_age, BE_age, Flux_age, Flux_TOC_swi, dxdy, depth_levels, age_levels] = benthic_test.test_benthic_BE_global(SA_value);


classdef benthic_test
    % test cases for benthic layer model
    
    properties
    end
    
    methods(Static)
        
        function swi = default_swi()
            % set default SWI conditions
            
            format longEng
            
            bsd = benthic_main();
            %bottom water concentrations
            swi.T = 8.0;                                        % temperature (degree C)
            
            % for 2G-model
            swi.C01_nonbio= 1.0*1e-2/12*bsd.rho_sed;            % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            swi.C02_nonbio= 1.0*1e-2/12*bsd.rho_sed;            % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            swi.Fnonbio1 = swi.C01_nonbio*(1-bsd.por)*bsd.w;    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
            swi.Fnonbio2 = swi.C02_nonbio*(1-bsd.por)*bsd.w;    % calculate flux [mol/(cm2 yr)] according non-bioturbated flux
            swi.C01 = swi.C01_nonbio;                           % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
            swi.C02 = swi.C02_nonbio;                           % resulting bioturbated SWI-concentration, to be calculated in benthic_zTOC.m
            
            % for nG-model
            swi.nG = 100;
            swi.p_a = 20.0;   % as in Dale ea 2015:  3e-4;
            swi.p_nu = 0.125;
            swi.C0 = 1.0 * 1e-2/12*bsd.rho_sed;                 % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)
            swi.TwoG_OM_model = false;
            
            %             swi.O20=300.0E-009;                                 % O2  concentration at SWI (mol/cm^3)
            %             swi.NO30=40.0e-9;                                   % NO3 concentration at SWI (mol/cm^3)
            %             swi.Nitrogen=true;                                  % calculate N (true/false)
            %             swi.NH40=10.0e-9;                                 	% NH4 concentration at SWI (mol/cm^3)
            %             swi.SO40=2.8E-005;                                	% SO4 concentration at SWI (mol/cm^3)
            %             swi.H2S0=0.0;                                       % H2S concentration at SWI (mol/cm^3)
            %             swi.PO40=40.0e-9;                                   % PO4 concentration at SWI (mol/cm^3)
            %             swi.Mflux0=365*0.2e-10*1/(1-bsd.por)*1/bsd.w;       % actually CONCENTRATION of M at the sediment [mol/cm3] : from flux input  365*0.2e-10 (mol/(cm2*yr))
            %             swi.DIC0=2.4E-006;                                 	% DIC concentration at SWI (mol/cm^3)
            %             swi.ALK0=2.4E-006;                                 	% ALK concentration at SWI (mol/cm^3)
            %             swi.S0=35;                                         	% Salinity at SWI (not used at the moment)
            %             swi.plot_PO4_DIC_ALK=true;
        end
        
        
        function test_a_w()
            wdepth_shallow = [25 50 75];
            wdepth = [wdepth_shallow (100:100:6000)];
            for i = 1:length(wdepth)
                w(i)=benthic_main.sedrate(wdepth(i));
                Dbio(i)=benthic_main.biorate(wdepth(i));
                a(i) = benthic_main.apparent_age(w(i));
                a_Boudreau(i) = benthic_main.apparent_age_Boudreau(w(i)*1000);
            end
            
            
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            fig1 = figure;
            plot(wdepth,w,'k-');
            xlabel('Depth (m)');
            ylabel('w (cm/yr)')
            txt1 = 'Using Depth-w relation from Burwicz et al. (2011)';
            text(500,0.1,txt1)
            print(fig1,'-depsc2', 'output/01_w_vs_Depth.eps');
            
            fig2 = figure;
            plot(w,a,'k-',w,a_Boudreau,'b--');
            xlabel('w (cm/yr)')
            ylabel('Initial age (yr)');
            txt2 = 'a-w relation from Arndt et al. (2013)';
            txt3 = 'a-w relation from Boudreau & Ruddick (1991)';
            hleg=legend(txt2, txt3,'Location','best');
            print(fig2,'-depsc2', 'output/02_Age_vs_w_log10.eps');
            
            fig3 = figure;
            plot(wdepth,a,'k-',wdepth, a_Boudreau,'b--');
            xlabel('Depth (m)')
            ylabel('Initial age (yr)');
            txt1 = 'Using Depth-w relation from Burwicz et al. (2011)';
            txt2 = 'a-w relation from Arndt et al. (2013)';
            txt3 = 'a-w relation from Boudreau & Ruddick (1991)';
            hleg=legend(txt2, txt3,'Location','best');
            set(hleg,'FontSize',16)
            text(500,700,txt1)
            print(fig3,'-depsc2', 'output/03_Age_vs_Depth_log10.eps');
            
            
            fig6 = figure;
            plot(wdepth,Dbio,'k-');
            xlabel('Depth (m)');
            ylabel('Dbio (cm^2/yr)')
            txt1 = 'Using Depth-Dbio relation of Middelburg et al. (1997)';
            text(500,25,txt1)
            print(fig6,'-depsc2', 'output/04_Dbio_vs_Depth.eps');
            
        end
        
         function run_OMEN_RCM_twice()
            % run OMEN-SED with default SWI conditions as in default_swi()
            % plot two TOC profiles: 2nd with 2 x k for z < zbio
            clear
            swi = benthic_test.default_swi();
            swi.TwoG_OM_model = false;
            swi.flux = false;
            swi.IntConst_GMD = true;            %            swi.nG = 20;
            % parameter a given in years e.g. 10^2, not as log10
            %           swi.p_a = 50;
            %           swi.p_nu = 0.35;
            res=benthic_test.test_benthic(1,swi);
            % set date-time or string going into plot function
            str_date = [num2str(res.swi.nG) 'G_a=' num2str(res.swi.p_a) '_nu=' num2str(res.swi.p_nu)];
            %            benthic_test.plot_column(res, false, swi, str_date)
            benthic_test.plot_TOC_twice(res, false, swi, str_date);
            
            
         end       
        
        function run_OMEN_RCM()
            % run OMEN-SED with default SWI conditions as in default_swi()
            clear
            swi = benthic_test.default_swi();
            swi.TwoG_OM_model = false;
            swi.flux = false;
            swi.IntConst_GMD = true;            %            swi.nG = 20;
            % parameter a given in years e.g. 10^2, not as log10
            %           swi.p_a = 50;
            %           swi.p_nu = 0.35;
            res=benthic_test.test_benthic(1,swi);
            % set date-time or string going into plot function
            str_date = [num2str(res.swi.nG) 'G_a=' num2str(res.swi.p_a) '_nu=' num2str(res.swi.p_nu)];
            %            benthic_test.plot_column(res, false, swi, str_date)
            benthic_test.plot_TOC(res, false, swi, str_date);
            
            % calculate depth integrated OM degradation rates
            swi.C0i = res.swi.C0i;
            Cox_rate.Cox_total = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1, res.bsd, swi, res);
            
            % calculate mean OM concentration in upper x cm
            [C_10, C1_11] = res.zTOC_RCM.calcC( 10, res.bsd, res.swi, res);
            OM_10=C_10* 100*12/res.bsd.rho_sed;
            x = 10;
            Mean_OM = 1/x * 100*12/res.bsd.rho_sed*res.zTOC_RCM.calcOM(0.0, x, 1, res.bsd, swi, res);
            
            % calcuate flux at SWI:
            [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
            [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(100, res.bsd, res.swi, res);
            fprintf('Flux at SWI %g \n',  F_TOC_swi);
            fprintf('Flux at 1m %g \n',  F_TOC_inf);
            fprintf('BE(1m) %g \n',  F_TOC_inf/F_TOC_swi*100);
            
        end
        
        function run_OMEN_RCM_Obs(pObs)
            % run OMEN-SED with default SWI conditions as in default_swi()
            
            swi = benthic_test.default_swi();
            swi.TwoG_OM_model = false;
            swi.flux = false;
            swi.IntConst_GMD = true;            %            swi.nG = 20;
            % parameter a given in years e.g. 10^2, not as log10
            %           swi.p_a = 50;
            %           swi.p_nu = 0.35;
            
            %sediment characteristics - run test_benthic() with specific BCs
            switch pObs
                case 1  % Self: OMEXDIA_2809_108m
                    wdepth = 108;
                    res.bsd = benthic_main(1,wdepth);
                    res.bsd.usescalarcode = true;                    
                    res.bsd.zbio=1.0;   % because it is an anoxic site
                    res.bsd.por=0.599;  % shelf default
                    swi.p_a = 0.1;
                    swi.C0 = 4.44 * 1e-2/12*res.bsd.rho_sed;                 % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)              
                    
                    %% set date-time or string going into plot function
                    str_date = [num2str(swi.nG) 'G_a=' num2str(swi.p_a) '_Shelf_108m'];

                case 2  % Margin: Sanata Barbara 585m
                    wdepth = 585;
                    res.bsd = benthic_main(1,wdepth);
                    res.bsd.usescalarcode = true;                    
                    res.bsd.zbio=1.0;   % because it is an anoxic site
                    res.bsd.por=0.695;  % shelf default
                    swi.p_a = 10.0;
                    swi.C0 = 5.5 * 1e-2/12*res.bsd.rho_sed;                 % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)                         
                   
                    %% set date-time or string going into plot function
                    str_date = [num2str(swi.nG) 'G_a=' num2str(swi.p_a) '_Margin1_585m'];
                    
                case 3  % Margin: OMEXDIA_2809_2213m
                    wdepth = 2213;
                    res.bsd = benthic_main(1,wdepth);
                    res.bsd.usescalarcode = true;                    
                    res.bsd.zbio=10.0;   % because it is an anoxic site
                    res.bsd.por=0.695;  % shelf default
                    swi.p_a = 1.0;
                    swi.C0 = 0.9452 * 1e-2/12*res.bsd.rho_sed;                 % mean of upper 1cm TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)                         
                   
                    %% set date-time or string going into plot function
                    str_date = [num2str(swi.nG) 'G_a=' num2str(swi.p_a) '_Margin2_2213m'];
                    
                case 4  % Abyss: Estes et al. (2019, NatGeosc) NA 11 & NA 12
                    wdepth = (5557+5367)/2;
                    res.bsd = benthic_main(1,wdepth);
%                    res.bsd.w = 0.12/1000;
                    res.bsd.usescalarcode = true;                    
                    res.bsd.zbio=10.0;   % because it is an anoxic site
                    res.bsd.zinf=3000.0;   % because it is an anoxic site
                    res.bsd.por=0.85;  % shelf default
                    swi.p_a = 20.0;
                    swi.C0 = 0.26 * 1e-2/12*res.bsd.rho_sed;                 % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)                         
                   
                    %% set date-time or string going into plot function
                    str_date = [num2str(swi.nG) 'G_a=' num2str(swi.p_a) '_Abyss1_NA'];
                    
               case 5  % Abyss: Estes et al. (2019, NatGeosc) SPG 1 & SPG 9
                    wdepth = (5699+4924)/2;
                    res.bsd = benthic_main(1,wdepth);
%                    res.bsd.w = 0.04/1000;
                    res.bsd.usescalarcode = true;                    
                    res.bsd.zbio=10.0;   % because it is an anoxic site
                    res.bsd.zinf=1000.0;   % because it is an anoxic site
                    res.bsd.por=0.85;  % shelf default
                    swi.p_a = 20.0;
                    swi.C0 = 0.29 * 1e-2/12*res.bsd.rho_sed;                 % TOC concentration at SWI (wt%) -> (mol/cm^3 bulk phase)                         
                   
                    %% set date-time or string going into plot function
                    str_date = [num2str(swi.nG) 'G_a=' num2str(swi.p_a) '_Abyss2_SP'];
                    
            end
            
            %            res=benthic_test.test_benthic(1,swi);
            %% test_benthic()                    
            res.swi = swi;
            
            % calculate
            if(swi.TwoG_OM_model)
                res.zTOC = benthic_zTOC(res.bsd);
            else
                res.zTOC_RCM = benthic_zTOC_RCM(res.bsd);
            end
            
            if(swi.TwoG_OM_model)
                res = res.zTOC.calc(res.bsd,res.swi, res);
                O2_demand_flux = -(res.swi.Fnonbio1+res.swi.Fnonbio2)*res.bsd.OC/((1-res.bsd.por)./res.bsd.por)
            else
                % Adding into on RCM for MultiG approach
                [res.zTOC_RCM.k, res.swi.C0i, res.swi.Fnonbioi] = benthic_test.RCM(res.bsd, res.swi);
                res = res.zTOC_RCM.calc(res.bsd,res.swi, res);
            end
                    
            benthic_test.plot_TOC_Obs(res, false, swi, str_date, pObs);
            
            % calculate depth integrated OM degradation rates
            swi.C0i = res.swi.C0i;
            Cox_rate.Cox_total = res.zTOC_RCM.calcReac(0.0, res.bsd.zinf, 1, res.bsd, swi, res);
            
            % calculate mean OM concentration in upper x cm
            [C_10, C1_11] = res.zTOC_RCM.calcC( 10, res.bsd, res.swi, res);
            OM_10=C_10* 100*12/res.bsd.rho_sed;
            x = 10;
            Mean_OM = 1/x * 100*12/res.bsd.rho_sed*res.zTOC_RCM.calcOM(0.0, x, 1, res.bsd, swi, res);
            
            % calcuate flux at SWI:
            [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
            [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(100, res.bsd, res.swi, res);
            fprintf('Flux at SWI %g \n',  F_TOC_swi);
            fprintf('Flux at 1m %g \n',  F_TOC_inf);
            fprintf('BE(1m) %g \n',  F_TOC_inf/F_TOC_swi*100);
            
        end
        
        function [res, BE_depth, Flux_depth, Total_burial_depth, Total_burial_age, BE_age, Flux_age, Flux_TOC_swi, dxdy, depth_levels, age_levels, str_outdir] = test_benthic_BE_global(SA_factor_in)
            %% Calculate the BE using TOC observations from Seiter and BC as in Bradley ea. 2020
            SA_factor = SA_factor_in;   % for sensitivity ana multiply the default porosity, Dbio, sed-rate with this factor [1.0 : for default]
            SA_factor_str=sprintf('%.1f',SA_factor);
            
            use_por_map = true;     % use the spatial map for porosity?
            save_results = true;
            use_Lee = true;
            
            
            warning('query','all')
            warning('off','all')
            warning
            %__________________________________________________________________________
            %   load data
            %__________________________________________________________________________
            if(use_Lee)
                addpath('./data/Lee_et_al_2019/')
                load('Lee_toc_lr_weighted.mat')  % Lee data mean weighted by grid-size, meanhas NaN for terrestial cells
                toc = Lee_toc_lr_weighted;
                load('lat_lr.mat')
                lat = lat_lr;
                load('long_lr.mat')
                long = long_lr;
                load sed_holo.mat
                sed_holo = sed_holo(1:end-1, 1:end-1);  % delete extra row and column
                load ABYSS_MAP.mat
                ABYSS_MAP = ABYSS_MAP(1:end-1, 1:end-1);  % delete extra row and column
                load SHELF_MAP.mat
                SHELF_MAP = SHELF_MAP(1:end-1, 1:end-1);  % delete extra row and column
                load MARGIN_MAP.mat
                MARGIN_MAP = MARGIN_MAP(1:end-1, 1:end-1);  % delete extra row and column
                load zholo.mat
                zholo = zholo(1:end-1, 1:end-1);  % delete extra row and column
                
                %                load('water_depth_updated_Lee.mat')
                %                water_depth_updated = -water_depth_updated_Lee;
                load('./data/RECCAP2/bathymetry_matrix_new_ud.mat')
                water_depth_updated = -bathymetry_matrix_new_ud;
                
            	load('./data/O2_BW_WOA2018_hr_Qdegr.mat');   % BW  O2 [muM = 10^-6 mol/kg]
                loc_BW_O2 = O2_BW_WOA2018_hr_Qdegr;
                
                
                if(use_por_map)
                    load('./data/RECCAP2/porosity_matrix_new_ud.mat')
                    porosity = porosity_matrix_new_ud;
                end
            else
                addpath('./data/')
                load 'toc_NaN.mat'  % has NaN for terrestial cells
                load 'lat.dat'
                load 'long.dat'
                load sed_holo.mat
                load ABYSS_MAP.mat
                load SHELF_MAP.mat
                load MARGIN_MAP.mat
                load zholo.mat
                load('water_depth_updated.mat')
                water_depth_updated = -water_depth_updated;
                
                if(use_por_map)
                    load('./data/RECCAP2/porosity_matrix_ex.mat')
                    porosity = porosity_matrix_ex;
                end
            end
            
            % specify depth-levels to calculate BE and so on
            depth_levels = [11 15 20 30 40 50 100 200 300 400 500 1000];  % in cm
            %  	depth_levels = [100];  % in cm
            % specify size of results
            BE_depth = cell(length(depth_levels), 1);
            Flux_depth = cell(length(depth_levels), 1);
            % Flux_depth_perArea = cell(length(depth_levels), 1);
            
            % specify age-levels to calculate BE and so on
            age_levels = [0.01 0.1 0.5 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 1000];  % in kyr
            if(false)   % if the combined matrix of age layers does not exist yet
                Age_layers = cell(length(age_levels), 1);
                i=1;
                load('data/z_10yr.mat');
                Age_layers{i}=z_10yr;
                load('data/z_100yr.mat');
                Age_layers{i+1}=z_100yr;
                load('data/z_500yr.mat');
                Age_layers{i+2}=z_500yr;
                load('data/z_1kyr.mat');
                Age_layers{i+3}=z_1kyr;
                load('data/z_2kyr.mat');
                Age_layers{i+4}=z_2kyr;
                load('data/z_3kyr.mat');
                Age_layers{i+5}=z_3kyr;
                load('data/z_4kyr.mat');
                Age_layers{i+6}=z_4kyr;
                load('data/z_5kyr.mat');
                Age_layers{i+7}=z_5kyr;
                load('data/z_6kyr.mat');
                Age_layers{i+8}=z_6kyr;
                load('data/z_7kyr.mat');
                Age_layers{i+9}=z_7kyr;
                load('data/z_8kyr.mat');
                Age_layers{i+10}=z_8kyr;
                load('data/z_9kyr.mat');
                Age_layers{i+11}=z_9kyr;
                load('data/z_10kyr.mat');
                Age_layers{i+12}=z_10kyr;
                load('data/z_20kyr.mat');
                Age_layers{i+13}=z_20kyr;
                load('data/z_30kyr.mat');
                Age_layers{i+14}=z_30kyr;
                load('data/z_40kyr.mat');
                Age_layers{i+15}=z_40kyr;
                load('data/z_50kyr.mat');
                Age_layers{i+16}=z_50kyr;
                load('data/z_100kyr.mat');
                Age_layers{i+17}=z_100kyr;
                load('data/z_1Myr.mat');
                Age_layers{i+18}=z_1Myr;
                
                save('data/Age_layers_combined.mat' , 'Age_layers')
            else
                load('data/Age_layers_combined.mat');
            end
            
            % specify size of results
            BE_age = cell(length(age_levels), 1);
            Flux_age = cell(length(age_levels), 1);
            % Flux_depth_perArea = cell(length(depth_levels), 1);
            
            
            zbio = 10;
            [m,n]=size(toc);
            
            zmax_sed_holo = max(max(zholo))+zbio;     % define maximum depth of sediment layer (in cm)
            zmax_sed_holo_round = ceil(zmax_sed_holo/100)*100;  % round to nearest meter
            
            ncl = 1;
            res.bsd = benthic_main(ncl);
            res.bsd.zinf = zmax_sed_holo_round;
            
            
            res.bsd.usescalarcode = ncl==1;
            
            swi = benthic_test.default_swi();
            swi.TwoG_OM_model = false;
            swi.flux = false;
            swi.IntConst_GMD = true;
            %__________________________________________________________________________
            %   set specific parameters
            %__________________________________________________________________________
            holocene=11700;                                     %age holocene (yrs)
            pleistocene=2.58e6;                                 %age pleistocene (yrs)
            %        zbio=10;
            
            res.swi = swi;
            
            xstart=1; %55;
            xstop=m;  %55;        %lat
            ystart=1;  % 6;
            ystop=n;  % 6;        %long (or the other way around- who knows!)
            
            
            for x = 1:m
                %         	for x = xstart:xstop
                x
                % calculate volume
                % convert deg to cm CODAS package (by E.Firing,et al.)
                rlat = lat(x) * pi/180;
                m = 111132.09  - 566.05 * cos(2 * rlat)+ 1.2 * cos(4 * rlat);
                dy = 0.25*m*100.0; %cm
                p = 111415.13 * cos(rlat) - 94.55 * cos(3 * rlat);
                dx = 0.25*p*100.0; %cm
                
                for y = 1:n
                    %                for y = ystart:ystop
                    %                   y
                    % set terrestrial lat-long to NaN
                    %                    if isnan(sed_holo(x,y))
                    if(~use_por_map)
                        porosity(x,y)=999;  % set to arbitrary non-nan value if no spatial map for porosity is used for next if-check
                    end
                    
                    if ((isnan(toc(x,y))))    % now check for sed_holo or toc = NaN
                        
                        toc(x,y)=NaN;
                        dxdy(x,y)   = NaN;
                        
                        
                        Flux_TOC_swi(x,y) = NaN;
                        Flux_TOC_swi_perArea(x,y) = NaN;
                        %Flux_depth_perArea{i}(x,y) = NaN;
                        
                        for i=1:length(depth_levels)
                            BE_depth{i}(x,y) =NaN;
                            Flux_depth{i}(x,y)=NaN;
                        end
                        
                        for j=1:length(age_levels)
                            BE_age{j}(x,y) =NaN;
                            Flux_age{j}(x,y)=NaN;
                        end
                    else
                        
                        res.swi.p_a=0.1*SHELF_MAP(x,y)+1.0*MARGIN_MAP(x,y)+20.0*ABYSS_MAP(x,y);    %reactive continuum model a=10?-1-10?5
                      	res.bsd.zbio = 10.0;
                        % check for low oxygen, i.e. O2 < 60 muM 
                        if(loc_BW_O2(x,y)<60)
                            res.swi.p_a=1.0*SHELF_MAP(x,y)+10.0*MARGIN_MAP(x,y)+20.0*ABYSS_MAP(x,y);    
                            res.bsd.zbio = 1.0;
                        end
                        %         res.swi.p_nu =0.125*SHELF_MAP(i,j)+0.125*MARGIN_MAP(i,j)+0.125*ABYSS_MAP(i,j); %reactive continuum model nu=0.01-1, mostly around 0.125

                        
                        if(~isnan(water_depth_updated(x,y)))
                            res.bsd.Dbio=benthic_main.biorate(water_depth_updated(x,y));
                        else
                            
                            res.bsd.Dbio=(27.6*SHELF_MAP(x,y)+8.1*MARGIN_MAP(x,y)+0.59*ABYSS_MAP(x,y));  % Bioturbation coefficient
                        end
                        if(use_por_map)
                            if(~isnan(porosity(x,y)))
                                res.bsd.por=SA_factor*porosity(x,y)/100;
                            else
                                res.bsd.por=SA_factor*(0.599*SHELF_MAP(x,y)+0.695*MARGIN_MAP(x,y)+0.755*ABYSS_MAP(x,y));   % use mean value
                            end
                        else
                            %                            res.bsd.por=0.45*SHELF_MAP(x,y)+0.74*MARGIN_MAP(x,y)+0.7*ABYSS_MAP(x,y);   % porosity at SWI as in Sci Adv
                            %                            res.bsd.por=0.53*SHELF_MAP(x,y)+0.613*MARGIN_MAP(x,y)+0.705*ABYSS_MAP(x,y);    % porosity at SWI low values
                            res.bsd.por=0.66*SHELF_MAP(x,y)+0.73*MARGIN_MAP(x,y)+0.75*ABYSS_MAP(x,y);    % porosity at SWI high values
                        end
                        if(~isnan(sed_holo(x,y)))
                            res.bsd.w=(2 - SA_factor)*sed_holo(x,y);       % 2 - factor, so 0.9 and 1.1 is  swaped as higher sed-rate causes more burial
                        else    % use default
                            res.bsd.w=(2 - SA_factor)*benthic_main.sedrate(water_depth_updated(x,y));
                        end
                                                
                        % initialize parameters of analytical solution
                        conv=2.5/100*(1-res.bsd.por);   % convert wt% -> g/cm3 (total sediment)   -- this is what Sandra did for James
                        conv=2.5/100;   % convert wt% -> g/cm3 (total sediment)     -- I think in my model I need it as bulk sediment because I apply *(1-por) when I calculate the flux in calcCflx() !!??
                        res.swi.C0=toc(x,y)*conv;          % POC at SWI (wt% -> g/cm3(total sediment))
                                               
                        % initialize & calculate
                        res.zTOC_RCM = benthic_zTOC_RCM(res.bsd);
                        % Adding into on RCM for MultiG approach
                        [res.zTOC_RCM.k, res.swi.C0i, res.swi.Fnonbioi] = benthic_test.RCM(res.bsd, res.swi);
                        res = res.zTOC_RCM.calc(res.bsd,res.swi, res);
                        
                        dxdy(x,y)   = dx*dy;    % cm^2 per grid-cell
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% calculate TOC fluxes
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        %%%%%%%%%%%%%%%%%
                        % flux through SWI (just for BE calculation - total value is way too high -- maybe a few very wrong grid-cells!?)
                        [F_TOC_swi, F_TOC1_swi] = res.zTOC_RCM.calcCflx(0, res.bsd, res.swi, res);
                        Flux_TOC_swi(x,y) = nansum(F_TOC1_swi)*dx*dy;  %   in g/yr
                        Flux_TOC_swi_perArea(x,y) = nansum(F_TOC1_swi);  %
                        
                        
                        %%%%%%%%%%%%%%%%%
                        % through different depth levels - calculate maps of BE and Fluxes
                        for i=1:length(depth_levels)
                            %  depth_levels = [11 15 20 30 40 50 100 200 300 400 500 1000];  % in cm
                            [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(depth_levels(i), res.bsd, res.swi, res);
                            Flux_depth{i}(x,y) = nansum(F_TOC1_inf)*dx*dy;  %   in g/yr
                            % Flux_depth_perArea{i}(x,y) = nansum(F_TOC1_inf);  %   in g /(cm^2 yr)
                            BE_depth{i}(x,y) = nansum(F_TOC1_inf)/nansum(F_TOC1_swi)*100;
                        end
                        
                        %%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%
                        % through different age-levels - calculate maps of BE and Fluxes
                        for j=1:length(age_levels)
                            %  age_levels = [0.1 0.5 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 1000];  % in kyr
                            % NOTE: add zbio as the age related depth-layers start with 0 = zbio
                            [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(Age_layers{j}(x,y) +  zbio, res.bsd, res.swi, res);
                            Flux_age{j}(x,y) = nansum(F_TOC1_inf)*dx*dy;  %   in g/yr
                            % Flux_depth_perArea{i}(x,y) = nansum(F_TOC1_inf);  %   in g /(cm^2 yr)
                            BE_age{j}(x,y) = nansum(F_TOC1_inf)/nansum(F_TOC1_swi)*100;
                        end                        
                    end
                    
                    
                end
            end
            
            % calculate global total burial per depth level
            for i=1:length(depth_levels)
                Total_burial_depth(i) = nansum(nansum(Flux_depth{i}));
            end
            
            % calculate global total burial per age level
            for j=1:length(age_levels)
                Total_burial_age(j) = nansum(nansum(Flux_age{j}));
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % save results
            if(save_results)
                % find current path
                str_current_path = pwd;
                str_outdir = ['SpatialPor_Lee_wdepth_combinedSA_', SA_factor_str ,'_por'];
                par_pathout = ['output/Exp_results/' str_outdir];
                par_pathout = [str_current_path '/' par_pathout];
                if ~(exist(par_pathout,'dir') == 7), mkdir(par_pathout);  end
                
                save([par_pathout, '/BE_depth.mat'] , 'BE_depth')
                save([par_pathout, '/Flux_depth.mat'] , 'Flux_depth')
                
                save([par_pathout, '/BE_age.mat' ], 'BE_age')
                save([par_pathout, '/Flux_age.mat'] , 'Flux_age')
                
                save([par_pathout, '/Flux_TOC_swi.mat' ], 'Flux_TOC_swi')
                
                save([par_pathout, '/dxdy.mat' ], 'dxdy')
                
                save([par_pathout, '/depth_levels.mat'] , 'depth_levels')
                save([par_pathout, '/age_level.mat'] , 'age_levels')
                
                save([par_pathout, '/Total_burial_depth.mat'] , 'Total_burial_depth')
                save([par_pathout, '/Total_burial_age.mat' ], 'Total_burial_age')
                

            end
            
            % analyze and plot results
            fun_analyze_plot_results()
        end
        
                
        function res = test_benthic( ncl, swi )
            loc_BW_O2_anoxia = 5.0e-9;       	% set to 5.0 nanomol/cm^3
            if nargin < 1
                ncl = 1;
            end
            
            res.bsd = benthic_main(ncl);
            res.bsd.usescalarcode = ncl==1;
            
            
            if nargin < 2 || isempty(swi)
                swi = benthic_test.default_swi();
            end
            
            if ncl > 1  % set up O2 gradient for testing
                O20 = swi.O20;
                for i = 1:ncl
                    swi.O20(i) = 10*(i-1)/(ncl-1)*O20;
                end
            end
            
            res.swi = swi;
            
            % calculate
            if(swi.TwoG_OM_model)
                res.zTOC = benthic_zTOC(res.bsd);
            else
                res.zTOC_RCM = benthic_zTOC_RCM(res.bsd);
            end
            
            if(swi.TwoG_OM_model)
                res = res.zTOC.calc(res.bsd,res.swi, res);
                O2_demand_flux = -(res.swi.Fnonbio1+res.swi.Fnonbio2)*res.bsd.OC/((1-res.bsd.por)./res.bsd.por)
            else
                % Adding into on RCM for MultiG approach
                [res.zTOC_RCM.k, res.swi.C0i, res.swi.Fnonbioi] = benthic_test.RCM(res.bsd, res.swi);
                res = res.zTOC_RCM.calc(res.bsd,res.swi, res);
            end
            
        end
        
        function plot_TOC(res, ~, swi, str_date)
            fig1=figure;
            
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            bsd = res.bsd;
            zgrid = 0:1.0:bsd.zinf;
            % for smooth 1m sediment column: zgrid = 0:0.1:bsd.zinf;
            
            %%% TOC
            for i=1:length(zgrid)
                [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
            end
            % Plot TOC fractions with colorpalette linspecer
            subplot(121)
            color = linspecer(swi.nG);
            for G = 1:swi.nG
                plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                hold on
            end
            
            % Plot sum (TOC)
            subplot(122)
            plot(100*C*12/bsd.rho_sed, -zgrid, ':k')
            hold on
            t=xlim;         % to draw penetration depths the correct lengths
            plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
            xlim(xlim)
            hold off
            xlabel ('TOC (wt%)')
            ylabel('Depth (cm)')
            %            title('Total TOC (wt%)')
            
            % save Figure
            print(fig1,'-depsc2', ['0_TOC_PROFILES_' str_date '_Dbio10.eps']);
            
            
        end
        
        function plot_TOC_twice(res, ~, swi, str_date)
            fig1=figure;
            
            set(0,'defaultLineLineWidth', 2)
            set(0,'DefaultAxesFontSize',12)
            
            bsd = res.bsd;
            zgrid = 0:1.0:bsd.zinf;
            % for smooth 1m sediment column: zgrid = 0:0.1:bsd.zinf;
            
            %%% TOC
            for i=1:length(zgrid)
                [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
            end
            % change a1, b1, A1, B1, A2 for bioturbated layer, i.e. k = k x 2
            res.rTOC_RCM.a1=(bsd.w-sqrt(bsd.w.^2+4.*res.zTOC_RCM.DC1.*res.zTOC_RCM.k.*2))./(2.*res.zTOC_RCM.DC1);
            res.rTOC_RCM.b1=(bsd.w+sqrt(bsd.w.^2+4.*res.zTOC_RCM.DC1.*res.zTOC_RCM.k.*2))./(2.*res.zTOC_RCM.DC1);
            
            res.rTOC_RCM.A1 =-(res.swi.C0i.*res.rTOC_RCM.b1.*exp(res.rTOC_RCM.b1.*bsd.zbio))./(res.rTOC_RCM.a1.*exp(res.rTOC_RCM.a1.*bsd.zbio)-res.rTOC_RCM.b1.*exp(res.rTOC_RCM.b1.*bsd.zbio)+bsd.tol_const);
            res.rTOC_RCM.B1 = res.swi.C0i-res.rTOC_RCM.A1;
            res.rTOC_RCM.A2 =(res.rTOC_RCM.A1.*(exp(res.rTOC_RCM.a1.*bsd.zbio)-exp(res.rTOC_RCM.b1.*bsd.zbio))+res.swi.C0i.*exp(res.rTOC_RCM.b1.*bsd.zbio))./(exp(res.rTOC_RCM.a2.*bsd.zbio)+bsd.tol_const);
            
            
            %%% TOC
            for i=1:length(zgrid)
                [C_2(i), C1_2(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                [Cflx_2(i), C1flx_2(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
            end
            
            
            if false
                % Plot TOC fractions with colorpalette linspecer
                subplot(121)
                color = linspecer(swi.nG);
                for G = 1:swi.nG
                    plot(100*C1(:,G)*12/bsd.rho_sed, -zgrid,'Color',color(G,:))
                    hold on
                end
            end
            
            % Plot sum (TOC)
            subplot(122)
            plot(100*C*12/bsd.rho_sed, -zgrid, '--k')
            hold on
            plot(100*C_2*12/bsd.rho_sed, -zgrid, ':k')
            hold on
            t=xlim;         % to draw penetration depths the correct lengths
            plot([0,t(1,2)], [-bsd.zbio,-bsd.zbio], 'k--')
            xlim(xlim)
            hold off
            xlabel ('TOC (wt%)')
            ylabel('Depth (cm)')
            %            title('Total TOC (wt%)')
            
            % save Figure
            print(fig1,'-depsc2', ['0_TOC_PROFILES_' str_date '_Dbio10.eps']);
            
            
        end        
        
        function plot_TOC_Obs(res, ~, swi, str_date, pObs)
            fig1=figure;
            
            switch pObs
                case 1  % Shelf 108m OMEXDIA_2809_108m 
                  	data.TOC=xlsread('./data/POC_Profiles/5_PE138_99-06_108m.xlsx','Corg','C2:D24');     % in wt%
                case 2
                    data.TOC=load('./data/POC_Profiles/Reimers_585m_BC68_Corg.dat','ascii');
                case 3
                    data.TOC=xlsread('./data/POC_Profiles/2_PE121_98-4_2213m.xlsx','Corg','C2:D39');
                case 4
                    data.TOC=xlsread('./data/POC_Profiles/Estes_ea_2019_Gyres.xlsx','plot','J44:K62');
                    data.TOC2=xlsread('./data/POC_Profiles/Estes_ea_2019_Gyres.xlsx','plot','J63:K80');
                case 5
                    data.TOC=xlsread('./data/POC_Profiles/Estes_ea_2019_Gyres.xlsx','plot','J132:K141');
                    data.TOC2=xlsread('./data/POC_Profiles/Estes_ea_2019_Gyres.xlsx','plot','J179:K187');

            end
            
            set(0,'defaultLineLineWidth', 3)
            set(0,'DefaultAxesFontSize',14)
            
            bsd = res.bsd;
            zgrid = 0:1.0:bsd.zinf;
            % for smooth 1m sediment column: zgrid = 0:0.1:bsd.zinf;
            
            %%% TOC
            for i=1:length(zgrid)
                [C(i), C1(i,:)] = res.zTOC_RCM.calcC( zgrid(i), bsd, res.swi, res);
                [Cflx(i), C1flx(i,:)] = res.zTOC_RCM.calcCflx( zgrid(i), bsd, res.swi, res);
            end
            
            % Plot sum (TOC)
            subplot(122)
            plot(100*C*12/bsd.rho_sed, -zgrid, 'b')
            hold on
            % Observationa TOC wt %
          	scatter(data.TOC(:,2), -data.TOC(:,1),'k','filled')
            if(pObs>=4)
                scatter(data.TOC2(:,2), -data.TOC2(:,1),'^','filled','MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5])                
            end
            t=xlim;         % to draw penetration depths the correct lengths
            plot([0,ceil(t(1,2))], [-bsd.zbio,-bsd.zbio], 'g--','linewidth', 1)
            if(pObs<3)
                 xlim([0 ceil(t(2))])
                ticks = [0:1:ceil(t(2))];
            elseif(pObs==3)
                xlim([0 1.5])
                ticks = [0:0.5:1.5];               
            else
                 xlim([0 0.3])
                ticks = [0 0.1 0.2 0.3];
               
            end
            xticks(ticks)
            hold off
            xlabel ('TOC (wt%)')
            ylabel('Depth (cm)')
            %            title('Total TOC (wt%)')
            
            % save Figure
            print(fig1,'-depsc2', ['0_TOC_PROFILES_' str_date '.eps']);
            
            
        end

        
        
        function [k, C0i, Fnonbioi] = RCM(bsd, swi)
            % For comparison with Dominik's 2G results
            if swi.nG == 2
                C0i(1:2) = 0.1 * 1e-2/12*bsd.rho_sed;       % TOC@SWI (wt%) -> (mol/cm^3 bulk phase), 2.5 sed.density (g/cm3) 0.1
                Fnonbioi = swi.C0*(1-bsd.por)*bsd.w;        % [mol/(cm2 yr)] according non-bioturbated flux
                k = [0.01 0.0001];
                swi.p_a = NaN;
                swi.p_nu = NaN;
            else
                %                 emin = log10(...
                %                     swi.p_nu./(swi.p_a+(bsd.zinf./bsd.w))...
                %                     ) - 1;                    % lower k limit for k-bins of multi-G approximation, i.e. k=1e-15 yr-1
                emin = -15;      % as in Dale ea 2015: -10;
                % subtracted 1 to be on the conservative side.
                emax = -log10(swi.p_a)+2;       % upper k limit for k-bins of multi-G approximation
                %                emax = 2;       % as in Dale ea 2015
                if emin >= emax;error('emin >= emax, this cannot be!');end
                %                 if emax >= log10(200);emax=log10(200);end
                
                k(1)= 10^(emin);
                kk(1)=10^(emin);
                F(1) = gammainc(swi.p_a*10^emin,swi.p_nu,'lower');
                kk(swi.nG)=10^(emax);
                k(swi.nG)=10^(emax);
                F(swi.nG) = gammainc(swi.p_a*10^emax,swi.p_nu,'upper');
                
                % Define the b.c. for all the intermediate fractions
                
                G=2:swi.nG-1;
                ne=emin+(1:swi.nG-2).*(emax-emin)./(swi.nG-1);
                kk(2:swi.nG-1)=10.^ne;
                G_inc_0 = gammainc(swi.p_a*kk(1:swi.nG-2),swi.p_nu,'upper'); % G-1 = 1:end-2
                G_inc_1 = gammainc(swi.p_a*kk(2:swi.nG-1),swi.p_nu,'upper'); % G = 2:end-1
                F(2:swi.nG-1) = (G_inc_0 - G_inc_1);
                % calculate the mean degradation rate for the intermediate OM fractions
                k(G)=kk(1:swi.nG-2)+(kk(2:swi.nG-1)-kk(1:swi.nG-2))/2;
                F(F<=eps)=eps;
                if abs(sum(F)-1) > 0.0001;warning('F~=1!!');end
                Fnonbioi = F.* ( swi.C0*(1-bsd.por)*bsd.w ); % NonBioturbated SWI
                C0i = F.*swi.C0;
            end
        end
        
 
    end
    
end

