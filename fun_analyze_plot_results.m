%% Plottng script for the results of the global simulations
% select via true/false what you want to plot
% You can either directly plot the return variables of the function
% benthic_test.test_benthic_BE_global -- If they are still in the workspace
% or select load_results to load output .mat files saved in the folder
% output/Exp_results

%% make output directory

str_current_path = pwd;
par_pathout = [str_current_path '/output/global/' str_outdir];
if ~(exist(par_pathout,'dir') == 7), mkdir(par_pathout);  end

use_Lee = true;

load_result = false;
new_plotting = true; % use Robinson projection with M_Map: A Mapping package for Matlab
path(path,'/home/domhu/Documents/MATLAB/M_Map');

plot_global_total_TOC_burial=true;         % global values vs age and depth level
plot_global_BE_depth=true;              % BE per sediment depth for 3 depth regimes


plot_BE_age = true;                    % plot BE for different depth layers
plot_BE_depth = true;                  % plot BE for different depth layers

% misc other plots
plot_Bflux_depth_total = false;    % plot total Bflux (in g/yr ) for different depth layers
plot_Bflux_depth_total_difference = false;   % plot difference burial flux at depth layers to 11cm
plot_Bflux_age_total = false;    % was true              % plot total Bflux (in g/yr ) for different age layers
plot_Difference_Bflux_age_total = false;     % plot difference burial flux at age-laters to 10years (below zbio)
plot_BE_depth_difference = false;        % plot BE difference to 11cm
plot_BE_age_difference = false;          % plot BE for difference to depth layer 10years
plot_BC = false;                         % plot zholo, zpleiso, SWI-TOC wt%
plot_Bflux_depth_perArea = false;       % plot Bflux (in g/(cm^2 yr) ) for different depth layers


% depth-levels (cm): [11,15,20,30,40,50,100,200,300,400,500,1000]
depth_levels = [11 15 20 30 40 50 100 200 300 400 500 1000];  % in cm
% age levels (kyrs): [0.01,0.1,0.5,1,2,3,4,5,6,7,8,9,10,20,30,40,50,100,1000]
age_levels = [0.01 0.1 0.5 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 1000];  % in kyr

%  load 'toc.dat'
    if(use_Lee)
                addpath('./data/Lee_et_al_2019/')
                load('lat_lr.mat')
                lat = lat_lr;
                load('long_lr.mat')
                long = long_lr;
                load ABYSS_MAP_Lee.mat
                ABYSS_MAP = ABYSS_MAP_Lee;
                load SHELF_MAP_Lee.mat
                SHELF_MAP = SHELF_MAP_Lee;
                load MARGIN_MAP_Lee.mat
                MARGIN_MAP = MARGIN_MAP_Lee;
                load zholo.mat            
                MARGIN_MAP = MARGIN_MAP;
        
    else
        load 'data/lat.dat'
        load 'data/long.dat'
        load ABYSS_MAP.mat
        load SHELF_MAP.mat
        load MARGIN_MAP.mat
    end

%
%  load('output/BE_depth.mat')
%
%  load('output/BE_age.mat')
%
%  load('output/Flux_TOC_swi.mat')
%
%  load('output/dxdy.mat')
%
%
% load('output/Total_burial_depth.mat' , 'Total_burial_depth')
% load('output/Total_burial_age.mat' , 'Total_burial_age')
%

str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];


if(plot_global_total_TOC_burial)
    
    % calculate Total burial per ocean area:
     if(load_result)
        load(['output/Exp_results/', str_outdir ,'/Flux_depth.mat'])
        load(['output/Exp_results/', str_outdir ,'/Flux_age.mat'])
        load(['output/Exp_results/', str_outdir ,'/Total_burial_depth.mat'] , 'Total_burial_depth')
        load(['output/Exp_results/', str_outdir ,'/Total_burial_age.mat'] , 'Total_burial_age')
        load(['output/Exp_results/', str_outdir ,'/depth_levels.mat'])
        load(['output/Exp_results/', str_outdir ,'/age_level.mat'])
    end
    
    % calculate global total burial per depth level
    % and total rate of change
    for i=1:length(depth_levels)  % [11,15,20,30,40,50,100,200,300,400,500,1000] cm
        Total_burial_depth_abyss(i) = nansum(nansum(Flux_depth{i}.*ABYSS_MAP));
        Total_burial_depth_shelf(i) = nansum(nansum(Flux_depth{i}.*SHELF_MAP));
        Total_burial_depth_margin(i) = nansum(nansum(Flux_depth{i}.*MARGIN_MAP));
    end
    
    % calculate global total burial per age level
    for j=1:length(age_levels)
        Total_burial_age_abyss(j) = nansum(nansum(Flux_age{j}.*ABYSS_MAP));
        Total_burial_age_shelf(j) = nansum(nansum(Flux_age{j}.*SHELF_MAP));
        Total_burial_age_margin(j) = nansum(nansum(Flux_age{j}.*MARGIN_MAP));
    end
    
    set(0,'defaultLineLineWidth', 2)
    set(0,'DefaultAxesFontSize',4)
    
    fig1=figure;
    set(gca,'FontSize',15)
    hold on
    box on
    plot(depth_levels, -Total_burial_depth./10^(12),'xk-')
    plot(depth_levels, -Total_burial_depth_abyss./10^(12),'xb--')
    plot(depth_levels, -Total_burial_depth_margin./10^(12),'xr--')
    plot(depth_levels, -Total_burial_depth_shelf./10^(12),'xg--')
    ylim([0 220])
    xlabel('Sediment depth (cm)')
    ylabel('Total burial (Tg C yr-1 = 10^{12}g )')
    hleg=legend('Total', 'Abyss', 'Margin', 'Shelf','Location','best');
    print(fig1,'-depsc2', [par_pathout, '/Global_total_TOC_burial_vs_DEPTH_with10yrs_' str_date '.eps']);
    
  	save(['output/Exp_results/', str_outdir ,'/Total_burial_depth_abyss.mat' ], 'Total_burial_depth_abyss')
  	save(['output/Exp_results/', str_outdir ,'/Total_burial_depth_margin.mat'] , 'Total_burial_depth_margin')
 	save(['output/Exp_results/', str_outdir ,'/Total_burial_depth_shelf.mat'] , 'Total_burial_depth_shelf')
    
    
    fig2=figure;
    set(gca,'FontSize',15)
    hold on
    box on
    plot(age_levels(1:end-1), -Total_burial_age(1:end-1)./10^(12),'xk-')
    plot(age_levels(1:end-1), -Total_burial_age_abyss(1:end-1)./10^(12),'xb--')
    plot(age_levels(1:end-1), -Total_burial_age_margin(1:end-1)./10^(12),'xr--')
    plot(age_levels(1:end-1), -Total_burial_age_shelf(1:end-1)./10^(12),'xg--')
    %caxis([-1e+9 0.0])
    ylim([0 220])
    xlabel('Sediment age (kyrs)')
    ylabel('Total burial (Tg C yr-1 = 10^{12}g )')
    hleg=legend('Total', 'Abyss', 'Margin', 'Shelf','Location','best');
    print(fig2,'-depsc2', ['output/global/', str_outdir ,'/Global_total_TOC_burial_vs_AGE_with10yrs_' str_date '.eps']);
    
   	save(['output/Exp_results/', str_outdir ,'/Total_burial_age_abyss.mat' ], 'Total_burial_age_abyss')
  	save(['output/Exp_results/', str_outdir ,'/Total_burial_age_margin.mat' ], 'Total_burial_age_margin')
 	save(['output/Exp_results/', str_outdir ,'/Total_burial_age_shelf.mat'] , 'Total_burial_age_shelf')
       
    
    calc_rate_of_change = false;
    if(calc_rate_of_change)
    % calc avg rate of change per cm
    Total_burial_depth_RateofChange_diff = diff(Total_burial_depth);
    Total_burial_depth_abyss_RateofChange_diff = diff(Total_burial_depth_abyss);
    Total_burial_depth_margin_RateofChange_diff = diff(Total_burial_depth_margin);
    Total_burial_depth_shelf_RateofChange_diff = diff(Total_burial_depth_shelf);
    
    depth_levels_diff = diff(depth_levels);
    
    for i=1:(length(depth_levels)-1)
        depth_levels_mid(i) = (depth_levels(i+1)-depth_levels(i))/2 + depth_levels(i);
    end
    txt = '\DeltaTOC-burial / \Deltadepth';
    fig12=figure;
    set(gca,'FontSize',15)
    hold on
    box on
    plot(depth_levels_mid, (Total_burial_depth_RateofChange_diff./10^(12))./depth_levels_diff,'xk-')
    plot(depth_levels_mid, (Total_burial_depth_abyss_RateofChange_diff./10^(12))./depth_levels_diff,'xb--')
    plot(depth_levels_mid, (Total_burial_depth_margin_RateofChange_diff./10^(12))./depth_levels_diff,'xr--')
    plot(depth_levels_mid, (Total_burial_depth_shelf_RateofChange_diff./10^(12))./depth_levels_diff,'xg--')
    %caxis([-1e+9 0.0])
    %     ylim([0 120])
    text(200, 2.5, txt)
    xlabel('Sediment depth (cm)')
    ylabel('Avg. Rate of TOC burial Change (per cm)')
    hleg=legend('Total', 'Abyss', 'Margin', 'Shelf','Location','best');
    print(fig12,'-depsc2', ['output/global/', str_outdir ,'/Global_total_TOC_burial_vs_depth_RATEofCHANGE_DIFF_with10yrs_' str_date '.eps']);
    
    
    
    % calc avg rate of change per kyr
    Total_burial_age_RateofChange_diff = diff(Total_burial_age);
    Total_burial_age_abyss_RateofChange_diff = diff(Total_burial_age_abyss);
    Total_burial_age_margin_RateofChange_diff = diff(Total_burial_age_margin);
    Total_burial_age_shelf_RateofChange_diff = diff(Total_burial_age_shelf);
    
    age_levels_diff = diff(age_levels);
    
    for i=1:(length(age_levels)-1)
        age_levels_mid(i) = (age_levels(i+1)-age_levels(i))/2 + age_levels(i);
    end
    txt = '\DeltaTOC-burial / \Deltaage';
    fig22=figure;
    set(gca,'FontSize',15)
    hold on
    box on
    plot(age_levels_mid(1:end-1), (Total_burial_age_RateofChange_diff(1:end-1)./10^(12))./age_levels_diff(1:end-1),'xk-')
    plot(age_levels_mid(1:end-1), (Total_burial_age_abyss_RateofChange_diff(1:end-1)./10^(12))./age_levels_diff(1:end-1),'xb--')
    plot(age_levels_mid(1:end-1), (Total_burial_age_margin_RateofChange_diff(1:end-1)./10^(12))./age_levels_diff(1:end-1),'xr--')
    plot(age_levels_mid(1:end-1), (Total_burial_age_shelf_RateofChange_diff(1:end-1)./10^(12))./age_levels_diff(1:end-1),'xg--')
    %caxis([-1e+9 0.0])
    %     ylim([0 120])
    text(20, 20, txt)
    xlabel('Sediment age (kyrs)')
    ylabel('Avg. Rate of TOC burial Change (per kyr)')
    hleg=legend('Total', 'Abyss', 'Margin', 'Shelf','Location','best');
    print(fig22,'-depsc2', ['output/global/', str_outdir ,'/Global_total_TOC_burial_vs_AGE_100kyr_RATEofCHANGE_DIFF_with10yrs_' str_date '.eps']);
    end
end



if(plot_global_BE_depth)
    
    % calculate Total burial per ocean area:
    if(load_result)
        load(['output/Exp_results/', str_outdir ,'/Flux_depth.mat'])          % spatial map of fluxes in g/yr (fluxes are grid-cell size weighted)
        load(['output/Exp_results/', str_outdir ,'/Flux_age.mat'])            % spatial map of fluxes in g/yr (fluxes are grid-cell size weighted)
        load(['output/Exp_results/', str_outdir ,'/Flux_TOC_swi.mat'])        % spatial map of fluxes in g/yr (fluxes are grid-cell size weighted)
    end
    
    % calculate BE per depth level
    for i=1:length(depth_levels)  % [11,15,20,30,40,50,100,200,300,400,500,1000] cm
        Total_BE_depth(i) = nansum(nansum(Flux_depth{i}))/nansum(nansum(Flux_TOC_swi))*100;  % Global BE (fluxes are grid-cell weighted)
        Total_BE_depth_abyss(i) = nansum(nansum(Flux_depth{i}.*ABYSS_MAP))/nansum(nansum(Flux_TOC_swi.*ABYSS_MAP))*100;
        Total_BE_depth_shelf(i) = nansum(nansum(Flux_depth{i}.*SHELF_MAP))/nansum(nansum(Flux_TOC_swi.*SHELF_MAP))*100;
        Total_BE_depth_margin(i) = nansum(nansum(Flux_depth{i}.*MARGIN_MAP))/nansum(nansum(Flux_TOC_swi.*MARGIN_MAP))*100;
    end
    
    
    % calculate global total burial per age level
    for j=1:length(age_levels)
        Total_BE_age(j) = nansum(nansum(Flux_age{j}))/nansum(nansum(Flux_TOC_swi))*100;
        Total_BE_age_abyss(j) = nansum(nansum(Flux_age{j}.*ABYSS_MAP))/nansum(nansum(Flux_TOC_swi.*ABYSS_MAP))*100;
        Total_BE_age_shelf(j) = nansum(nansum(Flux_age{j}.*SHELF_MAP))/nansum(nansum(Flux_TOC_swi.*SHELF_MAP))*100;
        Total_BE_age_margin(j) = nansum(nansum(Flux_age{j}.*MARGIN_MAP))/nansum(nansum(Flux_TOC_swi.*MARGIN_MAP))*100;
    end
    
    set(0,'defaultLineLineWidth', 2)
    set(0,'DefaultAxesFontSize',4)
    
    fig1=figure;
    set(gca,'FontSize',15)
    hold on
    box on
    plot(depth_levels, Total_BE_depth,'xk-')
    plot(depth_levels, Total_BE_depth_abyss,'xb--')
    plot(depth_levels, Total_BE_depth_margin,'xr--')
    plot(depth_levels, Total_BE_depth_shelf,'xg--')
    %    ylim([0 220])
    xlabel('Sediment depth (cm)')
    ylabel('Burial Efficiency (in %)')
    hleg=legend('Total', 'Abyss', 'Margin', 'Shelf','Location','best');
    print(fig1,'-depsc2', ['output/global/', str_outdir ,'/Global_total_BE_vs_DEPTH_' str_date '.eps']);
    
  	save(['output/Exp_results/', str_outdir ,'/Total_BE_depth.mat' ], 'Total_BE_depth')
  	save(['output/Exp_results/', str_outdir ,'/Total_BE_depth_abyss.mat'] , 'Total_BE_depth_abyss')
  	save(['output/Exp_results/', str_outdir ,'/Total_BE_depth_margin.mat'] , 'Total_BE_depth_margin')
 	save(['output/Exp_results/', str_outdir ,'/Total_BE_depth_shelf.mat'] , 'Total_BE_depth_shelf')

    
    fig2=figure;
    set(gca,'FontSize',15)
    hold on
    box on
    plot(age_levels(1:end-1), Total_BE_age(1:end-1),'xk-')
    plot(age_levels(1:end-1), Total_BE_age_abyss(1:end-1),'xb--')
    plot(age_levels(1:end-1), Total_BE_age_margin(1:end-1),'xr--')
    plot(age_levels(1:end-1), Total_BE_age_shelf(1:end-1),'xg--')
    %caxis([-1e+9 0.0])
    %    ylim([0 220])
    xlabel('Sediment age (kyrs)')
    ylabel('Burial Efficiency (in %)')
    hleg=legend('Total', 'Abyss', 'Margin', 'Shelf','Location','best');
    print(fig2,'-depsc2', ['output/global/', str_outdir ,'/Global_total_BE_vs_AGE_' str_date '.eps']);
    
  	save(['output/Exp_results/', str_outdir ,'/Total_BE_age.mat' ], 'Total_BE_age')
  	save(['output/Exp_results/', str_outdir ,'/Total_BE_age_abyss.mat' ], 'Total_BE_age_abyss')
  	save(['output/Exp_results/', str_outdir ,'/Total_BE_age_margin.mat' ], 'Total_BE_age_margin')
 	save(['output/Exp_results/', str_outdir ,'/Total_BE_age_shelf.mat' ], 'Total_BE_age_shelf')
    
end





if(plot_Bflux_depth_total)

    set(0,'DefaultAxesFontSize',8)

    if(load_result)
        load(['output/Exp_results/', str_outdir ,'/Flux_depth.mat'])
        load(['output/Exp_results/', str_outdir ,'/dxdy.mat'])    % cm^2 per grid-cell
        load(['output/Exp_results/', str_outdir ,'/depth_levels.mat'])
        %      	depth_levels = [11 15 20 30 40 50 100 200 300 400 500 1000];  % in cm
        load(['output/Exp_results/', str_outdir ,'/Total_burial_depth.mat'] , 'Total_burial_depth')
    end
    %
    % fig_SWI_Bflux = figure;
    % set(gca,'FontSize',30)
    % pcolor(long, lat, Flux_TOC_swi);
    % title('TOC flux at the SWI (g/yr)')
    % hold on;
    % shading interp;
    % %contour(long, lat, toc,'LineColor','k')
    % colorbar ('horizontal')
    % colormap(flipud(parula))
    % caxis([-2e+9 0.0])
    % xlabel('Longitude')
    % ylabel('Latitude')
    % print(fig_SWI_Bflux,'-depsc2', ['output/global/Bflux_SWI_AdjustScale_x2_' str_date '.eps']);
    
    
    for i=7 %1:length(depth_levels)
        
        str_title = ['TOC flux at ' num2str(depth_levels(i)) 'cm (in g/(m^2yr))'];
        
        Flux_depth_gm2yr = Flux_depth{i}./dxdy.*100^2;
        
        fig(i) = figure;
        if(new_plotting)
            font_size = 12;

            m_proj('Robinson','longitudes',[-180 179.99], ...
                'latitudes',[-90 90]);
            hold on;
            levels = [0.0:0.5:15.0];
            limits = [0.0 5.0];
            m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
            [C,h] = m_contourf(long, lat,  -Flux_depth_gm2yr, levels);
            set(h,'LineColor','none')
%            title(str_title, 'FontSize', font_size )
            m_coast('linewidth',1,'color','k');
            m_coast('patch',[.7 .7 .7]);
            colorbar ('horizontal')
            colormap(parula)
            caxis(limits)
 %           xlabel('Longitude', 'FontSize', font_size)
 %           ylabel('Latitude', 'FontSize', font_size)
            
            save(['output/Exp_results/', str_outdir ,'/Flux_depth_gm2yr_' num2str(depth_levels(i)) 'cm_depth_.mat'] , 'Flux_depth_gm2yr')
%             print(fig(i),'-depsc2', ['output/global/', str_outdir ,'/Bflux_perm2yr_' num2str(depth_levels(i)) 'cm_depth_' str_date '_until5.eps']);
            print(fig(i),'-dpng',   ['output/global/', str_outdir ,'/Bflux_perm2yr_' num2str(depth_levels(i)) 'cm_depth_' str_date '_until5.png']);  % as .png
            saveas(fig(i),['output/global/', str_outdir ,'/Bflux_perm2yr_' num2str(depth_levels(i)) 'cm_depth_' str_date '_until5.pdf'],'pdf')

        else
            set(gca,'FontSize',30)
            pcolor(long, lat, Flux_depth{i}./dxdy.*100^2);
            title(str_title)
            hold on;
            shading interp;
            %contour(long, lat, toc,'LineColor','k')
            colorbar ('horizontal')
            colormap(flipud(parula))
            caxis([-10 0.0])
            xlabel('Longitude')
            ylabel('Latitude')
            print(fig(i),'-depsc2', ['output/global/', str_outdir ,'/Bflux_perm2yr_' num2str(depth_levels(i)) 'cm_depth_' str_date '_oldProj.eps']);
        end
        Total_burial_depth(i)./10^(12)  % in TgC yr1-
    end
    
end

if(plot_Bflux_age_total)
    
    if(load_result)
        load(['output/Exp_results/', str_outdir ,'/Flux_age.mat'])    % in g/yr (was multiplied with dxdy)
        load(['output/Exp_results/', str_outdir ,'/dxdy.mat'])    % cm^2 per grid-cell
        load(['output/Exp_results/', str_outdir ,'/age_level.mat']);  % in kyr
        %      	age_levels = [0.01 0.1 0.5 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 1000];  % in kyr
        load(['output/Exp_results/', str_outdir ,'/Total_burial_age.mat'] , 'Total_burial_age')
        
    end
    
    for i= [2] %1:length(age_levels)
        
        str_title = ['TOC flux at ' num2str(age_levels(i)) 'kyr (in g/(m^2yr))'];
        
        Flux_age_gm2yr = Flux_age{i}./dxdy.*100^2;
        fig(i) = figure;
        if(new_plotting)
            font_size = 12;

            m_proj('Robinson','longitudes',[-180 179.99], ...
                'latitudes',[-90 90]);
            hold on;
            levels = [0.0:0.5:15.0];
            limits = [0.0 5.0];
            m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
            [C,h] = m_contourf(long, lat, -Flux_age_gm2yr , levels);
            set(h,'LineColor','none')
 %           title(str_title, 'FontSize', font_size )
            m_coast('linewidth',1,'color','k');
            m_coast('patch',[.7 .7 .7]);
            colorbar ('horizontal')
            colormap(parula)
            caxis(limits)
%            xlabel('Longitude', 'FontSize', font_size)
%            ylabel('Latitude', 'FontSize', font_size)
            
        	save(['output/Exp_results/', str_outdir ,'/Flux_age_gm2yr_' num2str(age_levels(i)) 'kyr_age.mat'] , 'Flux_age_gm2yr')
%             print(fig(i),'-depsc2', ['output/global/', str_outdir ,'/Bflux_perm2yr_' num2str(age_levels(i)) 'kyr_age_' str_date '_until5.eps']);
            print(fig(i),'-dpng',   ['output/global/', str_outdir ,'/Bflux_perm2yr_' num2str(age_levels(i)) 'kyr_age_' str_date '_until5.png']);  % as .png
         	saveas(fig(i),['output/global/', str_outdir ,'/Bflux_perm2yr_' num2str(age_levels(i)) 'kyr_age_' str_date '_until5.pdf'],'pdf')

        else
        
        
        set(gca,'FontSize',30)
        pcolor(long, lat, Flux_age{i}./dxdy.*100^2);
        title(str_title)
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(flipud(parula))
        caxis([-10 0.0])
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig(i),'-depsc2', ['output/global/', str_outdir ,'/Bflux_perm2yr_' num2str(age_levels(i)) 'kyr_age_' str_date '.eps']);
        end
        Total_burial_age(i)./10^(12)  % in TgC yr1-
        
    end
    
end

if(plot_Difference_Bflux_age_total)
    
    if(load_result)
        load('output/Exp_results/Flux_age.mat')
    end
    
    for i=2:length(age_levels)
        
        str_title = ['TOC flux at ' num2str(age_levels(i)) 'kyr - 10 years(in g/yr)'];
        
        fig(i) = figure;
        set(gca,'FontSize',30)
        pcolor(long, lat, Flux_age{i}-Flux_age{1});
        title(str_title)
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(parula)
        caxis([0.0 5e+8])
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig(i),'-depsc2', ['output/global/', str_outdir ,'/Difference_Bflux_' num2str(age_levels(i)) 'kyr_to_10yrs_' str_date '.eps']);
        
    end
    
end




if(plot_Bflux_depth_total_difference)
    
    if(load_result)
        load('output/Exp_results/Flux_depth.mat')
        load('output/Exp_results/depth_levels.mat')
    end
    %
    % fig_SWI_Bflux = figure;
    % set(gca,'FontSize',30)
    % pcolor(long, lat, Flux_TOC_swi);
    % title('TOC flux at the SWI (g/yr)')
    % hold on;
    % shading interp;
    % %contour(long, lat, toc,'LineColor','k')
    % colorbar ('horizontal')
    % colormap(flipud(parula))
    % caxis([-2e+9 0.0])
    % xlabel('Longitude')
    % ylabel('Latitude')
    % print(fig_SWI_Bflux,'-depsc2', ['output/global/Bflux_SWI_AdjustScale_x2_' str_date '.eps']);
    
    
    for i=2:length(depth_levels)
        
        str_title = ['TOC flux at ' num2str(depth_levels(i)) 'cm - 11cm (in g/yr)'];
        
        fig(i) = figure;
        set(gca,'FontSize',30)
        pcolor(long, lat, Flux_depth{i}-Flux_depth{1});
        title(str_title)
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(parula)
        caxis([0.0 5e+8])
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig(i),'-depsc2', ['output/global/', str_outdir ,'/Difference_Bflux_' num2str(depth_levels(i)) 'cm_to_11cm_' str_date '.eps']);
        
    end
    
end

if(plot_BE_age)
    font_size = 12;

    if(load_result)
        load(['output/Exp_results/', str_outdir ,'/BE_age.mat'])
        load(['output/Exp_results/', str_outdir ,'/age_level.mat'])
    end    % age_levels = [0.01 0.1 0.5 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 1000];  % in kyr
    
    for i= [2] % 13 18] %1:length(age_levels)
        
        str_title = [num2str(age_levels(i)) ' kyrs'];
        
        fig(i) = figure;
        if(new_plotting)
            m_proj('Robinson','longitudes',[-180 179.99], ...
                'latitudes',[-90 90]);
            hold on;
            levels = [0:0.5:15.0];
            limits = [0 10.0];
            m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
            [C,h] = m_contourf(long, lat,  BE_age{i}, levels);
            set(h,'LineColor','none')
            title(str_title, 'FontSize', font_size )
            m_coast('linewidth',1,'color','k');
            m_coast('patch',[.7 .7 .7]);
            colorbar ('horizontal')
            colormap(parula)
            caxis(limits)
%             xlabel('Longitude', 'FontSize', font_size)
%             ylabel('Latitude', 'FontSize', font_size)
            print(fig(i),'-depsc2', ['output/global/', str_outdir ,'/BE_' num2str(age_levels(i)) 'kyr_age_' str_date '.eps']);  % as .eps
            print(fig(i),'-dpng', ['output/global/', str_outdir ,'/BE_' num2str(age_levels(i)) 'kyr_age_' str_date '.png' ]);  % as .png
            saveas(fig(i),['output/global/', str_outdir ,'/BE_' num2str(age_levels(i)) 'kyr_age_' str_date '.pdf'],'pdf')

        else
            set(gca,'FontSize',30)
            pcolor(long, lat, BE_age{i});
            title(str_title)
            hold on;
            shading interp;
            %contour(long, lat, toc,'LineColor','k')
            colorbar ('horizontal')
            colormap(parula)
            caxis([0 10.0])
            xlabel('Longitude')
            ylabel('Latitude')
            print(fig(i),'-depsc2', ['output/global/', str_outdir ,'/BE_' num2str(age_levels(i)) 'kyr_age_' str_date '_oldProj.eps']);
        end
        
    end
    
end


if(plot_BE_age_difference)
    
    
    % age_levels = [0.01 0.1 0.5 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 1000];  % in kyr
    
    for i=2:length(age_levels)
        
        str_title = ['BE ' num2str(age_levels(i)) 'kyrs - 10 yrs (in %)'];
        
        fig(i) = figure;
        set(gca,'FontSize',30)
        pcolor(long, lat, BE_age{i}-BE_age{1});
        title(str_title)
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(flipud(parula))
        caxis([-5 0.0])
        xlabel('Longitude')
        ylabel('Latitude')
        print(fig(i),'-depsc2', ['output/global/', str_outdir ,'/Difference_BE_' num2str(age_levels(i)) 'kyr_to_10yrs_' str_date '.eps']);
        print(fig(i),'-dpng', ['output/global/', str_outdir ,'/Difference_BE_' num2str(age_levels(i)) 'kyr_to_10yrs_' str_date '.png']);
       	saveas(fig(i),['output/global/', str_outdir ,'/Difference_BE_' num2str(age_levels(i)) 'kyr_to_10yrs_' str_date '.pdf'],'pdf')

    end
    
end

if(plot_BE_depth)
    
    font_size = 12;
    if(load_result)        
        load(['output/Exp_results/', str_outdir ,'/BE_depth.mat'])
        BE = BE_depth;
        load(['output/Exp_results/', str_outdir ,'/depth_levels.mat'])
    end
    %           	depth_levels = [11 15 20 30 40 50 100 200 300 400 500 1000];  % in cm
    BE = BE_depth;
    
    
    if(false)
        fig_11cm = figure;
        set(gca,'FontSize',30)
        pcolor(long, lat, BE{1});
        % pcolor(long(1:end-1), lat(1:end-1), BE);
        title('Burial efficiency at 11cm')
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(parula)
        caxis([0 10.0])
        xlabel('Longitude')
        ylabel('Latitude')
        % File_text = ['Tansect_BE_@_' num2str(age_horizon/1000) 'kyrs_fix_a'];
        print(fig_11cm,'-depsc2', ['output/global/', str_outdir ,'/BE_11cm_depth_' str_date '.eps']);
        
        
        fig_20cm = figure;
        set(gca,'FontSize',30)
        pcolor(long, lat, BE{3});
        % pcolor(long(1:end-1), lat(1:end-1), BE);
        title('Burial efficiency at 20cm')
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(parula)
        caxis([0 10.0])
        xlabel('Longitude')
        ylabel('Latitude')
        
        % File_text = ['Tansect_BE_@_' num2str(age_horizon/1000) 'kyrs_fix_a'];
        print(fig_20cm,'-depsc2', ['output/global/', str_outdir ,'/BE_20cm_depth_' str_date '.eps']);
        
        fig_50cm = figure;
        set(gca,'FontSize',30)
        pcolor(long, lat, BE{6});
        % pcolor(long(1:end-1), lat(1:end-1), BE);
        title('Burial efficiency at 50cm')
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(parula)
        caxis([0 10.0])
        xlabel('Longitude')
        ylabel('Latitude')
        
        % File_text = ['Tansect_BE_@_' num2str(age_horizon/1000) 'kyrs_fix_a'];
        print(fig_50cm,'-depsc2', ['output/global/', str_outdir ,'/BE_50cm_depth_' str_date '.eps']);
        
        
        fig_5m = figure;
        set(gca,'FontSize',30)
        pcolor(long, lat, BE{11});
        % pcolor(long(1:end-1), lat(1:end-1), BE);
        title('Burial efficiency at 5m')
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(parula)
        caxis([0 10.0])
        xlabel('Longitude')
        ylabel('Latitude')
        % File_text = ['Tansect_BE_@_' num2str(age_horizon/1000) 'kyrs_fix_a'];
        print(fig_5m,'-depsc2', ['output/global/', str_outdir ,'/BE_5m_depth_' str_date '.eps']);
        
    end
    

    
    fig_1m = figure;
    if(new_plotting)
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:15.0];
        limits = [0 10.0];
     	m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
        [C,h] = m_contourf(long, lat,  BE{7}, levels);
        set(h,'LineColor','none')
        title('1 mbsf', 'FontSize', font_size)
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
%        m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
%         xlabel('Longitude', 'FontSize', font_size)
%         ylabel('Latitude', 'FontSize', font_size)
        %   print(fig_1m,'-depsc2',
        %   ['output/global/', str_outdir ,'/BE_1m_depth_' str_date '.eps']);  % save as .eps
         print(fig_1m,'-depsc2', ['output/global/', str_outdir ,'/BE_1m_depth_' str_date '.eps']);  % as .eps
        print(fig_1m,'-dpng', ['output/global/', str_outdir ,'/BE_1m_depth_' str_date '.png' ]);  % save as .png
       	saveas(fig_1m,['output/global/', str_outdir ,'/BE_1m_depth_' str_date '.pdf'],'pdf')

    else
        
        set(gca,'FontSize',30)
        pcolor(long, lat, BE{7});
        % pcolor(long(1:end-1), lat(1:end-1), BE);
        title('Burial efficiency at 1m')
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(parula)
        caxis([0 10.0])
        xlabel('Longitude')
        ylabel('Latitude')
        
        % File_text = ['Tansect_BE_@_' num2str(age_horizon/1000) 'kyrs_fix_a'];
        print(fig_1m,'-depsc2', ['output/global/', str_outdir ,'/BE_1m_depth_' str_date '_oldProj.eps']);
    end
    
if false      
    fig_20cm = figure;
    if(new_plotting)
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:15.0];
        limits = [0 10.0];
        m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
        [C,h] = m_contourf(long, lat,  BE{3}, levels);
        set(h,'LineColor','none')
        title('20 cmbsf', 'FontSize', font_size)
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
%         xlabel('Longitude', 'FontSize', font_size)
%         ylabel('Latitude', 'FontSize', font_size)
        %   print(fig_1m,'-depsc2',
        %   ['output/global/', str_outdir ,'/BE_1m_depth_' str_date '.eps']);  % save as .eps
     	print(fig_20cm,'-depsc2', ['output/global/', str_outdir ,'/BE_20cm_depth_' str_date '.eps']);  % as .eps
        print(fig_20cm,'-dpng', ['output/global/', str_outdir ,'/BE_20cm_depth_' str_date ]);  % save as .png
        
    else
        
    end   
    fig_10m = figure;
    if(new_plotting)
        m_proj('Robinson','longitudes',[-180 179.99], ...
            'latitudes',[-90 90]);
        hold on;
        levels = [0:0.5:15.0];
        limits = [0 10.0];
        m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
        [C,h] = m_contourf(long, lat,  BE{12}, levels);
        set(h,'LineColor','none')
        title('10 mbsf', 'FontSize', font_size)
        m_coast('linewidth',1,'color','k');
        m_coast('patch',[.7 .7 .7]);
%         m_grid('fontsize',8);
        colorbar ('horizontal')
        colormap(parula)
        caxis(limits)
%         xlabel('Longitude', 'FontSize', font_size)
%         ylabel('Latitude', 'FontSize', font_size)
        %   print(fig_10m,'-depsc2',
        %   ['output/global/', str_outdir ,'/BE_10m_depth_' str_date '.eps']); % save as .eps
       print(fig_1m,'-depsc2', ['output/global/', str_outdir ,'/BE_10m_depth_' str_date '.eps']);  % as .eps
        print(fig_10m,'-dpng', ['output/global/', str_outdir ,'/BE_10m_depth_' str_date ]);  % save as .png
        
    else
        set(gca,'FontSize',30)
        pcolor(long, lat, BE{12});
        % pcolor(long(1:end-1), lat(1:end-1), BE);
        title('Burial efficiency at 10m')
        hold on;
        shading interp;
        %contour(long, lat, toc,'LineColor','k')
        colorbar ('horizontal')
        colormap(parula)
        caxis([0 10.0])
        xlabel('Longitude')
        ylabel('Latitude')
        
        % File_text = ['Tansect_BE_@_' num2str(age_horizon/1000) 'kyrs_fix_a'];
        print(fig_10m,'-depsc2', ['output/global/', str_outdir ,'/BE_10m_depth_' str_date '_oldProj.eps']);
    end
end
end


if(plot_BE_depth_difference)
    %
    if(load_result)
        
        load('output/Exp_results/BE_depth.mat')
        BE = BE_depth;
        load('output/Exp_results/depth_levels.mat')
    end
    %  depth_levels = [11 15 20 30 40 50 100 200 300 400 500 1000];  % in cm
    
    
    
    fig_20cm = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, BE{3}-BE{1});
    % pcolor(long(1:end-1), lat(1:end-1), BE);
    title('BE 20cm - 11cm')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula) )
    caxis([-5 0])
    xlabel('Longitude')
    ylabel('Latitude')
    
    % File_text = ['Tansect_BE_@_' num2str(age_horizon/1000) 'kyrs_fix_a'];
    print(fig_20cm,'-depsc2', ['output/global/', str_outdir ,'/Diff_to11cm_BE_20cm_depth_' str_date '.eps']);
    
    fig_50cm = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, BE{6}-BE{1});
    % pcolor(long(1:end-1), lat(1:end-1), BE);
    title('BE 50cm - 11cm')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula) )
    caxis([-5 0])
    xlabel('Longitude')
    ylabel('Latitude')
    
    % File_text = ['Tansect_BE_@_' num2str(age_horizon/1000) 'kyrs_fix_a'];
    print(fig_50cm,'-depsc2', ['output/global/', str_outdir ,'/Diff_to11cm_BE_50cm_depth_' str_date '.eps']);
    
    
    fig_1m = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, BE{7}-BE{1});
    % pcolor(long(1:end-1), lat(1:end-1), BE);
    title('BE 1m - 11cm')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula) )
    caxis([-5 0])
    xlabel('Longitude')
    ylabel('Latitude')
    
    % File_text = ['Tansect_BE_@_' num2str(age_horizon/1000) 'kyrs_fix_a'];
    print(fig_1m,'-depsc2', ['output/global/', str_outdir ,'/Diff_to11cm_BE_1m_depth_' str_date '.eps']);
    
    
    fig_10m = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, BE{12}-BE{1});
    % pcolor(long(1:end-1), lat(1:end-1), BE);
    title('BE 10m - 11cm')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula) )
    caxis([-5 0])
    xlabel('Longitude')
    ylabel('Latitude')
    
    % File_text = ['Tansect_BE_@_' num2str(age_horizon/1000) 'kyrs_fix_a'];
    print(fig_10m,'-depsc2', ['output/global/', str_outdir ,'/Diff_to11cm_BE_10m_depth_' str_date '.eps']);
end


if(plot_BC)
    load zholo.mat
    load zpleisto.mat
    load 'age.mat'
    
    toc(isnan(zholo))=NaN;
    
    fig_age= figure;
    set(gca,'FontSize',30)
    pcolor(long(1:end-1), lat(1:end-1), age_zz);
    title('Age.mat (years?)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    %caxis([0 0.001])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_age,'-depsc2', ['output/global/Age_' str_date '.eps']);
    
    
    fig_zholo= figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, zholo);
    title('Depth of Holocene layer (cm below zbio)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    %caxis([0 0.001])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_zholo,'-depsc2', ['output/global/zHolo_' str_date '.eps']);
    
    fig_zpleisto= figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, zpleisto);
    title('Depth of Pleistocene layer (cm below zbio)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    %caxis([0 0.001])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_zpleisto,'-depsc2', ['output/global/zPleisto_' str_date '.eps']);
    
    fig_SWI_TOC = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, toc);
    title('TOC (wt%)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    caxis([0 5.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_TOC,'-depsc2', ['output/global/SWI_TOC_' str_date '.eps']);
    
    fig_dxdy = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, dxdy);
    title('dx*dy')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(parula)
    %    caxis([0 5.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_dxdy,'-depsc2', ['output/global/DxDy_' str_date '.eps']);
end

if(plot_Bflux_depth_perArea)
    
    fig_SWI_Bflux = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, Flux.TOC_swi);
    title('TOC flux at the SWI (g/(cm^2 yr)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula))
    caxis([-10e-4 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_SWI_Bflux,'-depsc2', ['output/global/Bflux_preArea_SWI_AdjustScale_x10_' str_date '.eps']);
    
    
    fig_Holocene_Bflux = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, Flux.TOC_Holo);
    title('TOC flux at the Holocene (g/(cm^2 yr)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula))
    caxis([-1e+9 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_Holocene_Bflux,'-depsc2', ['output/global/Bflux_preArea_Holocene_AdjustScale_' str_date '.eps']);
    
    
    fig_10cm_Bflux = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, Flux.TOC_10cm);
    title('TOC flux at 10cm (g/(cm^2 yr)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula))
    caxis([-1e+9 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_10cm_Bflux,'-depsc2', ['output/global/Bflux_preArea_10cm_depth_' str_date '.eps']);
    
    
    fig_20cm_Bflux = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, Flux.TOC_20cm);
    title('TOC flux at 20cm (g/(cm^2 yr)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula))
    caxis([-1e+9 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_20cm_Bflux,'-depsc2', ['output/global/Bflux_preArea_20cm_depth_' str_date '.eps']);
    
    
    fig_50cm_Bflux = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, Flux.TOC_50cm);
    title('TOC flux at 50cm (g/(cm^2 yr)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula))
    caxis([-1e+9 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_50cm_Bflux,'-depsc2', ['output/global/Bflux_preArea_50cm_depth_' str_date '.eps']);
    
    
    fig_1m_Bflux = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, Flux.TOC_1m);
    title('TOC flux at 1m (g/(cm^2 yr)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula))
    caxis([-1e+9 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_1m_Bflux,'-depsc2', ['output/global/Bflux_preArea_1m_depth_' str_date '.eps']);
    
    
    fig_5m_Bflux = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, Flux.TOC_5m);
    title('TOC flux at 5m (g/(cm^2 yr)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula))
    caxis([-1e+9 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_5m_Bflux,'-depsc2', ['output/global/Bflux_preArea_5m_depth_' str_date '.eps']);
    
    
    fig_10m_Bflux = figure;
    set(gca,'FontSize',30)
    pcolor(long, lat, Flux.TOC_10m);
    title('TOC flux at 10m (g/(cm^2 yr)')
    hold on;
    shading interp;
    %contour(long, lat, toc,'LineColor','k')
    colorbar ('horizontal')
    colormap(flipud(parula))
    caxis([-1e+9 0.0])
    xlabel('Longitude')
    ylabel('Latitude')
    print(fig_10m_Bflux,'-depsc2', ['output/global/Bflux_preArea_10m_depth_' str_date '.eps']);
    
end
