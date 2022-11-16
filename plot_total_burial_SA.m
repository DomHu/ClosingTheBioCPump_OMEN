%% Plotting script for the results of the global simulations
% plot toc burial per ocean region with uncertainty intervals
% specify were different SA results are saved under "SET MODEL RESUTS FOLDERS"
use_Lee = true;

load_result = true;  % true: load results from existing runs; false: if results are still in workspace variables

N = 3; % use three different colors
C = linspecer(N);

depth_levels = [11 15 20 30 40 50 100 200 300 400 500 1000];  % in cm
age_levels = [0.01 0.1 0.5 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 1000];  % in kyr

%  load 'toc.dat'
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
    

str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];


%% SET MODEL RESUTS FOLDERS
str_df_burial = 'SpatialPor_Lee_wdepth_combinedSA_1.0_por';

% change por and sedrate
str_low_burial = 'SpatialPor_Lee_wdepth_combinedSA_1.1_por';
str_high_burial = 'SpatialPor_Lee_wdepth_combinedSA_0.9_por';
% change porosity only
% str_low_burial = '0909_BWO2_SpatialPor_Lee_wdepth_SA_por1.1';
% str_high_burial = '0909_BWO2_SpatialPor_Lee_wdepth_SA_por0.9';
% change sedrate only
% str_low_burial = '1708_SpatialPor_Lee_wdepth_SA_sedrate0.9';
% str_high_burial = '1708_SpatialPor_Lee_wdepth_SA_sedrate1.1';

load(['output/Exp_results/', str_df_burial, '/depth_levels.mat'])
load(['output/Exp_results/', str_df_burial, '/age_level.mat'])

if true
%% TOTAL TOC BURIAL

if(load_result)
    
    % load default
    Total_burial_depth_df = load(['output/Exp_results/', str_df_burial '/Total_burial_depth.mat']).Total_burial_depth;    
    Total_burial_depth_abyss_df = load(['output/Exp_results/', str_df_burial '/Total_burial_depth_abyss.mat' ], 'Total_burial_depth_abyss').Total_burial_depth_abyss;
    Total_burial_depth_margin_df = load(['output/Exp_results/', str_df_burial '/Total_burial_depth_margin.mat'] , 'Total_burial_depth_margin').Total_burial_depth_margin;
    Total_burial_depth_shelf_df = load(['output/Exp_results/', str_df_burial '/Total_burial_depth_shelf.mat'] , 'Total_burial_depth_shelf').Total_burial_depth_shelf;
    
    Total_burial_age_df = load(['output/Exp_results/', str_df_burial '/Total_burial_age.mat'] , 'Total_burial_age').Total_burial_age;
    Total_burial_age_abyss_df = load(['output/Exp_results/', str_df_burial '/Total_burial_age_abyss.mat'] , 'Total_burial_age_abyss').Total_burial_age_abyss;
    Total_burial_age_margin_df = load(['output/Exp_results/', str_df_burial '/Total_burial_age_margin.mat'] , 'Total_burial_age_margin').Total_burial_age_margin;
    Total_burial_age_shelf_df = load(['output/Exp_results/', str_df_burial '/Total_burial_age_shelf.mat'] , 'Total_burial_age_shelf').Total_burial_age_shelf;
    
    % load low
    Total_burial_depth_low = load(['output/Exp_results/', str_low_burial '/Total_burial_depth.mat']).Total_burial_depth;    
    Total_burial_depth_abyss_low = load(['output/Exp_results/', str_low_burial '/Total_burial_depth_abyss.mat' ], 'Total_burial_depth_abyss').Total_burial_depth_abyss;
    Total_burial_depth_margin_low = load(['output/Exp_results/', str_low_burial '/Total_burial_depth_margin.mat'] , 'Total_burial_depth_margin').Total_burial_depth_margin;
    Total_burial_depth_shelf_low = load(['output/Exp_results/', str_low_burial '/Total_burial_depth_shelf.mat'] , 'Total_burial_depth_shelf').Total_burial_depth_shelf;
    
    Total_burial_age_low = load(['output/Exp_results/', str_low_burial '/Total_burial_age.mat'] , 'Total_burial_age').Total_burial_age;
    Total_burial_age_abyss_low = load(['output/Exp_results/', str_low_burial '/Total_burial_age_abyss.mat'] , 'Total_burial_age_abyss').Total_burial_age_abyss;
    Total_burial_age_margin_low = load(['output/Exp_results/', str_low_burial '/Total_burial_age_margin.mat'] , 'Total_burial_age_margin').Total_burial_age_margin;
    Total_burial_age_shelf_low = load(['output/Exp_results/', str_low_burial '/Total_burial_age_shelf.mat'] , 'Total_burial_age_shelf').Total_burial_age_shelf;
    
    % load high
    Total_burial_depth_high = load(['output/Exp_results/', str_high_burial '/Total_burial_depth.mat']).Total_burial_depth;    
    Total_burial_depth_abyss_high = load(['output/Exp_results/', str_high_burial '/Total_burial_depth_abyss.mat' ], 'Total_burial_depth_abyss').Total_burial_depth_abyss;
    Total_burial_depth_margin_high = load(['output/Exp_results/', str_high_burial '/Total_burial_depth_margin.mat'] , 'Total_burial_depth_margin').Total_burial_depth_margin;
    Total_burial_depth_shelf_high = load(['output/Exp_results/', str_high_burial '/Total_burial_depth_shelf.mat'] , 'Total_burial_depth_shelf').Total_burial_depth_shelf;
    
    Total_burial_age_high = load(['output/Exp_results/', str_high_burial '/Total_burial_age.mat'] , 'Total_burial_age').Total_burial_age;
    Total_burial_age_abyss_high = load(['output/Exp_results/', str_high_burial '/Total_burial_age_abyss.mat'] , 'Total_burial_age_abyss').Total_burial_age_abyss;
    Total_burial_age_margin_high = load(['output/Exp_results/', str_high_burial '/Total_burial_age_margin.mat'] , 'Total_burial_age_margin').Total_burial_age_margin;
    Total_burial_age_shelf_high = load(['output/Exp_results/', str_high_burial '/Total_burial_age_shelf.mat'] , 'Total_burial_age_shelf').Total_burial_age_shelf;
    
end


set(0,'defaultLineLineWidth', 2)
set(0,'DefaultAxesFontSize',4)

fig1=figure;
set(gca,'FontSize',15)
hold on
box on
fill([depth_levels';flipud(depth_levels')],[-Total_burial_depth_low'./10^(12);flipud(-Total_burial_depth_high')./10^(12)],[.75 .75 .75],'linestyle','none');
fill([depth_levels';flipud(depth_levels')],[-Total_burial_depth_abyss_low'./10^(12);flipud(-Total_burial_depth_abyss_high')./10^(12)],[.75 .75 .75],'linestyle','none');
fill([depth_levels';flipud(depth_levels')],[-Total_burial_depth_margin_low'./10^(12);flipud(-Total_burial_depth_margin_high')./10^(12)],[.75 .75 .75],'linestyle','none');
fill([depth_levels';flipud(depth_levels')],[-Total_burial_depth_shelf_low'./10^(12);flipud(-Total_burial_depth_shelf_high')./10^(12)],[.75 .75 .75],'linestyle','none');

h(1) = plot(depth_levels, -Total_burial_depth_df./10^(12),'xk-');
h(2) = plot(depth_levels, -Total_burial_depth_abyss_df./10^(12),'x--', 'color', C(1,:));
h(3) = plot(depth_levels, -Total_burial_depth_margin_df./10^(12),'x--', 'color', C(2,:));
h(4) = plot(depth_levels, -Total_burial_depth_shelf_df./10^(12),'x--', 'color', C(3,:));
alpha(.6)

ylim([0 200])
xlabel('Sediment depth (cm)')
ylabel('Total burial (Tg C yr^{-1})')
hleg=legend(h([1 2 3 4]),'Total', 'Abyss', 'Margin', 'Shelf','Location','best');
% print(fig1,'-depsc2', ['output/global/', str_outdir ,'/Global_total_TOC_burial_vs_DEPTH_Unvertainty_POR_' str_date '.eps']);  % as .eps
% print(fig1,'-dpng', ['output/global/', str_outdir ,'/Global_total_TOC_burial_vs_DEPTH_Unvertainty_POR_' str_date '.png']);
saveas(fig1,['output/global/Global_total_TOC_burial_vs_DEPTH_Unvertainty_POR_SEDRATE_' str_date '.pdf'],'pdf')


fig2=figure;
set(gca,'FontSize',15)
hold on
box on
fill([age_levels';flipud(age_levels')],[-Total_burial_age_low'./10^(12);flipud(-Total_burial_age_high')./10^(12)],[.75 .75 .75],'linestyle','none');
fill([age_levels';flipud(age_levels')],[-Total_burial_age_abyss_low'./10^(12);flipud(-Total_burial_age_abyss_high')./10^(12)],[.75 .75 .75],'linestyle','none');
fill([age_levels';flipud(age_levels')],[-Total_burial_age_margin_low'./10^(12);flipud(-Total_burial_age_margin_high')./10^(12)],[.75 .75 .75],'linestyle','none');
fill([age_levels';flipud(age_levels')],[-Total_burial_age_shelf_low'./10^(12);flipud(-Total_burial_age_shelf_high')./10^(12)],[.75 .75 .75],'linestyle','none');

h(1) = plot(age_levels(1:end-1), -Total_burial_age_df(1:end-1)./10^(12),'xk-');
h(2) = plot(age_levels(1:end-1), -Total_burial_age_abyss_df(1:end-1)./10^(12),'x--', 'color', C(1,:));
h(3) = plot(age_levels(1:end-1), -Total_burial_age_margin_df(1:end-1)./10^(12),'x--', 'color', C(2,:));
h(4) = plot(age_levels(1:end-1), -Total_burial_age_shelf_df(1:end-1)./10^(12),'x--', 'color', C(3,:));
alpha(.6)
%caxis([-1e+9 0.0])
xlim([0 100])
ylim([0 200])
xlabel('Sediment age (kyrs)')
ylabel('Total burial (Tg C yr^{-1})')
hleg=legend(h([1 2 3 4]),'Total', 'Abyss', 'Margin', 'Shelf','Location','best');
saveas(fig2, ['output/global/Global_total_TOC_burial_vs_AGE_Unvertainty_POR_SEDRATE_' str_date '.pdf'],'pdf');

end

%% BE vs DEPTH


    % load default
    Total_BE_depth_df = load(['output/Exp_results/', str_df_burial '/Total_BE_depth.mat']).Total_BE_depth;    
    Total_BE_depth_abyss_df = load(['output/Exp_results/', str_df_burial '/Total_BE_depth_abyss.mat' ], 'Total_BE_depth_abyss').Total_BE_depth_abyss;
    Total_BE_depth_margin_df = load(['output/Exp_results/', str_df_burial '/Total_BE_depth_margin.mat'] , 'Total_BE_depth_margin').Total_BE_depth_margin;
    Total_BE_depth_shelf_df = load(['output/Exp_results/', str_df_burial '/Total_BE_depth_shelf.mat'] , 'Total_BE_depth_shelf').Total_BE_depth_shelf;
    
    Total_BE_age_df = load(['output/Exp_results/', str_df_burial '/Total_BE_age.mat'] , 'Total_BE_age').Total_BE_age;
    Total_BE_age_abyss_df = load(['output/Exp_results/', str_df_burial '/Total_BE_age_abyss.mat'] , 'Total_BE_age_abyss').Total_BE_age_abyss;
    Total_BE_age_margin_df = load(['output/Exp_results/', str_df_burial '/Total_BE_age_margin.mat'] , 'Total_BE_age_margin').Total_BE_age_margin;
    Total_BE_age_shelf_df = load(['output/Exp_results/', str_df_burial '/Total_BE_age_shelf.mat'] , 'Total_BE_age_shelf').Total_BE_age_shelf;

    % load low
    Total_BE_depth_low = load(['output/Exp_results/', str_low_burial '/Total_BE_depth.mat']).Total_BE_depth;    
    Total_BE_depth_abyss_low = load(['output/Exp_results/', str_low_burial '/Total_BE_depth_abyss.mat' ], 'Total_BE_depth_abyss').Total_BE_depth_abyss;
    Total_BE_depth_margin_low = load(['output/Exp_results/', str_low_burial '/Total_BE_depth_margin.mat'] , 'Total_BE_depth_margin').Total_BE_depth_margin;
    Total_BE_depth_shelf_low = load(['output/Exp_results/', str_low_burial '/Total_BE_depth_shelf.mat'] , 'Total_BE_depth_shelf').Total_BE_depth_shelf;
    
    Total_BE_age_low = load(['output/Exp_results/', str_low_burial '/Total_BE_age.mat'] , 'Total_BE_age').Total_BE_age;
    Total_BE_age_abyss_low = load(['output/Exp_results/', str_low_burial '/Total_BE_age_abyss.mat'] , 'Total_BE_age_abyss').Total_BE_age_abyss;
    Total_BE_age_margin_low = load(['output/Exp_results/', str_low_burial '/Total_BE_age_margin.mat'] , 'Total_BE_age_margin').Total_BE_age_margin;
    Total_BE_age_shelf_low = load(['output/Exp_results/', str_low_burial '/Total_BE_age_shelf.mat'] , 'Total_BE_age_shelf').Total_BE_age_shelf;
    
    % load high
    Total_BE_depth_high = load(['output/Exp_results/', str_high_burial '/Total_BE_depth.mat']).Total_BE_depth;    
    Total_BE_depth_abyss_high = load(['output/Exp_results/', str_high_burial '/Total_BE_depth_abyss.mat' ], 'Total_BE_depth_abyss').Total_BE_depth_abyss;
    Total_BE_depth_margin_high = load(['output/Exp_results/', str_high_burial '/Total_BE_depth_margin.mat'] , 'Total_BE_depth_margin').Total_BE_depth_margin;
    Total_BE_depth_shelf_high = load(['output/Exp_results/', str_high_burial '/Total_BE_depth_shelf.mat'] , 'Total_BE_depth_shelf').Total_BE_depth_shelf;
    
    Total_BE_age_high = load(['output/Exp_results/', str_high_burial '/Total_BE_age.mat'] , 'Total_BE_age').Total_BE_age;
    Total_BE_age_abyss_high = load(['output/Exp_results/', str_high_burial '/Total_BE_age_abyss.mat'] , 'Total_BE_age_abyss').Total_BE_age_abyss;
    Total_BE_age_margin_high = load(['output/Exp_results/', str_high_burial '/Total_BE_age_margin.mat'] , 'Total_BE_age_margin').Total_BE_age_margin;
    Total_BE_age_shelf_high = load(['output/Exp_results/', str_high_burial '/Total_BE_age_shelf.mat'] , 'Total_BE_age_shelf').Total_BE_age_shelf;
    
    
    set(0,'defaultLineLineWidth', 2)
    set(0,'DefaultAxesFontSize',4)
    
    fig3=figure;
    set(gca,'FontSize',15)
    hold on
    box on
    fill([depth_levels';flipud(depth_levels')],[Total_BE_depth_low';flipud(Total_BE_depth_high')],[.75 .75 .75],'linestyle','none');
    fill([depth_levels';flipud(depth_levels')],[Total_BE_depth_abyss_low';flipud(Total_BE_depth_abyss_high')],[.75 .75 .75],'linestyle','none');
    fill([depth_levels';flipud(depth_levels')],[Total_BE_depth_margin_low';flipud(Total_BE_depth_margin_high')],[.75 .75 .75],'linestyle','none');
    fill([depth_levels';flipud(depth_levels')],[Total_BE_depth_shelf_low';flipud(Total_BE_depth_shelf_high')],[.75 .75 .75],'linestyle','none');
    
    h(1) = plot(depth_levels, Total_BE_depth_df,'xk-');
    h(2) = plot(depth_levels, Total_BE_depth_abyss_df,'x--', 'color', C(1,:));
    h(3) = plot(depth_levels, Total_BE_depth_margin_df,'x--', 'color', C(2,:));
    h(4) = plot(depth_levels, Total_BE_depth_shelf_df,'x--', 'color', C(3,:));
    alpha(.6)
    ylim([0 10])
    xlabel('Sediment depth (cm)')
    ylabel('Burial or transfer efficiency (%)')
    hleg=legend(h([1 2 3 4]),'Total', 'Abyss', 'Margin', 'Shelf','Location','best');
    saveas(fig3, ['output/global/Global_total_BE_vs_DEPTH_Unvertainty_POR_SEDRATE_' str_date '.pdf'], 'pdf');
%     saveas(fig3, ['output/global/', str_outdir ,'/Global_total_BE_vs_DEPTH_Uncertainty_SEDRATE_ONLY_' str_date '.png'], 'png');
 
        

    fig4=figure;
    set(gca,'FontSize',15)
    hold on
    box on
    fill([age_levels';flipud(age_levels')],[Total_BE_age_low';flipud(Total_BE_age_high')],[.75 .75 .75],'linestyle','none');
    fill([age_levels';flipud(age_levels')],[Total_BE_age_abyss_low';flipud(Total_BE_age_abyss_high')],[.75 .75 .75],'linestyle','none');
    fill([age_levels';flipud(age_levels')],[Total_BE_age_margin_low';flipud(Total_BE_age_margin_high')],[.75 .75 .75],'linestyle','none');
    fill([age_levels';flipud(age_levels')],[Total_BE_age_shelf_low';flipud(Total_BE_age_shelf_high')],[.75 .75 .75],'linestyle','none');

    h(1) = plot(age_levels(1:end-1), Total_BE_age_df(1:end-1),'xk-');
    h(2) = plot(age_levels(1:end-1), Total_BE_age_abyss_df(1:end-1),'x--', 'color', C(1,:));
    h(3) = plot(age_levels(1:end-1), Total_BE_age_margin_df(1:end-1),'x--', 'color', C(2,:));
    h(4) = plot(age_levels(1:end-1), Total_BE_age_shelf_df(1:end-1),'x--', 'color', C(3,:));
    alpha(.6)
    %caxis([-1e+9 0.0])
    xlim([0 100])
    ylim([0 10])
    xlabel('Sediment age (kyrs)')
    ylabel('Burial or transfer efficiency (%)')
    hleg=legend(h([1 2 3 4]),'Total', 'Abyss', 'Margin', 'Shelf','Location','best');
    saveas(fig4, ['output/global/Global_total_BE_vs_AGE_Unvertainty_POR_SEDRATE_' str_date '.pdf'], 'pdf');
