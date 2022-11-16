%% plot the boundary conditions for OMEN-SED

% need to intall M_Map for plotting https://www.eoas.ubc.ca/~rich/map.html
path(path,'/home/domhu/Documents/MATLAB/M_Map');

plot_BC = true;                     % plot boundary connditions?
plot_depth_range_for_SA = true;     % plot porosity & Dbio values for sensitivity analysis?

str_date = [datestr(date,7), datestr(date,5), datestr(date,11)];


if(plot_BC)
    use_Lee = true;
    
    set(0,'defaultLineLineWidth', 2)
    set(0,'DefaultAxesFontSize',8)
    
    
    %__________________________________________________________________________
    %   load data
    %__________________________________________________________________________
    
    addpath('./data/Lee_et_al_2019/')
    load('Lee_toc_lr_weighted.mat')  % Lee data mean weighted by grid-size, meanhas NaN for terrestial cells
    toc = Lee_toc_lr_weighted;
    load('lat_lr.mat')
    lat = lat_lr;
    load('long_lr.mat')
    long = long_lr;
    load sed_holo.mat
    sed_holo = sed_holo(1:end-1, 1:end-1);  % delete extra row and column
    load ABYSS_MAP_Lee.mat
    ABYSS_MAP = ABYSS_MAP_Lee;
    load SHELF_MAP_Lee.mat
    SHELF_MAP = SHELF_MAP_Lee;  % delete extra row and column
    load MARGIN_MAP_Lee.mat
    MARGIN_MAP = MARGIN_MAP_Lee;  % delete extra row and column
    load zholo.mat
    
    load('water_depth_updated_Lee.mat')     % water-depth NASA
    water_depth_updated = -water_depth_updated_Lee;
    
    load('./data/RECCAP2/bathymetry_matrix_new_ud.mat') % from GEBCO (https://www.gebco.net/)
    
    load('./data/RECCAP2/porosity_matrix_new_ud.mat')
    porosity = porosity_matrix_new_ud;
    
    load('./data/O2_BW_WOA2018_hr_Qdegr.mat');   % BW  O2 [muM = 10^-6 mol/kg]  -- need to convert into mol/cm^3 , i.e. *10^-3
    
    
    
    %% calculate Dbio:
    % bioturbation coeff, cm^2/yr (after Middelburg et al. (1997))
    
    [m,n]=size(water_depth_updated);
    for x=1:m
        for y=1:n
            Dbio(x,y) = 5.2*(10.0^(0.7624-0.0003972.*(water_depth_updated(x,y))));
        end
    end
    
    load('lat_lr.mat')
    lat = lat_lr;
    load('long_lr.mat')
    long = long_lr;
    
    %% plot Dbio
    fig(99) = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:0.5:30];
    limits = [0 25];
    m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
    [C,h] = m_contourf(long, lat,  Dbio, levels);
    set(h,'LineColor','none')
    title('Dbio (cm^2/yr)','fontsize',14) % -- Lee et al. (2019)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    print(fig(99),'-dpng', ['./09_Mask_Dbio.png']);
    print(fig(99),'-depsc2', ['./09_Mask_Dbio.eps']);
    
    
    %% plot TOC
    fig(1) = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:0.1:5.0];
    limits = [0 3.0];
    m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
    [C,h] = m_contourf(long, lat,  toc, levels);
    set(h,'LineColor','none')
    title('TOC (wt%)','fontsize',14) % -- Lee et al. (2019)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    print(fig(1),'-dpng', ['./01_Mask_SWI_TOC_Lee_' str_date '.png']);
    print(fig(1),'-depsc2', ['./01_Mask_SWI_TOC_Lee_' str_date '.eps']);
    
    
    %% sedimentation rate
    fig(2) = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:0.001:0.13];
    limits = [0 0.12];
    m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
    [C,h] = m_contourf(long, lat,  sed_holo, levels);
    set(h,'LineColor','none')
    title('Sedimentation rate (cm/yr)','fontsize',14)
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    print(fig(2),'-dpng', ['./02_Mask_Sed_rate_' str_date '.png']);
    print(fig(2),'-depsc2', ['./02_Mask_Sed_rate_' str_date '.eps']);
    
    
    
    %% porosity
    fig(3) = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [25:1:100];
    limits = [40 90];
    m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
    [C,h] = m_contourf(long, lat,  porosity, levels);
    set(h,'LineColor','none')
    title('Porosity (%)','fontsize',14) % -- Lee et al. (2019)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    print(fig(3),'-dpng', ['./03_Mask_Porosity_' str_date '.png']);
    print(fig(3),'-depsc2', ['./03_Mask_Porosity_' str_date '.eps']);
    
    
    
    %% water depth from GEBCO (https://www.gebco.net/)
    fig(6) = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:0.05:6.0];
    limits = [0 5.0];
    m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
    [C,h] = m_contourf(long, lat,  -bathymetry_matrix_new_ud/1000, levels);
    set(h,'LineColor','none')
    title('Water depth (km)','fontsize',14) % -- Lee et al. (2019)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    colorbar ('horizontal')
    colormap(flipud(parula))
    caxis(limits)
    print(fig(6),'-dpng', ['./04_Mask_bathymetry_RECCAP2_' str_date '.png']);
    print(fig(6),'-depsc2', ['./04_Mask_bathymetry_RECCAP2_' str_date '.eps']);
    
    
    %% Bottom water O2
    fig(7) = figure;
    m_proj('Robinson','longitudes',[-180 179.99], ...
        'latitudes',[-90 90]);
    hold on;
    levels = [0:2:350];
    limits = [0 350.0];
    m_grid('fontsize',8,'xticklabels',[],'yticklabels',[]);
    [C,h] = m_contourf(long, lat,  O2_BW_WOA2018_hr_Qdegr, levels);
    set(h,'LineColor','none')
    title('Seafloor Oxygen (\muM)','fontsize',14) % -- Lee et al. (2019)')
    m_coast('linewidth',1,'color','k');
    m_coast('patch',[.7 .7 .7]);
    colorbar ('horizontal')
    colormap(parula)
    caxis(limits)
    print(fig(7),'-dpng', ['./05_Mask_BW_O2_' str_date '.png']);
    print(fig(7),'-depsc2', ['./05_Mask_BW_O2_' str_date '.eps']);
    
    
end

if(plot_depth_range_for_SA)
    
    %% porosity vs depth
    %   	set(0,'defaultLineLineWidth', 2)
    %     set(0,'DefaultAxesFontSize',14)
    %
    formatSpec = '%.2f';
    
    
    load('./data/RECCAP2/porosity_matrix_new_ud.mat')
    porosity_matrix_new = porosity_matrix_new_ud;
    load('./data/RECCAP2/bathymetry_matrix_new_ud.mat')     % from GEBCO (https://www.gebco.net/)
    
    water_depth_updated = bathymetry_matrix_new_ud;
    
    [m,n]=size(water_depth_updated);
    
    % calculate Dbio:
    % bioturbation coeff, cm^2/yr (after Middelburg et al. (1997))
    for x=1:m
        for y=1:n
            Dbio(x,y) = 5.2*(10.0^(0.7624-0.0003972.*(-water_depth_updated(x,y))));
        end
    end
    
    % fewer depth bins as for our NatComms manuscript
    depth_bins_new = [-200 -3500];
    
    %% porosity
    mean_por_local(1)=nanmean(porosity_matrix_new(water_depth_updated>=depth_bins_new(1)));
    mean_por_local(2)=nanmean(porosity_matrix_new(water_depth_updated>=depth_bins_new(2) & water_depth_updated<depth_bins_new(1)));
    mean_por_local(3)=nanmean(porosity_matrix_new(water_depth_updated<depth_bins_new(2)));
    
    std_por_local(1)=nanstd(porosity_matrix_new(water_depth_updated>=depth_bins_new(1)));
    std_por_local(2)=nanstd(porosity_matrix_new(water_depth_updated>=depth_bins_new(2) & water_depth_updated<depth_bins_new(1)));
    std_por_local(3)=nanstd(porosity_matrix_new(water_depth_updated<depth_bins_new(2)));
    
    Bins = length(mean_por_local);
    
    fig_depth_por = figure;
    plot(porosity_matrix_new(:), water_depth_updated(:)/1000,'kx')
    xlim([30 100])
    title('Porosity vs Depth')
    xlabel('Porosity (%)')
    ylabel('Water Depth (km)')
    
    formatSpec = '%.1f';
    
    
    text(32,-6 , ['Shelf: Depth > ', num2str(depth_bins_new(1)/1000), 'km: '], 'FontSize', 8)
    text(32,-6.5, ['Mean = ', num2str(mean_por_local(1),formatSpec) , '%; Std = ', num2str(std_por_local(1),formatSpec)], 'FontSize', 8)
    
    text(32,-7.5 , ['Margins: Depth > ', num2str(depth_bins_new(2)/1000), 'km: '], 'FontSize', 8)
    text(32,-8, ['Mean = ', num2str(mean_por_local(2),formatSpec) , '%; Std = ', num2str(std_por_local(2),formatSpec)], 'FontSize', 8)
    
    text(32,-9, ['Abyss:  Depth <= ', num2str(depth_bins_new(Bins-1)/1000), 'km:'], 'FontSize', 8)
    text(32,-9.5, ['Mean = ', num2str(mean_por_local(Bins),formatSpec) , '%; Std = ', num2str(std_por_local(Bins),formatSpec)], 'FontSize', 8)
    print(fig_depth_por,'-dpng', ['./101_DepthRange_Porosity_fewerbins_RECCAP2_bathym.png']);
    
    
    
    %% Dbio
    mean_Dbio_local(1)=nanmean(Dbio(water_depth_updated>=depth_bins_new(1)));
    mean_Dbio_local(2)=nanmean(Dbio(water_depth_updated>=depth_bins_new(2) & water_depth_updated<depth_bins_new(1)));
    mean_Dbio_local(3)=nanmean(Dbio(water_depth_updated<depth_bins_new(2)));
    
    std_Dbio_local(1)=nanstd(Dbio(water_depth_updated>=depth_bins_new(1)));
    std_Dbio_local(2)=nanstd(Dbio(water_depth_updated>=depth_bins_new(2) & water_depth_updated<depth_bins_new(1)));
    std_Dbio_local(3)=nanstd(Dbio(water_depth_updated<depth_bins_new(2)));
    
    Bins = length(mean_Dbio_local);
    
    fig_depth_Dbio = figure;
    plot(Dbio(:), water_depth_updated(:)/1000,'kx')
    title('Dbio vs Depth')
    xlabel('Dbio (cm^2/yr)')
    ylabel('Water Depth (km)')
    
    formatSpec = '%.2f';
    
    
    text(5,-6 , ['Shelf: Depth > ', num2str(depth_bins_new(1)/1000), 'km: '], 'FontSize', 8)
    text(5,-6.5, ['Mean = ', num2str(mean_Dbio_local(1),formatSpec) , 'cm^2/yr; Std = ', num2str(std_Dbio_local(1),formatSpec), 'cm^2/yr'], 'FontSize', 8)
    
    text(5,-7.5 , ['Margins: Depth > ', num2str(depth_bins_new(2)/1000), 'km: '], 'FontSize', 8)
    text(5,-8, ['Mean = ', num2str(mean_Dbio_local(2),formatSpec) , 'cm^2/yr; Std = ', num2str(std_Dbio_local(2),formatSpec), 'cm^2/yr'], 'FontSize', 8)
    
    text(5,-9, ['Abyss:  Depth <= ', num2str(depth_bins_new(Bins-1)/1000), 'km:'], 'FontSize', 8)
    text(5,-9.5, ['Mean = ', num2str(mean_Dbio_local(Bins),formatSpec) , 'cm^2/yr; Std = ', num2str(std_Dbio_local(Bins),formatSpec), 'cm^2/yr'], 'FontSize', 8)
    print(fig_depth_Dbio,'-dpng', ['./102_DepthRange_Dbio_fewerbins_RECCAP2_bathym.png']);
    
end