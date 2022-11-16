function [Y_out] = OMEN_BurEff_SA(x)
%
% This is called from ../safe_R1.1/workflow_eet_OMEN_BurEff
% This function runs the OMEN-SED sediment model
% and returns the associated Burial flux of TOC and BE

%% Value of uncertain parameter
sedrate_in = x(1);
Dbio_in = x(2);
zbio_in = x(3);
por_in = x(4);


%% Calculate the BE using TOC observations from Seiter and BC as in Bradley ea. 2020

warning('query','all')
warning('off','all')
warning
%__________________________________________________________________________
%   load data
%__________________________________________________________________________
addpath('../data/Lee_et_al_2019/')
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
load('../data/RECCAP2/bathymetry_matrix_new_ud.mat')
water_depth_updated = -bathymetry_matrix_new_ud;


load('../data/RECCAP2/porosity_matrix_new_ud.mat')
porosity = porosity_matrix_new_ud;


% specify depth-levels to calculate BE and so on
depth_levels = [12 15 20 50 100 1000];  % in cm

% specify size of results
BE_depth = cell(length(depth_levels), 1);
Flux_depth = cell(length(depth_levels), 1);

% specify age-levels to calculate BE and so on
% this would give similar results than the SA for depth-levels


zbio = zbio_in;

[m,n]=size(toc);

zmax_sed_holo = max(max(zholo))+zbio;     % define maximum depth of sediment layer (in cm)
zmax_sed_holo_round = ceil(zmax_sed_holo/100)*100;  % round to nearest meter

ncl = 1;
res.bsd = benthic_main(ncl);
res.bsd.zbio = zbio_in;
res.bsd.zinf = zmax_sed_holo_round;


res.bsd.usescalarcode = ncl==1;

swi = benthic_test.default_swi();
swi.TwoG_OM_model = false;
swi.flux = false;
swi.IntConst_GMD = true;
%__________________________________________________________________________
%   set specific parameters
%__________________________________________________________________________

res.swi = swi;

xstart=1; %55;
xstop=m;  %55;        %lat
ystart=1;  % 6;
ystop=n;  % 6;        %long (or the other way around- who knows!)


for x = 1:m
    %         	for x = xstart:xstop
%     x
    % calculate volume
    % convert deg to cm CODAS package (by E.Firing,et al.)
    rlat = lat(x) * pi/180;
    m = 111132.09  - 566.05 * cos(2 * rlat)+ 1.2 * cos(4 * rlat);
    dy = 0.25*m*100.0; %cm
    p = 111415.13 * cos(rlat) - 94.55 * cos(3 * rlat);
    dx = 0.25*p*100.0; %cm
    
    for y = 1:n
        %                for y = ystart:ystop
        %                   
        
        if ((isnan(toc(x,y))))    % now check for sed_holo or toc = NaN
            
            toc(x,y)=NaN;
            dxdy(x,y)   = NaN;
            
            
            for i=1:length(depth_levels)
                BE_depth{i}(x,y) =NaN;
                Flux_depth{i}(x,y)=NaN;
            end
            
        else
            
            res.swi.p_a=0.1*SHELF_MAP(x,y)+1.0*MARGIN_MAP(x,y)+20.0*ABYSS_MAP(x,y);    %reactive continuum model a=10?-1-10?5
            %         res.swi.p_nu =0.125*SHELF_MAP(i,j)+0.125*MARGIN_MAP(i,j)+0.125*ABYSS_MAP(i,j); %reactive continuum model nu=0.01-1, mostly around 0.125
            %         loga=3.35-14.81*sed_holo(i,j);                          %sedimentation rate dependent formulation for a (Arndt et al., 2013)
            %         swi.p_a=10.^(loga);
            if(~isnan(water_depth_updated(x,y)))
                res.bsd.Dbio=benthic_main.biorate(water_depth_updated(x,y))*Dbio_in;
            else
                
                res.bsd.Dbio=(27.6*SHELF_MAP(x,y)+8.1*MARGIN_MAP(x,y)+0.59*ABYSS_MAP(x,y))*Dbio_in;  % Bioturbation coefficient
            end
                if(~isnan(porosity(x,y)))
                    res.bsd.por=por_in*porosity(x,y)/100;
                    if(res.bsd.por>0.95)
                        res.bsd.por=0.95;
                    end
                else
                    res.bsd.por=por_in*(0.599*SHELF_MAP(x,y)+0.695*MARGIN_MAP(x,y)+0.755*ABYSS_MAP(x,y));   % use mean value
                end

            if(~isnan(sed_holo(x,y)))
                res.bsd.w=sedrate_in*sed_holo(x,y);       % 2 - factor, so 0.9 and 1.1 is  swaped as higher sed-rate causes more burial
            else    % use default
                res.bsd.w=sedrate_in*benthic_main.sedrate(water_depth_updated(x,y));
            end
            
            %                        beta=0.5e-5*SHELF_MAP(x,y)+1.7e-5*MARGIN_MAP(x,y)+0.85e-5*ABYSS_MAP(x,y);
            
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
            
            
            %%%%%%%%%%%%%%%%%
            % through different depth levels - calculate maps of BE and Fluxes
            for i=1:length(depth_levels)
                %  depth_levels = [11 15 20 30 40 50 100 200 300 400 500 1000];  % in cm
                [F_TOC_inf, F_TOC1_inf] = res.zTOC_RCM.calcCflx(depth_levels(i), res.bsd, res.swi, res);
                Flux_depth{i}(x,y) = nansum(F_TOC1_inf)*dx*dy;  %   in g/yr
                % Flux_depth_perArea{i}(x,y) = nansum(F_TOC1_inf);  %   in g /(cm^2 yr)
                BE_depth{i}(x,y) = nansum(F_TOC1_inf)/nansum(F_TOC1_swi)*100;
            end
            
        end        
    end
end

% calculate global total burial per depth level
for i=1:length(depth_levels)
    Total_burial_depth(i) = nansum(nansum(Flux_depth{i}));
    Y_out(i) = nansum(nansum(Flux_depth{i}))/10^(12);  % in Pg C /yr
end


end


