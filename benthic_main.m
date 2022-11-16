% OMEN-SED 1.0 BENTHIC-MODEL Stand-alone matlab code
% Hülse et al (2017) GMD paper

% benthic_main.m
% Global properties for benthic model

classdef benthic_main < handle
    % Global properties for benthic model
    
    properties
        
        ncl;                                    % number of sediment columns
        usescalarcode = true;                   % use scalar code
        
        tol_const = 1e-18;                      % non-zero constant to avoid numerical errors (e.g. division by zero)
        
        %sediment characteristics
        rho_sed=2.5;                            % sediment density (g/cm3)
        wdepth=3000.0;                           % water depth (m)
        w;                                      % burial velocity  (cm/yr) - calculated by internal fct. sedrate()
        z0  = 0;                                % surface
        zbio=10.0;                              % bioturbation depth (cm)
        zinf=100;                               %Inifinity (cm)
        
        Dbio;                                   % bioturbation coefficient (cm2/yr) - calculated by internal fct. biorate()
        por=0.85;                               % porosity (-)
        tort=3.0;                               %tortuosity (-)
        irrigationFactor=1.0;                   %irrigation factor (-)
        dispFactor;                             %dispersion factor (-)
        
        %stoichiometric factors
        X_C;                                    % Carbon Redfield stoichiometry
        Y_N;                                    % Nitrogen Redfield stoichiometry
        Z_P;                                    % Phosphorous Redfield stoichiometry
        SD;                                     % volume factor solid->dissolved phase
        
        
        zoxgf = 0.0;                            % cm, rolloff NH4, H2S oxidation for small zox depth (was 0.1)
        
        % Diagnostic output from root finder
        %fzerooptions;
        %fzerooptions = optimset('Display','iter');
        %fzerooptions = optimset('Display','final');
        %fzerooptions = optimset('TolX',0.001);
        fzerooptions = optimset('TolX',100*eps);
    end
    
    methods
        function obj = benthic_main(ncl, wdepth)
            % set default values for the sediment columns
            if nargin > 0
                obj.ncl = ncl;
                obj.wdepth = obj.wdepth*ones(1,obj.ncl);
                obj.zbio = obj.zbio*ones(1,obj.ncl);
                obj.Dbio = obj.Dbio*ones(1,obj.ncl);
                obj.zinf = obj.zinf*ones(1,obj.ncl);
                obj.z0 = obj.z0*ones(1,obj.ncl);
            else
                obj.ncl = 1;
            end
            
            if nargin > 1
                obj.wdepth = wdepth;
            end
            
            obj.usescalarcode = (obj.ncl == 1);
            
            obj.w=benthic_main.sedrate(obj.wdepth);
            obj.Dbio=benthic_main.biorate(obj.wdepth);
            obj.dispFactor=obj.por.^(obj.tort-1.0).*obj.irrigationFactor;	%dispersion factor (-)
 %           obj.SD=(1-obj.por)./obj.por;
            
            obj.X_C=106;                        % Carbon Redfield stoichiometry
            obj.Y_N=16;                        	% Nitrogen Redfield stoichiometry
            obj.Z_P=1;                        	% Phosphorous Redfield stoichiometry

        end
    end
    
    methods(Static)
        
        function w = sedrate(wdepth)
            % sedimentation rate, cm/yr (after Middelburg et al. (1997))
            w_Middelburg = 10.0.^(-0.87478367-0.00043512*wdepth)*3.3;
            % sedimentation rate, cm/yr (after Burwicz et al. (2011))
            w1 = 0.117;
            w2 = 0.006;
            z1 = 200;
            z2 = 4000;
            c1 = 3;
            c2 = 10;
            w = w1/(1+(wdepth/z1)^c1) + w2/(1+(wdepth/z2)^c2);
            %            w = 6.92856570653099e-003;  % just fix it for now
        end
        
        function Dbio = biorate(wdepth)
            % bioturbation coeff, cm^2/yr (after Middelburg et al. (1997))
            Dbio= 5.2*(10.0^(0.7624-0.0003972*wdepth));
            %            Dbio=30.0;
        end
        
        function a = apparent_age(w)
            % Continuum model parameter, apparent initial age , yr (after Boudreau and Ruddick, 1991)
            %            a= 4970*exp(-0.0296/1000*w);
            % This one is from Sandra's review paper:
            a= 10^(3.35-14.81*w);
        end
        
        function a = apparent_age_Boudreau(w)
            % Continuum model parameter, apparent initial age , yr (after Boudreau and Ruddick, 1991)
            a= 4970*exp(-0.0296*w);
        end
    end
    
end

