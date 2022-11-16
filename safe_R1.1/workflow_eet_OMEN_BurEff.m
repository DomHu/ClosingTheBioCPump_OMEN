% This script provides the application of the SAFE toolbox to calculate the Elementary Effects Test presented in the SI. 
% Get the entire toolbox from https://www.safetoolbox.info/ and copy this file in the folder 
% Please contact me (dominik.huelse@gmx.de) directly if you want to .    
%                                                                                 
% METHOD                                                                       
%                                                                              
% This script provides an example of application of the Elementary Effects     
% Test (EET) or 'method of Morris' (Morris, 1991; Saltelli et al., 2008).      
%                                                                              
% The EET is a One-At-the-Time method for global Sensitivity Analysis.         
% It computes two indices for each input:                                      
% i) the mean (mi) of the EEs, which measures the total effect of an input     
% over the output;                                                             
% ii) the standard deviation (sigma) of the EEs, which measures the degree     
% of interactions with the other inputs.                                       
% Both sensitivity indices are relative measures, i.e. their value does not    
% have any specific meaning per se but it can only be used in pair-wise        
% comparison (e.g. if input x(1) has higher mean EEs than input x(3) than      
% x(1) is more influential than x(3)).                                         
%                                                                              
% For an application example in the environmental domain, see for instance     
% Nguyen and de Kok (2007).                                                    
%                                                                              
% MODEL AND STUDY AREA                                                         
%                                                                              
% The model under study is the rainfall-runoff model Hymod                     
% (see help of function hymod_sim.m for more details)                          
% applied to the Leaf catchment in Mississipi, USA                             
% (Sorooshian et al., 1983).                                                   
% The inputs subject to SA are the 5 model parameters, and the scalar          
% output for SA is a metric of model performance.                              
%                                                                              
% INDEX                                                                        
%                                                                              
% Steps:                                                                       
% 1. Add paths to required directories                                         
% 2. Load data and set-up the HBV model                                        
% 3. Sample inputs space                                                       
% 4. Run the model against input samples                                       
% 5. Compute the elementary effects                                            
% 6. Example of how to repeat computions after adding up new                   
%    input/output samples.                                                     
%                                                                              
% REFERENCES                                                                   
%                                                                              
% Morris, M.D. (1991), Factorial sampling plans for preliminary                
% computational experiments, Technometrics, 33(2).                             
%                                                                              
% Nguyen, T.G. and de Kok, J.L. (2007). Systematic testing of an integrated    
% systems model for coastal zone management using sensitivity and              
% uncertainty analyses. Env. Mod. & Soft., 22, 1572-1587.                      
%                                                                              
% Saltelli, A., et al. (2008) Global Sensitivity Analysis, The Primer,         
% Wiley.                                                                       
%                                                                              
% Sorooshian, S., Gupta, V., Fulton, J. (1983). Evaluation of maximum          
% likelihood parameter estimation techniques for conceptual rainfall-runoff    
% models: Influence of calibration data variability and length on model        
% credibility. Water Resour. Res., 19, 251-259.                                
                                                                               
% This script prepared by Francesca Pianosi and Fanny Sarrazin                 
% University of Bristol, 2014                                                  
% mail to: francesca.pianosi@bristol.ac.uk     

plot_MS_results = true;
                                                                               
%% Step 1 (add paths)                                                          
my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab    
% current directory to the SAFE directory. Otherwise, you may define           
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:            
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';                      
                                                                               
% Set current directory to 'my_dir' and add path to sub-folders:               
cd(my_dir)                                                              
addpath(genpath(my_dir))
addpath([ my_dir '/..'])
                                                                 
%% Step 2: setup the OMEN model
plot_combined = true;
% Number of uncertain parameters subject to SA:                                
M    = 4 ;                                                                     
% Parameter ranges (from literature):                                          
% sed-rate (factor), Dbio (factor), zbio ( value [1, 15]), porostity (factor)
xmin = [  0.8  0.5 1 0.8 ];                                                     
xmax = [  1.2  1.5 20 1.2];                 
% xmin = [  0.5  1 ];                                                     
% xmax = [  1.5 15];
    % Parameter distributions:                                                     
DistrFun  = 'unif'  ;                                                          
DistrPar = cell(M,1);                                                          
for i=1:M; DistrPar{i} = [ xmin(i) xmax(i) ] ; end                             
% Name of parameters (will be used to costumize plots):                        
X_labels = {'SedRate','Dbio','zbio','por'} ;                                    
                                                                               
% Define output:                                                               
myfun = 'OMEN_BurEff_SA' ;                                                          
   
Titles = {'12', '15', '20', '50', '100', '1000'};  % Depth levels in cm
out = length(Titles); % number of output from OMEN (e.g. SWI-flux O2, NO3, ...)

%% Step 3 (sample inputs space)                                                
                                                                               
r = 50 ; % Number of Elementary Effects                                       
% [notice that the final number of model evaluations will be equal to          
% r*(M+1)]                                                                     
                                                                               
% option 1: use the sampling method originally proposed by Morris (1991):      
% L = 6  ; % number of levels in the uniform grid                              
% design_type  = 'trajectory'; % (note used here but required later)           
% X = Morris_sampling(r,xmin,xmax,L); % (r*(M+1),M)                            
                                                                               
% option 2: Latin Hypercube sampling strategy                                  
SampStrategy = 'lhs' ; % Latin Hypercube                                       
design_type = 'radial';                                                        
% other options for design type:                                               
%design_type  = 'trajectory';        
if plot_MS_results
    load('X_1409_r50.mat')
    load('Y_1409_r50.mat')
%     load('X_1309.mat')
%     load('Y_1309.mat')
else
X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_type);              
                                                                               
%% Step 4 (run the model)                                                      
Y = model_evaluation(myfun,X) ; % size (r*(M+1),1)   
save('X_1409_r50.mat','X')
save('Y_1409_r50.mat', 'Y')
end
%% Step 5 (Computation of the Elementary effects)                              

for depths=3:out    % 11cm and 12 cm lead EEs = Inf
    j=depths-2;
% % Compute Elementary Effects:                                                  
% [ mi, sigma ] = EET_indices(r,xmin,xmax,X,Y(:,depths),design_type);                      
                                                                               
% % % Plot results in the plane (mean(EE),std(EE)):                                
% % EET_plot(mi, sigma,X_labels )                                                  
                                                                               
% Use bootstrapping to derive confidence bounds:                               
Nboot=100;                                           
    [mi,sigma,EE,mi_sd,sigma_sd,mi_lb,sigma_lb,mi_ub,sigma_ub] = ...               
    EET_indices(r,xmin,xmax,X,Y(:,depths),design_type,Nboot);   
    if(plot_combined)
        % also save in 3D vector for combined plotting
        mi_mul{j} = mi;
        sigma_mul{j}= sigma; 
%        EE{j}=EE;
        mi_sd_mul{j}=mi_sd;
        sigma_sd_mul{j}=sigma_sd;
        mi_lb_mul{j} = mi_lb;
        sigma_lb_mul{j}=sigma_lb;    
        mi_ub_mul{j} = mi_ub;
        sigma_ub_mul{j} = sigma_ub;
    end
    % Plot bootstrapping results in the plane (mean(EE),std(EE)):                  
%    EET_plot(mi,sigma,X_labels,mi_lb,mi_ub,sigma_lb,sigma_ub, depths)     
end

if(plot_combined)
  	EET_plot_multiple(mi_mul,sigma_mul,X_labels,mi_lb_mul,mi_ub_mul,sigma_lb_mul,sigma_ub_mul)     
end

%% NOT USED

if(false)   % not necessary
% Repeat computations using a decreasing number of samples so as to assess     
% if convergence was reached within the available dataset:                     
rr = [ r/5:r/5:r ] ;                                                           
m_r = EET_convergence(EE,rr);                                                  
% Plot the sensitivity measure (mean of elementary effects) as a function      
% of model evaluations:                                                        
figure; plot_convergence(m_r,rr*(M+1),[],[],[],...                             
'no of model evaluations','mean of EEs',X_labels)                              
                                                                               
% Repeat convergence analysis using bootstrapping:                             
Nboot = 100;                                                                   
rr = [ r/5:r/5:r ] ;                                                           
[m_r,s_r,m_lb_r,m_ub_r] = EET_convergence(EE,rr,Nboot);                        
% Plot the sensitivity measure (mean of elementary effects) as a function      
% of model evaluations:                                                        
figure; plot_convergence(m_r,rr*(M+1),m_lb_r,m_ub_r,[],...                     
'no of model evaluations','mean of EEs',X_labels)                              
                                                                               
%% Step 6 (Adding up new samples)                                              
                                                                               
r2 = 200 ; % increase of base sample size                                      
[X2,Xnew] = OAT_sampling_extend(X,r2,DistrFun,DistrPar,design_type);% extended 
% sample (it includes the already evaluated sample 'X' and the new one)        
                                                                               
% Evaluate model against the new sample                                        
Ynew = model_evaluation(myfun,Xnew,rain,evap,flow) ; % size((r2-r)*(M+1),1)    
                                                                               
% Put new and old results together                                             
Y2=[Y;Ynew]; % size (r2*(M+1),1)                                               
                                                                               
% Recompute indices                                                            
Nboot=100;                                                                     
[mi_n,sigma_n,EEn,mi_sdn,sigma_sdn,mi_lbn,sigma_lbn,mi_ubn,sigma_ubn] = ...      
EET_indices(r2,xmin,xmax,X2,Y2,design_type,Nboot);                             
EET_plot(mi_n,sigma_n,X_labels,mi_lbn,mi_ubn,sigma_lbn,sigma_ubn)                
                                                                               
% Repeat convergence analysis                                                  
Nboot = 100;                                                                   
rr2 = [ r2/5:r2/5:r2 ] ;                                                       
[m_rn,s_rn,m_lb_rn,m_ub_rn] = EET_convergence(EEn,rr2,Nboot);                  
% Plot the sensitivity measure (mean of elementary effects) as a function      
% of model evaluations:                                                        
figure; plot_convergence(m_rn,rr2*(M+1),m_lb_rn,m_ub_rn,[],...                 
'no of model evaluations','mean of EEs',X_labels)                              
end