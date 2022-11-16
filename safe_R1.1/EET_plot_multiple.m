function EET_plot_multiple(mi_all,sigma_all,varargin)
%
% Plot the sensitivity indices computed by the Elementary Effects Test -
% mean (mi) of the EEs on the horizontal axis and standard deviation (sigma)
% of the EEs on the vertical axis for multiple EEs
% (see help of EET_indices for more details about the EET and references)
%
% Usage:
%
% EET_plot(mi,sigma)                               
% EET_plot(mi,sigma,labelinput)
% EET_plot(mi,sigma,labelinput,mi_lb,mi_ub,sigma_lb,sigma_ub)    
%
% NOTE: The input variables are 3D vectors now
%         mi = mean of the elementary effects               - vector (1,M)
%      sigma = standard deviation of the elementary effects - vector (1,M)
% labelinput = strings for the x-axis labels            - cell array (1,M)
%      mi_lb = lower bound of 'mi'                          - vector (1,M)
%      mi_ub = upper bound of 'mi'                          - vector (1,M)
%   sigma_lb = lower bound of 'sigma'                       - vector (1,M)
%   sigma_ub = upper bound of 'sigma'                       - vector (1,M)

% This function is part of the SAFE Toolbox by F. Pianosi, F. Sarrazin 
% and T. Wagener at Bristol University (2015). 
% SAFE is provided without any warranty and for non-commercial use only. 
% For more details, see the Licence file included in the root directory 
% of this distribution.
% For any comment and feedback, or to discuss a Licence agreement for 
% commercial use, please contact: francesca.pianosi@bristol.ac.uk
% For details on how to cite SAFE in your publication, please see: 
% bristol.ac.uk/cabot/resources/safe-toolbox/

% Options for the graphic:
fn = 'Helvetica' ; % font type of axes, labels, etc.
%fn = 'Courier' ;
fs = 20 ; % font size of axes, labels, etc.
ms = 14 ; % marker size

% Options for the legend:
sorting   = 1  ; % If 1, inputs will be displayed in the legend 
% according to their influence, i.e. from most sensitive to least sensitive
% (if 0 they will be displayed according to their original order)
nb_legend = 5  ; % number of input names that will be displayed in the legend

% Options for the colours:
% You can produce a coloured plot or a black and white one
% (printer-friendly). Furthermore, you can use matlab colourmaps or
% repeat 5 'easy-to-distinguish' colours (see http://colorbrewer2.org/).
% Option 1a - coloured using colorbrewer: uncomment the following line:
col = [[228,26,28];[55,126,184];[77,175,74];[152,78,163];[255,127,0]]/256;  cc = 'k' ;
% Option 1b - coloured using matlab colormap: uncomment the following line:
%col=hsv(length(mi));   cc = 'k' ;
% Option 1a - B&W using matlab colorbrewer: uncomment the following line:
%col = [[37 37 37];[90 90 90];[150 150 150];[189 189 189];[217 217 217]]/256; cc = 'w' ;
% Option 1b - B&W using matlab colormap: uncomment the following line:
%col=gray(length(mi)); cc = 'w' ; 


%%%%%%%%%%%%%%
% Check inputs
%%%%%%%%%%%%%%

if ~isnumeric(mi_all{1}); error('''mi'' must be a vector of size (1,M)'); end
if ~isnumeric(sigma_all{1}); error('''sigma'' must be a vector of size (1,M)'); end
[N,M] = size(mi_all{1})  ;
[n,m] = size(sigma_all{1}) ;
[depths]= length(mi_all);    % how many different depths?

if N~=1; error('''mi'' must be a row vector'); end
if n~=1; error('''sigma'' must be a row vector'); end
if M~=m; error('''mi'' and''sigma'' must have the same number of elements'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover and check optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Set optional arguments to their default values:
labelinput = cell(1,M); for i=1:M; labelinput{i}=['X' num2str(i)]; end
mi_lb    = zeros(1,M) ; 
mi_ub    = zeros(1,M) ; 
sigma_lb = zeros(1,M) ;
sigma_ub = zeros(1,M) ;

% OPEN figure
figure
hold on


% Recover and update optional arguments:
if nargin > 2
    if ~isempty(varargin{1})
        labelinput = varargin{1};
        if ~iscell(labelinput); error('''labelinput'' must be a cell array'); end
        if length(labelinput)~=M; error('''labelinput'' must have M=%d components',M); end
        for i=1:M; if ~ischar(labelinput{i}); error('all components of ''labelinput'' must be string'); end; end
    end
end
for k=1:depths
    
    mi = mi_all{k};
    sigma = sigma_all{k};
if nargin > 3
    if ~isempty(varargin{2})
        mi_lb_all = varargin{2} ;
        mi_lb=mi_lb_all{k};
        if ~isnumeric(mi_lb); error('''mi_lb'' must be a vector of size (1,M)'); end
        [n,m] = size(mi_lb) ;
        if n~=1; error('''mi_lb'' must be a row vector'); end
        if M~=m; error('''mi'' and''mi_lb'' must have the same number of elements'); end
    end
end
if nargin > 4
    if ~isempty(varargin{3})
        mi_ub_all = varargin{3} ;
        mi_ub=mi_ub_all{k};
        if ~isnumeric(mi_ub); error('''mi_ub'' must be a vector of size (1,M)'); end
        [n,m] = size(mi_ub) ;
        if n~=1; error('''mi_ub'' must be a row vector'); end
        if M~=m; error('''mi'' and''mi_ub'' must have the same number of elements'); end
    end
end
if nargin > 4
    if ~isempty(varargin{4})
        sigma_lb_all = varargin{4} ;
        sigma_lb = sigma_lb_all{k};
        if ~isnumeric(sigma_lb); error('''sigma_lb'' must be a vector of size (1,M)'); end
        [n,m] = size(sigma_lb) ;
        if n~=1; error('''sigma_lb'' must be a row vector'); end
        if M~=m; error('''sigma'' and''sigma_lb'' must have the same number of elements'); end
    end
end
if nargin > 5
    if ~isempty(varargin{5})
        sigma_ub_all = varargin{5} ;
        sigma_ub = sigma_ub_all{k};
        if ~isnumeric(sigma_ub); error('''sigma_ub'' must be a vector of size (1,M)'); end
        [n,m] = size(sigma_ub) ;
        if n~=1; error('''sigma_ub'' must be a row vector'); end
        if M~=m; error('''sigma'' and''sigma_ub'' must have the same number of elements'); end
    end
end

%%

%%%%%%%%%%%%%%%%%%%
% Produce plot
%%%%%%%%%%%%%%%%%%%

[A,B]=size(col);
L=ceil(M/A);
clrs=repmat(col,L,1);

labelinput_new=cell(1,M);

if sorting
    [mi,Sidx]=sort(mi,'descend') ;
    mi_ub = mi_ub(Sidx) ;
    mi_lb = mi_lb(Sidx) ;
    sigma = sigma(Sidx) ;
    sigma_ub = sigma_ub(Sidx) ;
    sigma_lb = sigma_lb(Sidx) ;
    for i=1:M; labelinput_new{i} = labelinput{Sidx(i)} ;end;
end

if nb_legend<M
    labelinput_new=labelinput_new(1:nb_legend);
    labelinput_new{end}=[labelinput_new{end},'...'];
end


% First plot EEs mean & std as 1: square, 2: circles, 3: triangle, 4: upside-down triangle:

for i=1:M        
    switch k    % change markers 
        case 1
            plot(mi(i),sigma(i),'sk','MarkerFaceColor',clrs(i,:),'MarkerSize',ms,'MarkerEdgeColor','k')
        case 2
         	plot(mi(i),sigma(i),'ok','MarkerFaceColor',clrs(i,:),'MarkerSize',ms,'MarkerEdgeColor','k')
        case 3
            plot(mi(i),sigma(i),'^k','MarkerFaceColor',clrs(i,:),'MarkerSize',ms,'MarkerEdgeColor','k')
        case 4
            plot(mi(i),sigma(i),'vk','MarkerFaceColor',clrs(i,:),'MarkerSize',ms,'MarkerEdgeColor','k')
    end
end

%plot first the larger confidence areas
size_bounds=mi_ub-mi_lb; 
[tmp,idx]=sort(size_bounds,'descend');

for i=1:M % add rectangular shade:    
    h = fill([mi_lb(idx(i)),mi_lb(idx(i)),mi_ub(idx(i)),mi_ub(idx(i))],[sigma_lb(idx(i)),sigma_ub(idx(i)),sigma_ub(idx(i)),sigma_lb(idx(i))],clrs(idx(i),:));
end

% Plot again the markers (in case some have been overriden by the rectangles
% representing confidence bounds)
for i=1:M        
    switch k    % change markers 
        case 1
            plot(mi(i),sigma(i),'sk','MarkerFaceColor',clrs(i,:),'MarkerSize',ms,'MarkerEdgeColor','k')
        case 2
         	plot(mi(i),sigma(i),'ok','MarkerFaceColor',clrs(i,:),'MarkerSize',ms,'MarkerEdgeColor','k')
        case 3
            plot(mi(i),sigma(i),'^k','MarkerFaceColor',clrs(i,:),'MarkerSize',ms,'MarkerEdgeColor','k')
        case 4
            plot(mi(i),sigma(i),'vk','MarkerFaceColor',clrs(i,:),'MarkerSize',ms,'MarkerEdgeColor','k')
    end
end

end
% Create legend:
legend(labelinput_new, 'Location','NorthWest')
% xlim([0 100])
% ylim([0 70])   
xlim([0 100])
ylim([0 12])
xlabel('Mean of EEs','FontSize',fs,'FontName',fn)
ylabel('Standard deviation of EEs','FontSize',fs,'FontName',fn)
grid on
set(gca,'FontSize',fs,'FontName',fn)
box on
end


