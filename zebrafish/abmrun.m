function S = abmrun(totd,varargin)

run_start = tic;

%% 1. Setting up parameters.
default_param = abmset('ntype',2, ...
                       'domx',2000, ...
                       'domy',1000, ...
                       'domxt',130, ...
                       'domyt',130, ...
                       'rm',20, ...
                       'rx',20, ...
                       'gamma_loc', 75, ...
                       'gamma_long_inner', 318, ...
                       'gamma_long_outer', 343, ...
                       'mu', 1, ...
                       'nu', 1, ...
                       'psi', 1.2, ...
                       'p_d', 1/30, ...
                       'alpha', 1, ...
                       'beta', 3.5, ...
                       'eta', 6, ...
                       'd_crowd', 82, ...
                       'phi_1', 1.3, ...
                       'phi_2', 1.2, ...
                       'kappa', 10, ...
                       'd_rand', 100, ...
                       'a', 3000, ...
                       'b', 300, ...
                       'hmwidth', 100, ...
                       'plot_init', false);

param = default_param;

if nargin == 0 % 'totd' not even given
    % Number of days in reality as well as number of timesteps in the
    % Forward Euler stepping scheme with dt = 1.
    totd = 20; % At least 20, due to zebrafish biology    
elseif nargin == 1 % assume only 'totd' given
    assert(isnumeric(totd) & totd < 1000)
else % assume the third argument onwards modify the default parameters
    param = abmset(default_param, varargin{:});
end

ntype = param.ntype;
domx = param.domx; 
domy = param.domy;

% Default 130 from (16-domx/1000)*10^3/totd, using domx = 3000, totd = 100
% Jia's note: 7 is the initial size, 16 is the destinate size
domxt = param.domxt; 
% Jia's note: (3-domx/1000)*10^3/totd. Inconsistent.
% 1 is the initial size, 16 is the destinate size
domyt = 130; 

% 'r' means radius, but in the ABM, cells are treated as points.
% the default cell sizes were only used initially to select cell positions
rm = param.rm; % radius of melanophores
rx = param.rx; % radius of xanthophores

rall = [rm, rx]; % for iteration.

% effective neighborhood of a single cell where it exerts influence
gamma_loc = param.gamma_loc;

% "long" for "long-distance" interaction through cytoneme, lower and upper
% bounds of interaction given by inner and outer
gamma_long_inner = param.gamma_long_inner;
gamma_long_outer = param.gamma_long_outer; % gamma_long_inner + 25;
gammas = [gamma_loc, gamma_long_inner,gamma_long_outer];

% life-death dynamics
mu = param.mu;
nu = param.nu;
psi = param.psi;
p_d = param.p_d; % baseline death rate regardless of cell type

dpar = [mu,nu,psi,p_d]; % death parameters

% birth related parameter
alpha = param.alpha; 
beta = param.beta;% concerns ratios of melano/xanthophores, see supp. of Volkening 2015
eta = param.eta; % for melanophores
d_crowd = param.d_crowd; % minimum empty space for cells to be born (in micrometers)
phi_1 = param.phi_1;
phi_2 = param.phi_2;
kappa = param.kappa; % for xanthophores
d_rand = param.d_rand;

% defines a linear daily growth of cells born 'ndiff' in fnbirth.m
a = param.a;
b = param.b;

par_birth = [alpha, beta, eta, d_crowd, phi_1, phi_2, kappa, d_rand, a, b];

%% 2. Computing initial cell distribution.

hmwidth = param.hmwidth; % width of horizontal myoseptum, the init. yellow strip

% See initial plot for meaning.
iymin = (domx-hmwidth)/2;
iymax = iymin + hmwidth;
hmloc = [0+rx domx-rx iymin iymax];

% position of each cell at the initial step
xpm_1 = (rm/2:(rm*2+20):domx)'; % 34 cells
xpm_2 = (rm/2:(rm*2+20):domx)';
xpm_3 = (rm/2:(rm*2+20):domx)';
xpm_4 = (rm/2:(rm*2+20):domx)';
pm = [xpm_1, ones(size(xpm_1))*(rm*2+1); xpm_2, ones(size(xpm_2))*(domy/2-rx*2-20); xpm_3, ones(size(xpm_3))*(domy/2+rx*2+20); xpm_4, ones(size(xpm_4))*(domy-rm*2);];

xpx = (rx/2:(rx*2)-0.5:domx)'; % 51 cells
px = [xpx, ones(size(xpx))*domy/2];
%size(xpx) %for debug

% Plotting initial pattern?
if param.plot_init
    colr = [0 0 0];
    f1 = axes;
    plotcells(pm,rm,colr,0.3,f1);
    hold on
    colr = [1 0.65 0];
    plotcells(px,rx,colr,0.3,f1);
    axis equal
    axis([0 domx 0 domy])
end

% Store all of the information for ploting purpose later
S(1).domsize = [domx, domy];
S(1).pos = {pm,px}; % pm = positions of melano, px = positions of xantho

%% 3. Main loop
for indt = 1:totd % it was totd+1 for unknown reasons.
%     % spot check whether argument changes took effect
%      fprintf('DEBUG: nu = %g\n',nu) 
    
    fprintf('Day %d begins... \n', indt)
    fprintf('Domain sizes: (%d, %d)\n', domx, domy)
    fprintf('Number of cells: (%d, %d)\n', length(pm), length(px))
    tstart_main = tic;
    
    % find out how many life cells of each type we have
    nm = size(pm,1);
    nx = size(px,1);
    nall = [nm, nx];
    
    % nposall, nbpos, indxdth stores intermediate values to prepare for real update
    nposall = {pm,px};
    
    fprintf("Computing cell movement...\n")
    s1 = tic;
    for indi = 1:ntype
        if ~isempty(nposall{indi}) % if there are nonzero number of type indi cells
            nposall{indi} = fnmove(domx,domy,nposall{1},nposall{2},rall,nall,indi); 
        end
    end
    toc(s1);
   
    fprintf("Computing newborn cell birth/division...\n")
    s2 = tic;
    nbpos = cell(1,ntype);
    for indi = 1:ntype
        nbpos{indi} = fnbirth(domx,domy,pm,px,rall,indi,par_birth,gammas,indt);
    end
    toc(s2);
    
    fprintf("Computing cell death ...\n")
    s3 = tic;
    indxdth = cell(1,ntype);
    for indi = 1:ntype
            % return the indexes of cell that undergo cell death
        indxdth{indi} = fndeath(pm,px,rall,indi,gammas,dpar); % return indices
    end
    toc(s3);

    
%     % after all calculations, update the status of each type of cells
%     indxcov_new = cell(1,ntype);
%     cov_pos = cell(1,ntype);
%     cov_to = [1 3 2 5 4];

    

    fprintf("Real update of birth, death and movement data...\n") 
    tic;
    for indi = 1:ntype
        tempposi = nposall{indi};
        tempposi = [tempposi; nbpos{indi}]; % add location of new born cells
        tempposi(indxdth{indi},:) = []; % eliminate dead cells
        nposall{indi} = tempposi;
        if isempty(nposall{indi})
            nposall{indi} = [];
        end
    end
    
    % assign new positions to the position vector
    pm = nposall{1};
    px = nposall{2};
    toc

    fprintf("Updating locations due to domain growth...\n")
    tic;
    domx = domx + domxt; % add 130um to width and height everyday
    domy = domy + domyt; % % add 130um to width and height everyday
    strx = domx/(domx-domxt); % stretch rate, minus the growth of one day, compared with the previous day
    stry = domy/(domy-domyt);
    
    % stretching effects of cell locations due to domain growth, does not
    % change number of cells
    pm = pm .* repmat([strx, stry],size(pm,1),size(pm,2)/2); % the the 3rd argument of repmat is artificial to ensure this works for empty pm
    px = px .* repmat([strx, stry],size(px,1),size(px,2)/2);
    toc;
    
    fprintf("Storing cell locatinos...\n")
    tic;
    S(indt+1).domsize = [domx domy]; % index 1 is used for initial condition 
    S(indt+1).pos = {pm,px};
    toc;

    fprintf("This time step took %.3g seconds.\n\n",toc(tstart_main))

end % simulation has finished

fprintf("The whole run for %i timesteps took %.3g seconds.\n",totd, toc(run_start))
fprintf("User inputs:\n")
disp(varargin)

% S contains all simulation data, totd, ntype can be inferred.

% Jia's notes:
% radius growth should be each to implement: just add the change in radius
% with time after every time step

% variation in the birth or death rate is easy to implement: we just need
% to change the limits for each birth and death step
end