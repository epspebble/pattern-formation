%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
clf
% main program for the agent-based model of zebrafish pattern formation
% L = 318;
% omega = 25;
% l = 75;
% 
% xi = 1.2;
% pdeath = 0.0333;
% Ndiff = 1; %
% 
% dcrowd = 82;
% phi = 1.3;
% psi = 1.2;
% 
% drand = 100;
% pm = 0.03;
% px = 0.005;


% other parameters
ntype = 2;


domx = 2000; % initialize domain
domy = 1000; %initialize domain

totd = 20; % today number of days, tried 20, but on Volkening 2015, Fig 14, up to 150 days

domxt = 130; %(16-domx/1000)*10^3/totd; % 7 is the initial size, 16 is the destinate size
domyt = 130;% (3-domy/1000)*10^3/totd; % 1 is the initial size, 16 is the destinate size

rm = 20; % mphos
rx =  20; % dense xphos
rall = [rm, rx];



gamma_loc = 75;
gamma_long_inner = 318;
gamma_long_outer = gamma_long_inner + 25;
gammas = [gamma_loc, gamma_long_inner,gamma_long_outer];

% par for death
mu = 1;
nu = 1;
psi = 1.2;
p_d = 0.0333;
dpar = [mu,nu,psi,p_d];

% par for birth/differentiation 
alpha = 1;
beta = 3.5;
eta = 6;

d_crowd = 82;

phi_1 = 1.3;
phi_2 = 1.2;
kappa = 10;

d_rand = 100;
par_birth = [alpha, beta, eta, d_crowd, phi_1, phi_2, kappa, d_rand];






% set random seed
% rng(2)
    
% location of horizontal myoseptum, the yellow strip
hmwidth = 100; % width of hm
iymin = (domx-hmwidth)/2;
iymax = iymin + hmwidth;

hmloc = [0+rx domx-rx iymin iymax];

%%% initialization

% position of each cell at the initial step
xpm_1 = (rall(1)/2:(rall(1)*2+20):domx)'; % 34 cells
xpm_2 = (rall(1)/2:(rall(1)*2+20):domx)';
xpm_3 = (rall(1)/2:(rall(1)*2+20):domx)';
xpm_4 = (rall(1)/2:(rall(1)*2+20):domx)';
pm = [xpm_1, ones(size(xpm_1))*(rall(1)*2+1); xpm_2, ones(size(xpm_2))*(domy/2-rall(2)*2-20); xpm_3, ones(size(xpm_3))*(domy/2+rall(2)*2+20); xpm_4, ones(size(xpm_4))*(domy-rall(1)*2);];

xpx = (rall(2)/2:(rall(2)*2)-0.5:domx)'; % 51 cells
px = [xpx, ones(size(xpx))*domy/2];
% size(xpx)

colr = [0 0 0];
f1 = axes;
plotcells(pm,rm,colr,0.3,f1);
hold on
colr = [1 0.65 0];
plotcells(px,rx,colr,0.3,f1);
axis equal
axis([0 domx 0 domy])
%%
clc
%%%%%%%% Store all of the information for ploting purpose later
S(1).domsize = [domx, domy];
S(1).pos = {pm,px};
domx = domx;
domy = domy;

celltype = {'m','x'};
for indt = 1:totd+1
    
    % initialize number of each type of cells in the domain
    nm = size(pm,1);
    nx = size(px,1);
    nall = [nm, nx];


    nposall = {pm,px};

    % update cell movement
    for indi = 1:ntype
        if ~isempty(nposall{indi}) % if there are nonzero number of type indi cells
            nposall{indi} = fnmove(domx,domy,nposall{1},nposall{2},rall,nall,indi); 
        end
    end
    

%     % update cell birth/division
    nbpos = cell(1,ntype);
    for indi = 1:ntype
        nbpos{indi} = fnbirth(domx,domy,pm,px,rall,indi,par_birth,gammas,indt);
    end
    
%     % update cell death 
%     % return the indexes of cell that undergo cell death
    indxdth = cell(1,ntype);
    for indi = 1:ntype
        indxdth{indi} = fndeath(pm,px,rall,indi,gammas,dpar);
    end
    
   
%     
%   % after all calculations, update the status of each type of cells
    indxcov_new = cell(1,ntype);
    cov_pos = cell(1,ntype);
    cov_to = [1 3 2 5 4];

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

    
    %%%%%%%%%%%%% update domain growth
    domx = domx + domxt; % add 130um to width and height everyday
    domy = domy + domyt; % % add 130um to width and height everyday
    strx = domx/(domx-domxt); % stretch rate, minus the growth of one day, compared with the previous day
    stry = domy/(domy-domyt);
    
    % stretching effects of cell locations due to domain growth 
    pm = pm .* repmat([strx, stry],size(pm,1),size(pm,2)/2); % the the 3rd argument of repmat is artificial to ensure this works for empty pm
    px = px .* repmat([strx, stry],size(px,1),size(px,2)/2);
    
    
    % for storage only 
    S(indt+1).domsize = [domx domy]; % index 1 is used for initial condition 
    S(indt+1).pos = {pm,px};
end

%% plot the cell positions as a function of time 
% clc
clf
plt = totd;
%colorlist = {[0 0 0],[1 0.65 0],[1 1 0],[0.95 0.97 1],[0.96 0.96 0.96]}; % black, orange, yellow, light blue, white
colorlist = {[0 0 0],[1 0.65 0],[1 1 0],[0.65 0.67 1],[0.67 0.78 0.94]}; % black, orange, yellow, light blue, white

translist = [0.7 0.3 0.3 0.3 0.3]; % list of transparency value
for indt = 1:plt
    %hold on
%     f1 = axes('Position',[0.1 0.1 0.77 0.81]);
%     hold(f1,'on');
    hold on
    for indi = 1:ntype
        plotcells(S(indt).pos{indi},rall(indi),colorlist{indi},translist(indi));
        title(['Time at ', num2str(indt),' days']);
        plot([0 S(indt).domsize(1)],[0 0],'b');
        plot([0 S(indt).domsize(1)],[S(indt).domsize(2) S(indt).domsize(2)],'b');
        axis equal
        %axis([0 5000 0 4000]);
    end
    pause(0.2)
    hold off
    % hold(f1,'off');
    if indt ~= plt
        clf;
    end
        
end
% radius grwoth should be each to implement: just add the change in radius
% with time after every time step

% variation in the birth or death rate is easy to implement: we just need
% to change the limits for each birth and death step






