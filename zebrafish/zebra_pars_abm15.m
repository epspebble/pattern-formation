% parameters for the zebrafish model

% 2 types of cells
ntype = 2;


% radii of the 5 types of cells
rm = 20; % mphos
rx =  20; % dense xphos


% 
rall = [rm, rx];


% characteristic distances
dneb = 10; % represent the neiborhood of a cell 
dxx = 36; % mu m,contact interactions 
dxm = 50; % short-range  
dmm = 82; % long-range

ds = dxx;



% time duration of consideration

% growth rate of domain
% domxt = (16-7)*10^3/totd; % 2mm is the initial size, 16 is the destinate size
% domyt = (3-1)*10^3/totd; % 1mm is the initial size, 3 is the destinate size

% initial size of the domain of consideration
Dlength = 2000; % 2mm
Dwidth = 1000;  % 1mm
domx = Dlength; % initialize domain
domy = Dwidth; %initialize domain
totd = 20;
domxt = 130; %(16-Dlength/1000)*10^3/totd; % 7 is the initial size, 16 is the destinate size
domyt = 130;% (3-Dwidth/1000)*10^3/totd; % 1 is the initial size, 16 is the destinate size


nm = 200; % number of random locations selected each dy for possible mphore birth (induced by other cells)

% number of locations for possible cell birth
nnew = 100;



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


