function newpos = fnmove(domx,domy,pm,px,rall,nall,indx)
% pm: position of all m cells
% px: position of all x cells
% indx: the indx-th cell of the current cell type


% Jia inferred from Supp. Fig. 2 by manual matching. There is no biological
% basis for such repulsive force. It is due to an empirical observation
% that these cells move away from each other.
Rxm = 11000;
rxm = 20;
Axm = 0;
axm = 1;

Rmm = 5000;
rmm = 20;
Amm = 0;
amm = 1;

Rmx = 6000;
Amx = Rmx-1600;
amx = 9;
rmx = (Rmx*amx)/(Amx-200*amx);

Rxx = 2000;
rxx = 10;
Axx = 0;
axx = 1;

dt = 1;
fprintf("   Working on cell type %i.\n", indx)
fprintf("   Step 1 - randomizing... ")
tic;

allposcell = {pm,px};
curtype = allposcell{indx}; % current type, n by 2
newpos = zeros(size(curtype));

% update cell movement in arbitary order
randorder = 1:size(curtype,1);
randorder = randorder(randperm(length(randorder)));

fprintf("   Took %.2g seconds.\n" ,toc);


fprintf("   Step 2 - computing force function and move... ")
tic;
for ind = 1:size(curtype,1)
    
    indi = randorder(ind);
    % The location of the indi-th cell of type indi
    cellipos = curtype(indi,:); % x and y position, 1 by 2
    
    % noipos contains all cells of the current type but not the current
    % cell; no ith position
    noipos = curtype;
    noipos(indi,:) = []; % n-1 by 2
    
    if indx == 1
        % movement of m cells
        % When distance = Lmm, force = 0; if dist < Lmm, force >0
        Fmm = Qkl(cellipos, noipos, Rmm, rmm, Amm, amm); 
        Fxm = Qkl(cellipos, px, Rxm, rxm, Axm, axm);
%         diffpos = -sum(Fmm,1) - sum(Fxm,1);
        diffpos = -sum(Fmm) - sum(Fxm);
        
    elseif indx == 2
        % movement of x cells...
        Fxx   = Qkl(cellipos, noipos, Rxx, rxx, Axx, axx);
        Fmx   = Qkl(cellipos, pm, Rmx, rmx, Amx, amx);
%        diffpos = -sum(Fxx,1) - sum(Fmx,1);             
        diffpos = -sum(Fxx) - sum(Fmx);             
    end
    
    % first update the location of the current cell
    ncelli = cellipos + diffpos*dt; % forward Euler, dt=1 is fixed.
  
    % if the new location overlaps with other cell, consider shorten the
    % distance of movement
    if ncelli(1)>domx
        ncelli(1) = domx-ncelli(1);
    elseif ncelli(1)<0
        ncelli(1) = -ncelli(1);
    end
    
    if ncelli(2)>domy
        ncelli(2) = domy-ncelli(2);
    elseif ncelli(2)<0
        ncelli(2) = -ncelli(2);
    end
    newpos(indi,:) = ncelli;
    
end

fprintf("   Took %.2g seconds.\n", toc);
