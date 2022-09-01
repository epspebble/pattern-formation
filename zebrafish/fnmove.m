function newpos = fnmove(domx,domy,pm,px,rall,nall,indx)
% pm: position of all m cells
% px: position of all x cells
% indx: the indx-th cell of the current cell type
% 
Rxm = 10000;
rxm = 15;
Axm = 10;
axm = 10;


Rmm = 5000;
rmm = 13;
Amm = 10;
amm = 10;


Rmx = 3000;
rmx = 20;
Amx = 1500;
amx = 3;


Rxx = 2000;
rxx = 10;
Axx = 10;
axx = 10;



dt = 1;

% length of nalls is the total number of cells
% nalls contains cell radius for all cells
%nalls 


allposcell = {pm,px};
curtype = allposcell{indx}; % current type, n by 2
newpos = zeros(size(curtype));

if indx == 1 
    allpos = [pm];
    % nalls = repelem(rall(1),nall(1)).'; % every elements in vector rall get repeated nall(1) times
elseif indx == 2 
    allpos = [px];
    % nalls = repelem(rall(2),nall(2).');
end


% location of each type of cells in each layer
locs = {1:nall(1),nall(1)+1:sum(nall(1:2))}; 

% update cell movement in arbitary order
randorder = 1:size(curtype,1);
randorder = randorder(randperm(length(randorder)));




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
        Fmm = Qkl(cellipos, noipos, Rmm, rmm, Amm, amm); % short range?
        Fxm = Qkl(cellipos, px, Rxm, rxm, Axm, axm);
        diffpos = -sum(Fmm,1) - sum(Fxm,1);
        
    elseif indx == 2
        Fxx   = Qkl(cellipos, noipos, Rxx, rxx, Axx, axx);
        Fmx   = Qkl(cellipos, pm, Rmx, rmx, Amx, amx);
        diffpos = -sum(Fxx,1) - sum(Fmx,1);             
    end
    
    % first update the location of the current cell
    ncelli = cellipos + diffpos*dt; % dt should be 1
  
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
