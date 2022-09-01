%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
clf
% main program for the agent-based model of zebrafish pattern formation
% load pameters
zebra_pars_abm15

% set random seed
% rng(2)
    
% location of horizontal myoseptum, the yellow strip
hmwidth = 100; % width of hm
iymin = (Dwidth-hmwidth)/2;
iymax = iymin + hmwidth;

hmloc = [0+rx Dlength-rx iymin iymax];

%%% initialization

% position of each cell at the initial step, 2022 a single cell?
xpm_1 = (rall(1)/2:(rall(1)*2+40):Dlength)'; 
xpm_2 = (rall(1)/2:(rall(1)*2+40):Dlength)';
xpm_3 = (rall(1)/2:(rall(1)*2+40):Dlength)';
xpm_4 = (rall(1)/2:(rall(1)*2+40):Dlength)';
pm = [xpm_1, ones(size(xpm_1))*(rall(1)*2+10); xpm_2, ones(size(xpm_2))*(Dwidth/2-rall(2)*2-20); xpm_3, ones(size(xpm_3))*(Dwidth/2+rall(2)*2+20); xpm_4, ones(size(xpm_4))*(Dwidth-rall(1)*2);];
%px = hmcellpos(hmloc); 
xpx = (rall(2)/2:(rall(2)*2+5):Dlength)';
px = [xpx, ones(size(xpx))*Dwidth/2];

% we put a strip of dense iphos in the horizontal myoseptum region
% maxinit = 3000;%floor(hmwidth/(2*rx)*Dlength/(2*rx)); % maximum number of xphores initialized in the hm
% initnid = randi([10 maxinit],1); % the two numbers are the minimum and maximum number of cells on the horizontal myoseptum
% for i=1:initnid
%     newpos = hmcellpos(hmloc);
%     if nooverlap_test(newpos,px,rx,rx*ones(size(px,1),1))
%         px = [px;newpos]; % n by 2
%     end
% end



% colr = [0 0 0];
% f1 = axes;
% plotcells(pm,rm,colr,0.3,f1);
% hold on
% colr = [1 0.65 0];
% plotcells(px,rx,colr,0.3,f1);
% axis equal
% axis([0 Dlength 0 Dwidth])
%%
clc
%%%%%%%% Store all of the information for ploting purpose later
S(1).domsize = [Dlength, Dwidth];
S(1).pos = {pm,px};
domx = Dlength;
domy = Dwidth;

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
            nposall{indi} = fnmove(domx,domy,nposall{1},nposall{2},rall,nall,indi); % need to be modified to take "nposall{1},nposall{2}" automatically
        end
    end
    

%     % update cell birth/division
    nbpos = cell(1,ntype);
    for indi = 1:ntype
        nbpos{indi} = fnbirth(domx,domy,pm,px,rall,indi,par_birth,dpar); % output are lists of positions of newly created cells
    end
% 
%     
%     % update cell death 
%     % return the indexes of cell that undergo cell death
    indxdth = cell(1,ntype);
    for indi = 1:ntype
        indxdth{indi} = fndeath(pm,px,rall,celltype{indi},gammas,dpar);
    end
    
   
%     
%     % after all calculations, update the status of each type of cells
    indxcov_new = cell(1,ntype);
    cov_pos = cell(1,ntype);
    cov_to = [1 3 2 5 4];
%     for indi = 2:ntype
%         indxcov_new{indi} = setdiff(indxcov{indi},indxdth{indi}); % to account for cell conversion, remove cells undergo death
%         nposall{cov_to(indi)} = [nposall{cov_to(indi)}; nposall{indi}(indxcov_new{indi},:)];
%         indxdth{indi} = union(indxdth{indi},indxcov{indi}); % remove all cells that undergo conversion and cell death
%     end

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
    
    % NOTE HERE WE ASSUME THE RADIUS OF EACH IS NOT GROWING AS DOMAIN
    % GROWS!!!
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
        axis([0 5000 0 4000]);
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






