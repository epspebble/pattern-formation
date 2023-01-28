function nposlist = fnbirth(domx,domy,pm,px,rall,indx,par_birth,gammas,indt)


alpha = par_birth(1);
beta = par_birth(2);
eta = par_birth(3);
d_crowd = par_birth(4);
phi_1 = par_birth(5);
phi_2 = par_birth(6);
kappa = par_birth(7);
d_rand = par_birth(8);

gamma_loc = gammas(1);
gamma_long_inner = gammas(2);
gamma_long_outer = gammas(3);

nposlist = [];

if indx == 1
    brate = 0.03;
elseif indx == 2
    brate = 0.005;
end

ndiff = 2000+indt*200; % not specified in the paper 
for ind = 1:ndiff
    xypos = zeros(1,2); % tentative position of new born cells
    xypos(1,1) = (domx-2*rall(indx))*rand(1,1) + rall(indx); % domx, length of the domain in the x direction
    xypos(1,2) = (domy-2*rall(indx))*rand(1,1) + rall(indx); % domy, length of the domain in the y direction    birth_m = dists(xypos, pm);
    if indx == 1
        dsmm = dists(xypos, pm);
        dsxm = dists(xypos, px);
        num_loc_mnb = sum(dsmm<gamma_loc);
        num_loc_xnb = sum(dsxm<gamma_loc); 
        num_podia_mnb = sum(dsmm>gamma_long_inner & dsmm<gamma_long_outer);
        num_podia_xnb = sum(dsxm>gamma_long_inner & dsxm<gamma_long_outer);
        cond1 = num_loc_mnb>alpha*num_loc_xnb;    
        cond2 = num_podia_xnb>beta*num_podia_mnb; % wrong expression in the supplement
        num_loc_mnb = sum(dsmm<d_crowd);
        num_loc_xnb = sum(dsxm<d_crowd);
        cond3 = num_loc_xnb + num_loc_mnb < eta;
        r = rand();
        cond4 = r<=brate;
    elseif indx == 2
        dsmm = dists(xypos, pm);
        dsxm = dists(xypos, px);
        num_loc_mnb = sum(dsmm<gamma_loc);
        num_loc_xnb = sum(dsxm<gamma_loc);  
        num_podia_mnb = sum(dsmm>gamma_long_inner & dsmm<gamma_long_outer);
        num_podia_xnb = sum(dsxm>gamma_long_inner & dsxm<gamma_long_outer);
        cond1 = num_loc_xnb>phi_1*num_loc_mnb;    
        cond2 = num_podia_mnb>phi_2*num_podia_xnb; % might be wrong here according to the description in the paper
        num_loc_mnb = sum(dsmm<d_crowd);
        num_loc_xnb = sum(dsxm<d_crowd);
        cond3 = num_loc_xnb + num_loc_mnb < kappa;
        r = rand();
        cond4 = r<=brate;
    end
    if (cond1 && cond2 && cond3 )|| (cond4  && sum(dsmm<d_rand)==0 && sum(dsxm<d_rand) == 0)
       nposlist = [nposlist;xypos];
    end
end
    
    
% nnew = 100;
% for i = 1:nnew
%     r = rand();
%     xypos = zeros(1,2); % tentative position of new born cells
%     xypos(1,1) = (domx-2*rall(indx))*rand(1,1) + rall(indx); % domx, length of the domain in the x direction
%     xypos(1,2) = (domy-2*rall(indx))*rand(1,1) + rall(indx); % domy, length of the domain in the y direction    birth_m = dists(xypos, pm);
%     dsmm = dists(xypos, pm);
%     dsxm = dists(xypos, px);
%     if indx == 1 && sum(dsmm<d_rand)==0 && sum(dsxm<d_rand) == 0% type m 
%        if r<=brate
%             % new cells because of cell division
%            nposlist = [nposlist;xypos];
%        end
%     elseif indx == 2 && sum(dsmm<d_rand)==0 && sum(dsxm<d_rand) == 0 % type x 
%        if r<=brate
%           nposlist = [nposlist;xypos];
%        end
%     end
% end

