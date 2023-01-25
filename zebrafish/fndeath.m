function indxtypei = fndeath(pm,px,rall,indc,gammas,dpar)
% input: 2 lists of cell locations, radius of each type of cell, flag
% indicates which cell death we are looking at

% return a list of cell index to be eliminated 
gamma_loc = gammas(1);
gamma_long_inner = gammas(2);
gamma_long_outer = gammas(3);


mu = dpar(1);
nu = dpar(2);
psi = dpar(3);
p_d = dpar(4);

indxtypei = [];

if indc == 1

    for indi=1:size(pm,1)
        % effect of xphores on mphore death, short distance
        dsxm = dists(pm(indi,:), px);
        temppm = pm;
        temppm(indi,:) = [];
        dsmm = dists(pm(indi,:), temppm);
        num_mnb = sum(dsmm<gamma_loc);
        num_xnb = sum(dsxm<gamma_loc);
        eff_sxm = num_xnb>mu*num_mnb;
        
        num_man = sum(dsmm>gamma_long_inner & dsmm<gamma_long_outer);
        num_xan = sum(dsxm>gamma_long_inner & dsxm<gamma_long_outer);
        eff_lxm = num_man > psi*num_xan;
        rdd = rand();
        if eff_sxm || (eff_lxm && rdd<p_d)
            indxtypei = [indxtypei,indi];
        end
    end
elseif indc == 2

    for indi=1:size(px,1)
        % effect of mphores on xphore death
        dsmx = dists(px(indi,:), pm);
        temppx = px;
        temppx(indi,:) = [];
        dsxx = dists(px(indi,:), temppx);
        num_mnb = sum(dsmx<gamma_loc);
        num_xnb = sum(dsxx<gamma_loc);
        eff_mx = num_mnb > nu*num_xnb;

        if eff_mx
            indxtypei = [indxtypei,indi];
        end
    end
end
