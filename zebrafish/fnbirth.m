function nposlist = fnbirth(domx,domy,pm,px,rall,indx,par_birth,gammas)


alpha = par_birth(1);
beta = par_birth(2);
eta = par_birth(3);
d_crowd = par_birth(4);
phi_1 = par_birth(5);
phi_2 = par_birth(6);
kappa = par_birth(7);
d_rand = par_birth(8);

gamma_loc = gammas(1);

nposlist = [];

if indx == 1
    brate = 0.03;
elseif indx == 2
    brate = 0.005;
end

ndiff = ceil((domx*domy)/1000); % not specified 
for ind = 1:ndiff
    xypos = zeros(1,2); % tentative position of new born cells
    xypos(1,1) = (domx-2*rall(indx))*rand(1,1) + rall(indx); % domx, length of the domain in the x direction
    xypos(1,2) = (domy-2*rall(indx))*rand(1,1) + rall(indx); % domy, length of the domain in the y direction    birth_m = dists(xypos, pm);
    if indx == 1
        birth_m = dists(xypos, pm);
        birth_x = dists(xypos, px);
        num_mnb = sum(birth_m<gamma_loc);
        num_xnb = sum(birth_x<gamma_loc);
        cond1 = num_mnb>alpha*num_xnb;    
        cond2 = num_xnb>beta*num_mnb;
        num_mnb = sum(birth_m<d_crowd);
        num_xnb = sum(birth_x<d_crowd);
        cond3 = num_xnb + num_mnb < eta;
    elseif indx == 2
        birth_m = dists(xypos, pm);
        birth_x = dists(xypos, px);
        num_mnb = sum(birth_m<gamma_loc);
        num_xnb = sum(birth_x<gamma_loc);
        cond1 = num_xnb>phi_1*num_mnb;    
        cond2 = num_mnb>phi_2*num_xnb;
        num_mnb = sum(birth_m<d_crowd);
        num_xnb = sum(birth_x<d_crowd);
        cond3 = num_xnb + num_mnb < kappa;
    end
    if cond1 && cond2 && cond3
       nposlist = [nposlist;xypos];
    end
    
    
    nnew = 50;
    for i = 1:nnew
        r = rand();
        xypos = zeros(1,2); % tentative position of new born cells
        xypos(1,1) = (domx-2*rall(indx))*rand(1,1) + rall(indx); % domx, length of the domain in the x direction
        xypos(1,2) = (domy-2*rall(indx))*rand(1,1) + rall(indx); % domy, length of the domain in the y direction    birth_m = dists(xypos, pm);

        if indx == 1 && sum(birth_m<d_rand)==0 && sum(birth_x<d_rand) == 0% type m 
           if r<=brate
                % new cells because of cell division
               nposlist = [nposlist;xypos];
           end
        elseif indx == 2 && sum(birth_m<d_rand)==0 && sum(birth_x<d_rand) == 0 % type x 
           if r<=brate
              nposlist = [nposlist;xypos];
           end
        end
    end
end
