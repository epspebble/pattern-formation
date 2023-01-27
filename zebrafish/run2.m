function run2(n)
% n to be provided by SLURM_ARRAY_TASK_ID

% let's do n = 11 x 21 = 231
% a = 2000, 2500, 3000, ..., 7000 (11 numbers)
% b = 300, 310, 320, ..., 500 (21 numbers)

a = 2000+mod(n-1,11)*500;
b = 300+floor((n-1)/11)*10;
S = abmrun(10,"a",a,"b",b);

fn = sprintf("run2_a%i_b%i",a,b);
fprintf("Saving data to %s.mat ...\n",fn);
save(fn,"S")
