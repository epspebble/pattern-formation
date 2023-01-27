% Test run on a remote cluster for a single job.
runst = tic;

S = abmrun(50,"a",3500);
fprintf("\n\n\nFinished simulation.\n")
toc(runst)

fn = "run2_a3500";
fprintf("Saving data to %s".fn)
save(fn,"S")
