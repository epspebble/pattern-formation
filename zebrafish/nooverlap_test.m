function flag = nooverlap_test(newx,plist,rad1,rad2)
% rad1: radius of newx
% rad2: could be a list, radius of plist
distances = dists(newx,plist);
if any(distances<(rad1+rad2)) % need to multiply by two to account for the distance between two cells
    flag = false;
else
    flag = true;
end