function pos = hmcellpos(hmloc)
xmin = hmloc(1);
xmax = hmloc(2);
ymin = hmloc(3);
ymax = hmloc(4);
pos = [0 0];
pos(1) = (xmax-xmin)*rand(1,1) + xmin;
pos(2) = (ymax-ymin)*rand(1,1) + ymin;
