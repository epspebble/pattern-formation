function plotcells(centers,radius,colr,falpha,f1)
%f1=figure(1);
% viscircles(centers,ones(1,size(centers,1))*radius)

%circles
k_circle_resolution=20;
indx = 1:k_circle_resolution+1;
x_unit_circle = cos(indx/k_circle_resolution*2*pi);
y_unit_circle = sin(indx/k_circle_resolution*2*pi);


hold on
for indi = 1:size(centers,1)
%    th = 0:pi/50:2*pi;
%     xunit = radius * cos(th) + centers(indi,1);
%     yunit = radius * sin(th) + centers(indi,2);
%     
%     fill(f1,xunit, yunit, colr,'EdgeColor',colr,'FaceALpha',falpha);
    %set(f1,'xticklabel',[],'yticklabel',[])
    
    x_cell = centers(indi,1) + radius*x_unit_circle;
    y_cell = centers(indi,2) + radius*y_unit_circle;
    plot(x_cell,y_cell,'Color',colr);
end