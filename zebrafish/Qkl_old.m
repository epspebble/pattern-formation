function force = Qkl(pt, plist, Rkl, rkl, Akl, akl)
% input the point to be investigated: pt
% and the list of points possibly in effect: plist
% rad1: radius of the point pt
% rad2: common radius of the list of points plist
% Lkl: characteristic length between the two types of cells
% epskl: force strength of the two types of cells

% Q-function 
force = zeros(1,2);
if ~isempty(plist)
    dkl = plist-pt; % n by 2 matrix
    x = dkl(:,1);
    y = dkl(:,2);
    Qx = -Rkl*exp(-sqrt(x.^2+y.^2)/rkl)*2.*x/(2*rkl*sqrt(x.^2+y.^2))+Akl*exp(-sqrt(x.^2+y.^2)/akl)*2.*x/(2*akl*sqrt(x.^2+y.^2)); % x-component of the gradient
    Qy = -Rkl*exp(-sqrt(x.^2+y.^2)/rkl)*2.*y/(2*rkl*sqrt(x.^2+y.^2))+Akl*exp(-sqrt(x.^2+y.^2)/akl)*2.*y/(2*akl*sqrt(x.^2+y.^2)); % y-component of the gradient
    force = [Qx(:),Qy(:)];
end