%
clf
clc
clearvars
d = 0:0.1:80;

hold on
% xm
Rxm = 10000;
rxm = 15;
Axm = 10;
axm = 10;
% Q = Rkl*exp(-sqrt(d.^2)/rkl)-Akl*exp(-sqrt(d.^2)/akl);
% plot(d,Q);
x = zeros(size(d));
y = d;
Qx = Rxm*exp(-sqrt(x.^2+y.^2)/rxm)*2.*x./(2*rxm*sqrt(x.^2+y.^2))-Axm*exp(-sqrt(x.^2+y.^2)/axm)*2.*x./(2*axm*sqrt(x.^2+y.^2)); % x-component of the gradient
Qy = -Rxm*exp(-sqrt(x.^2+y.^2)/rxm)*2.*y./(2*rxm*sqrt(x.^2+y.^2))+Axm*exp(-sqrt(x.^2+y.^2)/axm)*2.*y./(2*axm*sqrt(x.^2+y.^2));
Qxy = Qy;%sqrt(Qx.^2+Qy.^2);
plot(d,Qxy)

Rmm = 5000;
rmm = 13;
Amm = 10;
amm = 10;
% Q = Rkl*exp(-sqrt(d.^2)/rkl)-Akl*exp(-sqrt(d.^2)/akl);
% plot(d,Q);
Qx = Rmm*exp(-sqrt(x.^2+y.^2)/rmm)*2.*x./(2*rmm*sqrt(x.^2+y.^2))-Amm*exp(-sqrt(x.^2+y.^2)/amm)*2.*x./(2*amm*sqrt(x.^2+y.^2)); % x-component of the gradient
Qy = -Rmm*exp(-sqrt(x.^2+y.^2)/rmm)*2.*y./(2*rmm*sqrt(x.^2+y.^2))+Amm*exp(-sqrt(x.^2+y.^2)/amm)*2.*y./(2*amm*sqrt(x.^2+y.^2));
Qxy = Qy;%sqrt(Qx.^2+Qy.^2);
plot(d,Qxy)


Rmx = 3000;
rmx = 20;
Amx = 1500;
amx = 3;
% Q = Rkl*exp(-sqrt(d.^2)/rkl)-Akl*exp(-sqrt(d.^2)/akl);
% plot(d,Q);
Qx = Rmx*exp(-sqrt(x.^2+y.^2)/rmx)*2.*x./(2*rmx*sqrt(x.^2+y.^2))-Amx*exp(-sqrt(x.^2+y.^2)/amx)*2.*x./(2*amx*sqrt(x.^2+y.^2)); % x-component of the gradient
Qy = -Rmx*exp(-sqrt(x.^2+y.^2)/rmx)*2.*y./(2*rmx*sqrt(x.^2+y.^2))+Amx*exp(-sqrt(x.^2+y.^2)/amx)*2.*y./(2*amx*sqrt(x.^2+y.^2));
Qxy = Qy;%sqrt(Qx.^2+Qy.^2);
plot(d,Qxy)


Rxx = 2000;
rxx = 10;
Axx = 10;
axx = 10;

% Q = Rkl*exp(-sqrt(d.^2)/rkl)-Akl*exp(-sqrt(d.^2)/akl);
% plot(d,Q);


Qx = Rxx*exp(-sqrt(x.^2+y.^2)/rxx)*2.*x./(2*rxx*sqrt(x.^2+y.^2))-Axx*exp(-sqrt(x.^2+y.^2)/axx)*2.*x./(2*axx*sqrt(x.^2+y.^2)); % x-component of the gradient
Qy = -Rxx*exp(-sqrt(x.^2+y.^2)/rxx)*2.*y./(2*rxx*sqrt(x.^2+y.^2))+Axx*exp(-sqrt(x.^2+y.^2)/axx)*2.*y./(2*axx*sqrt(x.^2+y.^2));
Qxy = Qy;%sqrt(Qx.^2+Qy.^2);
plot(d,Qxy)
    
hold off