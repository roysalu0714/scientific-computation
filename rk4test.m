% test rk4 function
F_xy = @(t,r) 3.*exp(-t)-0.4*r;
a=0;
b=3;
y0=5;
iternum=10;
[y,x]=rk4(F_xy,a,b,y0,iternum);

plot(x,y)