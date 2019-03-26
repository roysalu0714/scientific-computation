function [y,x]=rk4(f,a,b,y0,iternum)
% Runge-Kutta 4th order methods for ODE initial value problems
% input
% f ODE function
% a,b interval
% y0 initial conditions
% iternum iterative times
% output
% y 
% x 
%
h=(b-a)/iternum;
y=zeros(iternum+1,length(y0));
x=a:h:b;
y(1,:)=y0;
for ii=1:iternum
    k1=f(x(ii),y(ii,:))';
    k2=f(x(ii)+h/2,y(ii,:)+h*k1/2)';
    k3=f(x(ii)+h/2,y(ii,:)+h*k2/2)';
    k4=f(x(ii)+h,y(ii,:)+h*k3)';
    y(ii+1,:)=y(ii,:)+h*(k1+2*k2+2*k3+k4)/6;
end