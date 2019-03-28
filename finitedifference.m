h=1/121;
epsilon=0.01;

l_0=0;
l=1;

alpha=1;
beta=1;

N=(l-l_0)/h;

B=zeros(N-1,1);
A=zeros(N-1,N-1);

a_i=zeros(1,N-1);
b_i=zeros(1,N-1);
c_i=zeros(1,N-1);
f_i=zeros(1,N-1);

for i=1:N-1
    a_i(1,i)=-2+h^2*1/epsilon*(-1);
    b_i(1,i)=1-h/2*1/epsilon*(i/N*(l-l_0))^2*(-1);
    c_i(1,i)=1+h/2*1/epsilon*(i/N*(l-l_0))^2*(-1);
end

k1=b_i(1,1);
b_i(1,1:N-2)=b_i(1,2:N-1);
b_i(1,N-1)=k1;

k2=c_i(1,N-1);
c_i(1,2:N-1)=c_i(1,1:N-2);
c_i(1,1)=k2;

f_i(1,1)=-alpha*b_i(1,1);
f_i(1,N-1)=-beta*c_i(1,N-1);

B=f_i;

A = spdiags([b_i',a_i',c_i'],[-1,0,1],N-1,N-1);

y=A\B';

x=l_0:h:l;

plot(x(2:N),y')

