%% compute the limit cycle in two dimension

clear,clc;
tic

%% initialization
mu = 250;
b = 4;
theta = 5;
c_x = 57;
c_p = 110;
nu = 0;
y0 = 2000;
gama = 4;
f_0 = 0.694;
c_f = 43.9;
beta = 1.70;

%% carbon system: vector filed (f,g)
f = @(x,y)    [mu*(1-b*x^gama/(x^gama+c_p^gama)-theta*(1-x^gama/(x^gama+c_x^gama))-nu)+y-y0]*f_0*x^beta/(x^beta+c_f^beta);
g = @(x,y)    mu*(1-b*x^gama/(x^gama+c_p^gama)+theta*(1-x^gama/(x^gama+c_x^gama))+nu)-y+y0;

%% Compute the trajectory using Euler scheme
T = 100;
N = 1e5;
dt = T/N;
t = 0 : dt: T;
%83.5819254216627;2290.50827400743
%[83.6339868155324;2290.44660430012]

 x(1) =80;
 y(1) =3000;
%92.6715592143817;2281.91346961487
for i = 1 : N
    x(i+1) = x(i) + f(x(i),y(i))*dt;
    y(i+1) = y(i) + g(x(i),y(i))*dt;
end
  plot(x,y,'b')
  hold on
  stable=[x(1:5000);y(1:5000)];
%  plot(t,x,'r');hold on
%  plot(t,y,'g')
p = 1;
LCV(1,1) = x(end);
LCV(2,1) = y(end);
while p > 0 || i==0
    i = i - 1;
    LCV(1,N-i) = x(i+1);
    LCV(2,N-i) = y(i+1);
    if norm(LCV(:,1)-LCV(:,N-i))<5 && t(end)-t(i)>0.5
       p = 0;
    end
end
LCV(:,N-i+1)=[x(end);y(end)];
plot(LCV(1,:),LCV(2,:),'r')
%save('LimitCycle.txt','u','-ascii','-double');
%save('LimitCycle_value_nu=0.67.mat','LCV');
toc