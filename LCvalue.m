function [u_stable,u]=LCvalue(x1,v)
%% Compute a sample of the original trajectory
%% initialization
mu = 250;
b = 4;
theta = 5;
c_x = 57;
c_p = 110;
nu =v;
y0 = 2000;
gama = 4;
f_0 = 0.694;
c_f = 43.9;
beta = 1.70;

%% carbon system: vector filed (f,g)
f = @(x,y)    (mu*(1-b*x^gama/(x^gama+c_p^gama)-theta*(1-x^gama/(x^gama+c_x^gama))-nu)+y-y0)*f_0*x^beta/(x^beta+c_f^beta);
g = @(x,y)    mu*(1-b*x^gama/(x^gama+c_p^gama)+theta*(1-x^gama/(x^gama+c_x^gama))+nu)-y+y0;

%% Compute the trajectory using Euler scheme
T = 1e2;
dd = 1e5;
dt = T/dd;
% t = 0 : dt: T;
x=zeros(1,dd+1);
y=zeros(1,dd+1);

x_stable=(b-1)^(-1/gama)*c_p;
y_stable=y0+mu*(theta+nu-(theta*c_p^gama)/((b-1)*c_x^gama+c_p^gama));
u_stable=[x_stable;y_stable];

 x(1)=x1(1);
 y(1)=x1(2);
 
for i = 1 : dd
    x(i+1) = x(i) + f(x(i),y(i))*dt;
    y(i+1) = y(i) + g(x(i),y(i))*dt;
end
u=[x(dd+1);y(dd+1)];


