function dydt=CbC(~,u)

x=u(1);
y=u(2);
nu =0.01;

mu = 250;
b = 4;
theta = 5;
c_x = 62;
c_p = 110;
y0 = 2000;
gama = 4;
f_0 = 0.694;
c_f = 43.9;
beta = 1.70;

%% carbon system: vector filed (f,g)
f =     (mu*(1-b*x^gama/(x^gama+c_p^gama)-theta*(1-x^gama/(x^gama+c_x^gama))-nu)+y-y0)*f_0*x^beta/(x^beta+c_f^beta);
g =     mu*(1-b*x^gama/(x^gama+c_p^gama)+theta*(1-x^gama/(x^gama+c_x^gama))+nu)-y+y0;

dydt=[f;g];