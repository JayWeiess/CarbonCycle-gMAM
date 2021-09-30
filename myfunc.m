function [H_1,H_2,H_21,H_22,Lambda,Theta,b,in_a]=myfunc(cx,v)
% generate the functions we need in the gMAM 
% input the parameter c_x and generate the function
% \nabla_x H, \nabla_\theta H, Hess_{\theta,x }H, Hess_{\theta \theta} H
% lambda, \hat{\vartheta}

%indentify the parameter
mu=250; be=4; theta=5; f0=0.694; w0=2000; gamma=4; beta=1.7;cp=110; cf=43.9;

s=@(x,y) x^gamma/(x^gamma+y^gamma);
f=@(x) f0*x.^beta/(x^beta+cf^beta);
ds=@(x,y) gamma*x^(gamma-1)*y^gamma/(x^gamma+y^gamma)^2;
df=@(x) beta*f0*cf^beta*x^(beta-1)/(x^beta+cf^beta)^2;

b=@(x,y) [(mu*(1-be*s(x,cp)-theta*(1-s(x,cx))-v)+y-w0)*f(x); mu*(1-be*s(x,cp)+theta*(1-s(x,cx))+v)-y+w0] ;
a=@(x,y) [mu^2*f(x)^2, 0; 0, mu^2] ;
in_a=@(x,y) [(1/(mu^2*f(x)^2)),0; 0, 1/(mu^2)];

grad_b1=@(x,y) [mu*(-be*ds(x,cp)+theta*ds(x,cx))*f(x)+(mu*(1-be*s(x,cp)-theta*(1-s(x,cx)-v))+y-w0)*df(x);f(x)];               %Caluculation of \nabla b_1
grad_b2=@(x,y) [mu*(-be*ds(x,cp)-theta*ds(x,cx));-1];                                                                                                                         %Caluculation of \nabla b_2
grad_a1=@(x,y) [2*mu^2*f(x)*df(x);0];                                                                                                                                                       %Caluculation of \nabla a_1
grad_a2=@(x,y) [0;0];                                                                                                                                                                                      %Caluculation of \nabla b_1


H_1= @(x,y,u,v) grad_b1(x,y).*u+grad_b2(x,y).*v+1/2*(grad_a1(x,y).*u^2+grad_a2(x,y).*v^2);                   % H_x, input the value of x-(x,y) and \vartheta-(u,v)
H_2= @(x,y,u,v) b(x,y)+a(x,y)*[u;v];                                                                                                                               %H_\theta, input the value of x - (x,y) and  \vartheta-(u,v)
H_21= @(x,y,u,v) [grad_b1(x,y)'+grad_a1(x,y)'.*u ; grad_b2(x,y)'+grad_a2(x,y)'.*u ];                                        %H_{\theta x}, input the value of x - (x,y) and  \vartheta-(u,v)
H_22= @(x,y,u,v) a(x,y);                                                                                                                                                    %H_{\theta \theta}, input the value of x - (x,y) and  \vartheta-(u,v)
Lambda= @(x,y,u,v) sqrt(b(x,y)'*in_a(x,y)*b(x,y))/sqrt([u,v]*in_a(x,y)*[u;v]);                                                                %The calculation of \lambda,  input the value of x - (x,y) and  \vartheta-(u,v)
Theta=@(x,y,u,v) in_a(x,y)*(Lambda(x,y,u,v)*[u;v]-b(x,y));                                                                                        %The calculation of \hat{\vartheta}, input the value of x - (x,y) and  \vartheta-(u,v)