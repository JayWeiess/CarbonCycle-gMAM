%% Compute the most likely transition path for a SDE

clear 
clc
% define the matrix function H_{\theta x},  H_x, H_\theta, \lambda,
% \hat{\varphi}, b(u,v)
% store the value of \varphi as a column vector
load('LimitCycle_value_nu=0.1.mat')
cx=62;
v=0.1;
[H_1,H_2,H_21,H_22,Lambda,Theta,b,in_a]=myfunc(cx,v);
%tau=1000;                                                                                                        %Define 1/\tau to be 10

K=1;                                       %The iteration times 
N=3000;                                      % The partition size
%T=1;                                            %The original time scale
tau=3000;
eta=0.1;
r=1;
C=1;
Action_mini=10;
P_mini=zeros(2,N+1);
idex=1;

i_start=1;
i_interval=10;
i_end=11;

P4=zeros(2*(i_end-i_start)/i_interval,N+1);
Actions=zeros((i_end-i_start)/i_interval,1);
iters=zeros((i_end-i_start)/i_interval,1);
js=zeros((i_end-i_start)/i_interval,1);

P=zeros(2,N+1);                      %The matrix that stores the function {\varphi_i^k}, where the (:,i+1)-th element of this matrix stores the value of {\varphi_i}
 %Define the end points
 [x1,~]=LCvalue([84;2500],v);
%x1=[83.5920067307416;2315.60787129321];  

 fprintf('%5s %15s %15s %15s\n', 'iter','j','Action','Action_mini');
for iter=i_start:i_interval:i_end
%for iter=1:1

%x1=[83;2800];
%x2=LCV(:,iter);     
x2=LCV(:,iter);     
%x2=[9.73348; 3153.56];
P=zeros(2,N+1); 
P1=zeros(2,N-1);                     % The matrix that stores  {\varphi_i^{\prime k}}
P2=zeros(2,N+1);                    % The matrix that stores  {\tilde{\varphi}_i}
P3=zeros(2,N+1);                    % An intermediate matrix
d=(x2(1)-x1(1)) /N;
 % generate a line with ends x1 and x2 as the initial value
P(1,:)=x1(1):d:x2(1);               
P(2,:)=interp1([x1(1),x2(1)],[x1(2),x2(2)],P(1,:));
%,     ,
d=2/N;

%% k given
%for j=1:K
ActionTemp=0;
theta=zeros(2,N-1);
lambda=zeros(1,N+1);
lambda1=zeros(1,N-1);
 for i=1:N-1
 P1(:, i)=(P(:,i+2)-P(:,i))/d;                                                                 %The calculation of  \varphi_i^{\prime k}
 
 if norm(P1(:,i))<=1e-9
     lambda(i+1)=0;
     theta(:,i)=[0;0];
 else
     theta(:,i)=Theta(P(1,i+1),P(2,i+1),P1(1,i),P1(2,i));                                                     %The calculation of \hat{\vartheta}_i^k
     lambda(i+1)=H_2(P(1,i+1),P(2,i+1),theta(1,i),theta(2,i))'*P1(:,i)/norm(P1(:,i))^2;   %The calculation of \lambda_i^k, the (i+1)-th element of this matrix stores the value of \lambda_i^k
    % lambda(i+1)=Lambda(P(1,i+1),P(2,i+1),P1(1,i),P1(2,i));
     %x=P(1,i+1); y=P(2,i+1); u=P1(1,i); v=P1(2,i);
     %lambda(i+1)=(b(x,y)'*in_a(x,y)*b(x,y))/([u,v]*in_a(x,y)*[u;v]);  
    % c(i)=b(x,y)'*in_a(x,y)*b(x,y);
 end
 end
lambda(1)=3*lambda(2)-3*lambda(3)+lambda(4);                       % define \lambda_0
lambda(N+1)=3*lambda(N)-3*lambda(N-1)+lambda(N-2);      % define \lambda_N
for i=1:N-1
lambda1(i)=(lambda(i+2)-lambda(i))/d;                                           %The calculation of \lambda^{\prime}
end
Action=3/2*P1(:,1)'*theta(:,1)+3/2*P1(:,N-1)'*theta(:,N-1);
for i=2:N-2
    Action=P1(:,i)'*theta(:,i)+Action;
end
Action=Action/N;
%% Solve the linear system Ax=b

A1=zeros(1,2*(N+1));                                                                           %define the elements of the matrix
A2=zeros(1,2*N-2);
A3=zeros(1,2*N+2);

for i=1:N-1
    A1(2*i+1:2*(i+1))=[-2*(lambda(i+1)*N)^2-tau,-2*(lambda(i+1)*N)^2-tau];
    A2(2*i-1:2*i)=[(lambda(i+1)*N)^2,(lambda(i+1)*N)^2] ;
    A3(2*i+1:2*i+2)=(lambda(i+1)*H_21(P(1,i+1),P(2,i+1),theta(1,i),theta(2,i))*P1(:,i)-H_22(P(1,i+1),P(2,i+1),theta(1,i),theta(2,i))*H_1(P(1,i+1),P(2,i+1),theta(1,i),theta(2,i))-lambda(i+1)*lambda1(i)*P1(:,i)-tau*P(:,i+1))';
end
A1(1:2)=[1,1];
A1(2*N+1:2*N+2)=[1,1];
A3(1:2)=x1'; A3(2*N+1:2*N+2)=x2';

v=[A1,A2,A2];                                                                                                            %Store the elements value of the sparse matrix S
row=[1:1:2*N+2, 3:1:2*N,3:1:2*N];                                                                     %Store the row index of the sparse matrix S
col=[1:1:2*N+2,1:1:2*N-2,5:1:2*N+2];                                                              %Store the coloumn index of the sparse matrix S
A=sparse(row,col,v);                                                                                                %Define the sparse matrix A
S=A\A3';                                                                                                                     %Solve the linear equation;
for i=1:N+1
    P2(:,i)=[S(2*i-1,1);S(2*i,1)];                                                                                 %Store the value of \tilde{\varphi}_i^k in P2
end

% %% Reparametrization though equidistance difference
% Length=sum(sqrt((P2(1,1:N)-P2(1,2:N+1)).^2+(P2(2,1:N)-P2(2,2:N+1)).^2));                                                              %Calculate the length of the interpolated curve
% delta=Length/N;
% [vec,res,k]=findp(P2(:,1),P2(:,2),delta,0);
% i=1;
% P3(:,1)=x1;
% if k~=0  
% i=i+1;
% P3(:,2:k+1)=vec;
% else
%     while k==0
%         i=i+1;
%         [vec,res,k]=findp(P2(:,i),P2(:,i+1),delta,res);
%      end
%     P3(:,2:k+1)=vec;
% end
% m=k+1;
%   while i<=N
%        [vec,res,k]=findp(P2(:,i),P2(:,i+1),delta,res);
%        if k~=0
%            P3(:,m+1:m+k)=vec;
%            m=m+k;
%        end
%       i=i+1;
%   end
%  %P3(:,N+1)=x2;
%  P=P3;
%  P(:,1)=x1;
%  P(:,N+1)=x2;
 
 %% Reparametrization though equi-real time 

one=ones(1,N+1);
T=zeros(2,N+1);
T1=zeros(2,N+1);
T(:,1)=[0;0];
T(:,2)=[1/N;1/(2*max(lambda(1),eta))+1/(2*max(lambda(2),eta))];
for i=3:N+1
    T(1,i)=(i-1)/N;
    T(2,i)=1/(2*max(eta,lambda(1)))+sum(one(2:i-1)./max(lambda(2:i-1),eta))+1/(2*max(eta,lambda(i)));
end
T1(1,:)=0:T(2,N+1)/N:T(2,N+1);
T1(2,:)=interp1(T(2,:),T(1,:),T1(1,:));
P3(1,:)=interp1(0:1/N:1,P2(1,:),T1(2,:));
P3(2,:)=interp1(0:1/N:1,P2(2,:),T1(2,:));
P=P3;
P(:,1)=x1;
P(:,N+1)=x2;
j=1;
while abs(Action-ActionTemp)>1e-7 && j<=K 
 ActionTemp=Action;
 j=j+1;
 theta=zeros(2,N-1);
lambda=zeros(1,N+1);
lambda1=zeros(1,N-1);
 for i=1:N-1
 P1(:, i)=(P(:,i+2)-P(:,i))/d;                                                                 %The calculation of  \varphi_i^{\prime k}
 
 if norm(P1(:,i))<=1e-9
     lambda(i+1)=0;
     theta(:,i)=[0;0];
 else
     theta(:,i)=Theta(P(1,i+1),P(2,i+1),P1(1,i),P1(2,i));                                                     %The calculation of \hat{\vartheta}_i^k
     lambda(i+1)=H_2(P(1,i+1),P(2,i+1),theta(1,i),theta(2,i))'*P1(:,i)/norm(P1(:,i))^2;   %The calculation of \lambda_i^k, the (i+1)-th element of this matrix stores the value of \lambda_i^k
    % lambda(i+1)=Lambda(P(1,i+1),P(2,i+1),P1(1,i),P1(2,i));
     %x=P(1,i+1); y=P(2,i+1); u=P1(1,i); v=P1(2,i);
     %lambda(i+1)=(b(x,y)'*in_a(x,y)*b(x,y))/([u,v]*in_a(x,y)*[u;v]);  
    % c(i)=b(x,y)'*in_a(x,y)*b(x,y);
 end
 end
lambda(1)=3*lambda(2)-3*lambda(3)+lambda(4);                       % define \lambda_0
lambda(N+1)=3*lambda(N)-3*lambda(N-1)+lambda(N-2);      % define \lambda_N
for i=1:N-1
lambda1(i)=(lambda(i+2)-lambda(i))/d;                                           %The calculation of \lambda^{\prime}
end
Action=3/2*P1(:,1)'*theta(:,1)+3/2*P1(:,N-1)'*theta(:,N-1);
for i=2:N-2
    Action=P1(:,i)'*theta(:,i)+Action;
end
Action=Action/N;
if abs(Action-0.0150)<1e-6
    break
end
%% Solve the linear system Ax=b

A1=zeros(1,2*(N+1));                                                                           %define the elements of the matrix
A2=zeros(1,2*N-2);
A3=zeros(1,2*N+2);

for i=1:N-1
    A1(2*i+1:2*(i+1))=[-2*(lambda(i+1)*N)^2-tau,-2*(lambda(i+1)*N)^2-tau];
    A2(2*i-1:2*i)=[(lambda(i+1)*N)^2,(lambda(i+1)*N)^2] ;
    A3(2*i+1:2*i+2)=(lambda(i+1)*H_21(P(1,i+1),P(2,i+1),theta(1,i),theta(2,i))*P1(:,i)-H_22(P(1,i+1),P(2,i+1),theta(1,i),theta(2,i))*H_1(P(1,i+1),P(2,i+1),theta(1,i),theta(2,i))-lambda(i+1)*lambda1(i)*P1(:,i)-tau*P(:,i+1))';
end
A1(1:2)=[1,1];
A1(2*N+1:2*N+2)=[1,1];
A3(1:2)=x1'; A3(2*N+1:2*N+2)=x2';

v=[A1,A2,A2];                                                                                                            %Store the elements value of the sparse matrix S
row=[1:1:2*N+2, 3:1:2*N,3:1:2*N];                                                                     %Store the row index of the sparse matrix S
col=[1:1:2*N+2,1:1:2*N-2,5:1:2*N+2];                                                              %Store the coloumn index of the sparse matrix S
A=sparse(row,col,v);                                                                                                %Define the sparse matrix A
S=A\A3';                                                                                                                     %Solve the linear equation;
for i=1:N+1
    P2(:,i)=[S(2*i-1,1);S(2*i,1)];                                                                                 %Store the value of \tilde{\varphi}_i^k in P2
end


%% Reparametrization though equidistance difference
Length=sum(sqrt((P2(1,1:N)-P2(1,2:N+1)).^2+(P2(2,1:N)-P2(2,2:N+1)).^2));                                                              %Calculate the length of the interpolated curve
delta=Length/N;
[vec,res,k]=findp(P2(:,1),P2(:,2),delta,0);
i=1;
P3(:,1)=x1;
if k~=0  
i=i+1;
P3(:,2:k+1)=vec;
else
    while k==0
        i=i+1;
        [vec,res,k]=findp(P2(:,i),P2(:,i+1),delta,res);
     end
    P3(:,2:k+1)=vec;
end
m=k+1;
  while i<=N
       [vec,res,k]=findp(P2(:,i),P2(:,i+1),delta,res);
       if k~=0
           P3(:,m+1:m+k)=vec;
           m=m+k;
       end
      i=i+1;
  end
 P3(:,N+1)=x2;
 P=P3;


% %% Reparametrization though equi-real time 
% %lambda=sum(abs(b(P2(1,:),P2(2,:))).^r,1);
% lambda=abs(lambda).^r;
% one=ones(1,N+1);
% T=zeros(2,N+1);
% T1=zeros(2,N+1);
% T(:,1)=[0;0];
% T(:,2)=[1/N;C/(2*max(lambda(1),eta))+C/(2*max(lambda(2),eta))];
% for i=3:N+1
%     T(1,i)=(i-1)/N;
%     T(2,i)=C/(2*max(eta,lambda(1)))+sum(C.*one(2:i-1)./max(lambda(2:i-1),eta))+C/(2*max(eta,lambda(i)));
% end
% T1(1,:)=0:T(2,N+1)/N:T(2,N+1);
% T1(2,:)=interp1(T(2,:),T(1,:),T1(1,:));
% P3(1,:)=interp1(0:1/N:1,P2(1,:),T1(2,:));
% P3(2,:)=interp1(0:1/N:1,P2(2,:),T1(2,:));
% P=P3;
% P(:,1)=x1;
% P(:,N+1)=x2;

end
if Action < Action_mini
    Action_mini=Action;
    P_mini=P;
    idex=iter;
end

P4(2*(iter-i_start)/i_interval+1:2*(iter-i_start)/i_interval+2,:)=P(:,1:N+1);
Actions((iter-i_start)/i_interval+1,1)=Action; 
iters((iter-i_start)/i_interval+1,1)=iter;
js((iter-i_start)/i_interval+1,1)=j;

 fprintf('%5d %15.4f %15.4e %15.4e\n',iter,j,Action,Action_mini);
 end



%  figure(1)
%  plot(P(1,2:N),P1(2,:));
% xlabel('x')
% ylabel('$\varphi^\prime$','interpreter','latex','FontSize', 12)
% title(' Plot of $\varphi^\prime$','interpreter','latex', 'FontSize', 18)

  figure(2)
  plot(P_mini(1,:),P_mini(2,:),'b');
  hold on 
  plot(LCV(1,:),LCV(2,:), 'r')
% %   hold on 
% %   plot(P(1,:),P(2,:),'go');
%   xlabel('x')
  ylabel('$\varphi$','interpreter','latex','FontSize', 12)
  title(' Plot of $\varphi$','interpreter','latex', 'FontSize', 18)
  
%   figure(3)
%   plot(P(1,1:N+1),lambda)
%   xlabel('x')
%    ylabel('$\lambda$','interpreter','latex','FontSize', 12)
%   title(' Plot of $\lambda$','interpreter','latex', 'FontSize', 18)
save('Results_nu=0.1_iter=30000_N=3000_tau=3000_no=1-10-3831.mat','Action_mini','Actions','P4','P_mini','idex','js');

