function X_multi=Generate_trajectory_turn(x0,T,q1,q2,Nsteps)

% Author: Ángel F. García-Fernández

% This code samples a random trajectory from a turning motion model 

M=[T^3/3 T^2/2; T^2/2 T];

Q=blkdiag(q1*M,q1*M,q2*T);

chol_Q=chol(Q)';

X_multi=zeros(5,Nsteps);
X_multi(:,1)=x0;
xk=x0;



for k=2:Nsteps
    Omega=xk(5);
    F=[1 sin(Omega*T)/Omega 0 -(1-cos(Omega*T))/Omega 0;...
        0 cos(Omega*T) 0 -sin(Omega*T) 0;...
        0 (1-cos(Omega*T))/Omega 1 sin(Omega*T)/Omega 0;...
        0 sin(Omega*T) 0 cos(Omega*T) 0;...
        0 0 0 0 1];
         
    xk1=F*xk+chol_Q*randn(5,1);
    
    xk=xk1;
    
    X_multi(:,k)=xk;
    
end

% figure(1)
% plot(X_multi(1,:),X_multi(3,:))
% axis equal
% grid on
    