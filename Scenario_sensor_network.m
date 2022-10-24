clear

%We load the trajectory
load('Trajectory')


%Dynamic model
T=0.5;
sigmaQ=0.2;
F=[1 T 0 0;0 1 0 0; 0 0 1 T;0 0 0 1];
Q=sigmaQ^2*kron(eye(2),[T^3/3 T^2/2; T^2/2 T]);



Nsteps=160;
Nmc=1000; %Number of Monte Carlo runs

%Initial variace at time step 0
P_ini=diag([49 4 1 2]);
chol_ini=chol(P_ini)';


%Measurement model parameters

var_z=0.1;
sqrt_var_z=sqrt(var_z);

%Spacing between sensors
delta=25;
%Positions of sensors
x_lin=0:delta:100;
y_lin=0:delta:100;
[x_s,y_s]=meshgrid(x_lin,y_lin);
x_s=x_s(:);
y_s=y_s(:);

Nsensors=length(x_s);


SNR0=1000;
d0_2=(1)^2; %d0 squared
R=diag(var_z*ones(1,Nsensors));


