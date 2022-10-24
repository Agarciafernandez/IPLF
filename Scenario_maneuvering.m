


Nsteps=130;
Nmc=1000; %Number of Monte Carlo runs
Nmc_trajectory=50; %Number of Monte Carlo simulations per sampled trajectories

N_trajectories=Nmc/Nmc_trajectory;

%Initial state
x0=[130,25,-20,1,-4*pi/180]';
P_ini=diag([5 5 20000 10 10^(-7)]);
chol_ini=chol(P_ini)';

%Dynamic model
T=0.2;
q1=0.5;
q2=1*10^(-6);

M=[T^3/3 T^2/2; T^2/2 T];

Q=blkdiag(q1*M,q1*M,q2*T);
chol_Q=chol(Q)';


%Measurement model parameters (we have two as there are two types of
%measurements)

var_r=(sqrt(1000))^2;
var_theta=(30*pi/180)^2;
var_fi=(30*pi/180)^2;
var_r_rate=(sqrt(100))^2;

var_r2=(sqrt(1000))^2;
var_theta2=(0.001*pi/180)^2;
var_fi2=(30*pi/180)^2;
var_r_rate2=(sqrt(0.0001))^2;

R1=diag([var_r, var_theta, var_fi, var_r_rate]);
chol_R1=chol(R1)';

R2=diag([var_r2, var_theta2,var_fi2, var_r_rate2]);
chol_R2=chol(R2)';


height=50;


%Nalter -> How many time steps there are between very accurate measurements
Nalter=60;
N0ffset_accurate=2;



%We sample the trajectories
X_multi_series=zeros(5*N_trajectories,Nsteps);

figure(1)
clf
hold on

for i=1:N_trajectories
    x_ini=x0+chol_ini*randn(5,1);
    X_multi_i=Generate_trajectory_turn(x_ini,T,q1,q2,Nsteps);
    X_multi_series(5*i-4:5*i,:)=X_multi_i;
    
    plot(X_multi_i(1,1:Nsteps-1),X_multi_i(3,1:Nsteps-1),'Linewidth',1.3)
    plot(X_multi_i(1,11:10:Nsteps-1),X_multi_i(3,11:10:Nsteps-1),'o','Linewidth',1.3)
end
plot(0,0,'xr','Linewidth',1.3,'MarkerSize',10)

hold off
axis equal
grid on
xlabel('x (m)')
ylabel('y (m)')


noise_z=randn(4,(Nsteps-1)*Nmc); %We obtain the measurement noise for all MC runs here
