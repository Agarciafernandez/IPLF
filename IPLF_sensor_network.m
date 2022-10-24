%Demo with the implementation of the iterated posterior linearisation
%filter (IPLF) proposed in 
% Á. F. García-Fernández, L. Svensson, M. R. Morelande and S. Särkkä, 
% "Posterior Linearization Filter: Principles and Implementation Using Sigma Points," 
% in IEEE Transactions on Signal Processing, vol. 63, no. 20, pp. 5561-5573, Oct.15, 2015,
% This file runs the target tracking scenario with a sensor network in this reference.

% Author: Ángel F. García-Fernández

clear
randn('seed',9)
rand('seed',9)


Scenario_sensor_network;
%Sigma-point parameters
Nx=4;
W0=1/3;
Wn=(1-W0)/(2*Nx);
weights=[W0,Wn*ones(1,2*Nx)];
square_error_t_tot=zeros(1,Nsteps-1);
initial_state_noise=randn(4,Nmc); %To sample the mean at time step 0 at each Monte Carlo run



%N_it: Number of iterations of the IPLF update
N_it=10;
%Threhold to determine convergence based on Kullback-Leibler divergence, it could be set to a negative value to
%always run the algorithm for N_it times
threshold=10^(-5);

N_it_t_tot=zeros(1,Nsteps-1);


randn('seed',9)
rand('seed',9)

for i=1:Nmc
    tic
    %We initiate the mean and covariance matrix at first time step
    Pk=P_ini; 
    meank=X_multi(:,1)+chol_ini*initial_state_noise(:,i);
    square_error_t=zeros(1,Nsteps-1);
    N_it_t=N_it*ones(1,Nsteps-1);
    
    for k=1:Nsteps-1
        %Position and velocity at time step k
        pos_x=X_multi(1,k);
        pos_y=X_multi(3,k);
        vel_x=X_multi(2,k);
        vel_y=X_multi(4,k);
        
        
        %We sample the measurement
        measurement_noise=sqrt_var_z*randn(Nsensors,1);
        d2_s=(pos_x-x_s).^2+(pos_y-y_s).^2;
        
        h_s=SNR0*d0_2./d2_s;
        saturated=d2_s<d0_2;
        h_s(saturated)=SNR0;
        h_s=sqrt(h_s);
        z_real=h_s+measurement_noise; %Real value of the measurement
        
        
        
        
        %IPLF update
        
        mean_pos_j=meank;
        cov_pos_j=Pk;
        
        for p=1:N_it
            
            
            %We calculate the sigma points w.r.t. mean_pos_j and cov_pos_j
            chol_var_mult=chol((Nx/(1-W0)*cov_pos_j));
            sigma_points=[zeros(1,4) ;chol_var_mult(1,:);-chol_var_mult(1,:);chol_var_mult(2,:);-chol_var_mult(2,:);...
                chol_var_mult(3,:);-chol_var_mult(3,:);chol_var_mult(4,:);-chol_var_mult(4,:)]';
            sigma_points=repmat(mean_pos_j,1,length(weights))+sigma_points;
            
            %We transform sigma points
            sigma_points_d2_s=(repmat(sigma_points(1,:),Nsensors,1)-repmat(x_s,1,2*Nx+1)).^2+...
                (repmat(sigma_points(3,:),Nsensors,1)-repmat(y_s,1,2*Nx+1)).^2;
            sigma_points_z=SNR0*d0_2./sigma_points_d2_s;
            saturated=sigma_points_d2_s<d0_2;
            sigma_points_z(saturated)=SNR0;
            
            sigma_points_z=sqrt(sigma_points_z);
            
            %We calculate moments for statistical linear regression
            z_pred=sigma_points_z*weights';
            var_pred=zeros(Nsensors);
            var_xz=zeros(4,Nsensors);
            
            for j=1:length(weights)
                sub_z_j=sigma_points_z(:,j)-z_pred;
                var_pred=var_pred+weights(j)*(sub_z_j*sub_z_j');
                var_xz=var_xz+weights(j)*(sigma_points(:,j)-mean_pos_j)*sub_z_j';
            end
            
            %Statistical linear regression parameters
            
            A_l=var_xz'/cov_pos_j;
            b_l=z_pred-A_l*mean_pos_j;
            Omega_l=var_pred-A_l*cov_pos_j*A_l';
            
            
            %KF moments (with respect to the prior)
            
            z_pred_j=A_l*meank+b_l;
            Psi_j=Pk*A_l';
            
            S_j=A_l*Pk*A_l'+R+Omega_l;
            
            sub_z=z_real-z_pred_j;
            
            mean_updated=meank+Psi_j/S_j*sub_z;
            P_updated=Pk-Psi_j/S_j*Psi_j';
            
            %KLD to measure change in distribution
            dist_k=dist_kullback(mean_pos_j,cov_pos_j,mean_updated,P_updated);
            if(threshold>0 && dist_k<threshold)
                N_it_t(k)=p;
                break;
            end
            
            
            mean_pos_j=mean_updated;
            cov_pos_j=P_updated;
        end
        
        
        square_error_t(k)=(mean_updated(1)-pos_x)^2+(mean_updated(3)-pos_y)^2;

        state_error=mean_updated-X_multi(:,k);
        pos_error=[mean_updated(1)-pos_x;mean_updated(3)-pos_y];
        var_pos_act=P_updated(1:2:4,1:2:4);
                

        
        
        %Predition step
        meank=F*mean_updated;
        Pk=F*P_updated*F'+Q;
        Pk=(Pk+Pk')/2;
        
    end
    square_error_t_tot=square_error_t_tot+square_error_t;
   

    N_it_t_tot=N_it_t_tot+N_it_t;
    
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' seconds'])
end

%Root mean square position error at each time step
RMSE_t_tot=sqrt(square_error_t_tot/Nmc);


%Average number of IPLF iterations to converge
N_it_t_tot=N_it_t_tot/Nmc;


%Uncomment to plot the scenario
% figure(1)
% plot(X_multi(1,1:Nsteps-1),X_multi(3,1:Nsteps-1),'Linewidth',1.3)
% hold on
% plot(X_multi(1,1),X_multi(3,1),'o','Markerfacecolor','b')
% plot(X_multi(1,11:10:Nsteps-1),X_multi(3,11:10:Nsteps-1),'o','Linewidth',1.3)
% plot(x_s,y_s,'xr','Linewidth',1.3,'MarkerSize',10);
% hold off
% grid on
% axis equal
% xlabel('x position (m)')
% ylabel('y position (m)')


figure(2)
plot(RMSE_t_tot)
ylabel ('RMS position error (m)')
xlabel('Time step')
title('IPLF')
grid on

figure(3)
plot(N_it_t_tot,'Linewidth',1.3)
ylabel('Average number of iterations')
xlabel('Time step')
title('IPLF')
grid on
axis([1, Nsteps-1, 0, N_it])


