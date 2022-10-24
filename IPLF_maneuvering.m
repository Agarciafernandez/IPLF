% Demo with the implementation of the iterated posterior linearisation
% filter (IPLF) proposed in 
% [1] Á. F. García-Fernández, L. Svensson, M. R. Morelande and S. Särkkä, 
% "Posterior Linearization Filter: Principles and Implementation Using Sigma Points," 
% in IEEE Transactions on Signal Processing, vol. 63, no. 20, pp. 5561-5573, Oct.15, 2015,
% This file runs the target tracking scenario with a maneuvering target and
% with a sensor that measures range, bearings, elevation and range rate

% NOTE: In this code and in the example in [1], we model a bearings
% measurement with additive Gaussian noise. In this case, if bearings
% are close to the -pi boundary, we may get errors with standard
% application of Kalman filters (iterated or not). 
% A principled approach to solve this problem is to model bearings
% measurements using circular statistics such as the von-Mises Fisher distribution
% or the Kent distribution. The IPLF for these types of distributions (as
% well as unscented/extended Kalman filters) are explained in 

% [2] Á. F. García-Fernández, F. Tronarp and S. Särkkä, "Gaussian Target Tracking With Direction-of-Arrival von Mises–Fisher Measurements," 
% in IEEE Transactions on Signal Processing, vol. 67, no. 11, pp. 2960-2972, 1 June1, 2019

% [3] . F. García-Fernández, S. Maskell, P. Horridge and J. Ralph, "Gaussian Tracking With Kent-Distributed Direction-of-Arrival Measurements," 
% in IEEE Transactions on Vehicular Technology, vol. 70, no. 7, pp. 7249-7254, July 2021

% Matlab code of these filters are found at https://github.com/Agarciafernandez/DOA-tracking


% A Less principled solution that can provide reasonable results consists of using circular sums and subtractions
% in the Kalman filter update (or IPLF update) as explained in 

% [4] D. F. Crouse, “Cubature/unscented/sigma point Kalman filtering with
% angular measurement models," in Proc. 18th Int. Conf. Inf. Fusion, 2015,
% pp. 1550–1557

% and used in

% Á. F. Garcia-Fernandez, M. R. Morelande and J. Grajal, "Truncated Unscented Kalman Filtering," 
% in IEEE Transactions on Signal Processing, vol. 60, no. 7, pp. 3372-3386, July 2012

% This code uses the standard sum and subtractions in the IPLF update so it will not provide
% suitable results for certain trajectories (adding the circular sum and
% subtractions in [4] is straightforward).


% Author: Ángel F. García-Fernández

clear
randn('seed',9)
rand('seed',9)

Scenario_maneuvering;

%Sigma-point parameters
Nx=5;
W0=1/3; %Weight at the mean

Wn=(1-W0)/(2*Nx);
weights=[W0,Wn*ones(1,2*Nx)];
square_error_t_tot=zeros(1,Nsteps-1);




rms_t_series=zeros(1,Nmc);


%N_it: Number of iterations of the IPLF update
N_it=10;
%Threhold to determine convergence based on Kullback-Leibler divergence, it could be set to a negative value to
%always run the algorithm for N_it times
threshold=10^(-5);

N_it_t_tot=zeros(1,Nsteps-1);

N_finaliza=0;
N_converge=0;
N_finaliza_series=zeros(1,Nmc);
N_converge_series=zeros(1,Nmc);


randn('seed',9)
rand('seed',9)

for i=1:Nmc
    tic
    Pk=P_ini;
    
    n_trajectory=fix((i-1)/Nmc_trajectory)+1; %Trajectory number
    X_multi=X_multi_series(5*n_trajectory-4:5*n_trajectory,:);    
    meank=x0;
    
    square_error_t=zeros(1,Nsteps-1);
   
    
    
    N_it_t=N_it*ones(1,Nsteps-1);
    
    for k=1:Nsteps-1
        pos_x=X_multi(1,k);
        pos_y=X_multi(3,k);
        vel_x=X_multi(2,k);
        vel_y=X_multi(4,k);
        omega_k=X_multi(5,k);
        
        
        %We check if we make a measurement of type 1 or type 2
        if(mod(k,Nalter)==N0ffset_accurate)
            noise=chol_R2*noise_z(:,k+(Nsteps-1)*(i-1));
        else
            noise=chol_R1*noise_z(:,k+(Nsteps-1)*(i-1));
        end
        %We obtain the measurement
        z_real=[sqrt(pos_x^2+pos_y^2+height^2);...
            atan2(pos_y,pos_x);...
            atan2(height,sqrt(pos_x^2+pos_y^2));...
            (pos_x*vel_x+pos_y*vel_y)/sqrt(pos_x^2+pos_y^2+height^2)]+noise;
        
        
        
        
        %IPLF update
        
        mean_pos_j=meank;
        cov_pos_j=Pk;
        
        
        for p=1:N_it
            
            
            %We calculate the sigma points w.r.t. mean_pos_j and cov_pos_j
            chol_var_mult=chol((Nx/(1-W0)*cov_pos_j));
            sigma_points=[zeros(1,5) ;chol_var_mult(1,:);-chol_var_mult(1,:);chol_var_mult(2,:);-chol_var_mult(2,:);...
                chol_var_mult(3,:);-chol_var_mult(3,:);chol_var_mult(4,:);-chol_var_mult(4,:);chol_var_mult(5,:);-chol_var_mult(5,:)]';
            sigma_points=repmat(mean_pos_j,1,length(weights))+sigma_points;

            %We transform sigma points
            sigma_points_r=sqrt(sigma_points(1,:).^2+sigma_points(3,:).^2+height^2);
            sigma_points_z=[sigma_points_r;...
                atan2(sigma_points(3,:),sigma_points(1,:));...
                atan2(height,sqrt(sigma_points(1,:).^2+sigma_points(3,:).^2));...
                (sigma_points(1,:).*sigma_points(2,:)+sigma_points(3,:).*sigma_points(4,:))./sigma_points_r];
            
                    
            z_pred=sigma_points_z*weights';
            
            
            var_pred=zeros(4);
            var_xz=zeros(5,4);
            
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
            
            %We add the covariance matrix of the corresponding measurement
            %noise
            
            if(mod(k,Nalter)==N0ffset_accurate)
                S_j=A_l*Pk*A_l'+R2+Omega_l;
            else
                S_j=A_l*Pk*A_l'+R1+Omega_l;
            end
            
            
           
            sub_z=z_real-z_pred_j;
            
            mean_updated=meank+Psi_j/S_j*sub_z;
            var_updated=Pk-Psi_j/S_j*Psi_j';
            
            %Diferencia KLD
            dist_k=dist_kullback(mean_pos_j,cov_pos_j,mean_updated,var_updated);
            
            if(threshold>0 && dist_k<threshold)
                N_it_t(k)=p;
                break;
            end
            
            N_it_t(k)=p;
            mean_pos_j=mean_updated;
            cov_pos_j=var_updated;
        end
        
  
        
        square_error_t(k)=(mean_updated(1)-pos_x)^2+(mean_updated(3)-pos_y)^2;
       
        
        %Predition step
        
        %We calculate the sigma points
        chol_var_mult=chol((Nx/(1-W0)*var_updated));
        sigma_points=[zeros(1,5) ;chol_var_mult(1,:);-chol_var_mult(1,:);chol_var_mult(2,:);-chol_var_mult(2,:);...
            chol_var_mult(3,:);-chol_var_mult(3,:);chol_var_mult(4,:);-chol_var_mult(4,:);chol_var_mult(5,:);-chol_var_mult(5,:)]';
        sigma_points=repmat(mean_updated,1,length(weights))+sigma_points;
        
        
        %We propagate the sigma points
        sigma_points_pred=zeros(size(sigma_points));
        var_pred=zeros(5,5);
        
        for j=1:length(weights)
            Omega=sigma_points(5,j);
            F=[1 sin(Omega*T)/Omega 0 -(1-cos(Omega*T))/Omega 0;...
                0 cos(Omega*T) 0 -sin(Omega*T) 0;...
                0 (1-cos(Omega*T))/Omega 1 sin(Omega*T)/Omega 0;...
                0 sin(Omega*T) 0 cos(Omega*T) 0;...
                0 0 0 0 1];
            
            sigma_points_pred(:,j)=F*sigma_points(:,j);
            var_pred=var_pred+weights(j)*sigma_points_pred(:,j)*sigma_points_pred(:,j)';
            
        end
        
        %The predicted mean and covariance matrix are
        mean_pred=sigma_points_pred*weights';
        var_pred=var_pred-mean_pred*mean_pred'+Q; 
        
        
        meank=mean_pred;
        Pk=var_pred;
        Pk=(Pk+Pk')/2;
     
        
        
    end
    square_error_t_tot=square_error_t_tot+square_error_t;
    N_it_t_tot=N_it_t_tot+N_it_t;


    
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' seconds'])
    
end

%Mean square error at each time step
square_error_t_tot=square_error_t_tot/Nmc;

%Root mean square position error at each time step
RMSE_t_tot=sqrt(square_error_t_tot/Nmc);

%Root mean square error across all time steps
RMSE_tot=sqrt(sum(square_error_t_tot)/(Nsteps-1))

%Average number of IPLF iterations to converge
N_it_t_tot=N_it_t_tot/Nmc;



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
grid on
axis([1, Nsteps-1, 0, N_it])


