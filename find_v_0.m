clear
%% This file tries to find the range of v_0_matrix for Scan_Compare_1D_2D.m
Date='1';
Obs_time_steps=10^4;

T_matrix=1;
dt=0.01;
nth_take=1;
intrinsic_delay=0;
D_0=1
delta_t_matrix=0.2
v_0_matrix=linspace(0,200,10)
% [num_transitions_matrix,theta_plus_matrix,theta_minus_matrix,R_matrix,D_omega_matrix,D_theta_matrix,theta_0_matrix,R_recip_matrix,v_0_matrix]= main(Date,nth_take,delta_t_matrix,T_matrix,v_0_matrix,dt,intrinsic_delay,Obs_time_steps,D_0);
[num_transitions_matrix,theta_plus_matrix,theta_minus_matrix,R_matrix,D_omega_matrix,D_theta_matrix,theta_0_matrix,R_recip_matrix]= main(Date,nth_take,delta_t_matrix,T_matrix,v_0_matrix,dt,intrinsic_delay,Obs_time_steps,D_0);
v_0_matrix
theta_0_matrix
%%
figure(2);clf;
hold on

plot(v_0_matrix,theta_0_matrix)
yline(pi/2)
yline(1)

%% Results
% D_0=1
% theta_0 ; 1 ; 1.5
% delta_t ; v_0_min ; ;v_0_max
% 0.1 ; 100; 150
% 0.2 ; 50 ; 75
% 0.5 ; 23 ; 40 
% 1 ; 11 ; 20
% 2 ; 6 ; 15
% 5 ; 3 ; 10
% 10 ; 1.5; 4

% v_0=theta_0/delta_t*10; 

%% The v_0_min are about inverse to delta_t

% D_0=10
% 1 ; 12; 52