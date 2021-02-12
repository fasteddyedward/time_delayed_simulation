%% 2020.1.27 This file is same as Compare_1D_2D_pc_Vik.m, but modified to 
%% scan over omega_0 and delta_t.
clear 
% close all
% cd('dt=0.1')
D_0_matrix=[1]
for D_0_index=1:length(D_0_matrix)
%% Running the compare functions
dt=0.01
D_0=D_0_matrix(D_0_index)
delta_t_matrix_vik=[0.1 0.2 0.5 1 2 5 10]'  %% Make delta_t_matrix as a column vector
delta_t_matrix=[0.1 0.2 0.5 1 2 5 10]'
intrinsic_delay=0.
theta_0_1D_matrix=(logspace(log10(1),log10(3),30));
% theta_0_1D_matrix=(logspace(log10(1),log10(3),15));

%%
file_name_pc=['2021.2.9_scan_pc_D_0=',num2str(D_0)];
file_name_vik=['2021.2.9_scan_vik_D_0=',num2str(D_0)];
%% Running the subfunctions
% compare_Viktor(D_0,dt,tV_max,theta_0_matrix,tau_matrix,runs)
% compare_Viktor(D_0,dt,500,linspace(1.1,1.6,15),tau_matrix,1)
% compare_Viktor(D_0,dt,50000,linspace(1.1,1.6,15),[1.8:0.1:4],1)

Scan_compare_Viktor(D_0,dt,50000,linspace(1.1,1.6,15),delta_t_matrix_vik,1,file_name_vik)


% compare_PC(D_0,dt,delta_t_matrix,v_0_matrix,Obs_time_steps)
% compare_PC(D_0,dt,tau_matrix,[11:0.5:20],10^4)
% compare_PC(D_0,dt,[2:0.2:4],5,10^6)

% Scan_compare_PC(D_0,dt,delta_t_matrix,theta_0_1D_matrix,10^6,file_name_pc)
Scan_compare_PC(D_0,dt,delta_t_matrix,theta_0_1D_matrix,10^6,file_name_pc,intrinsic_delay)

end