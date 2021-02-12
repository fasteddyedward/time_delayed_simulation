
%% Same as test.m. Functionalized to compare the transition rate with Viktor's in Compare_1D_2D_pc_Vik.m
function Scan_compare_PC(D_0,dt,delta_t_matrix,theta_0_1D_matrix,Obs_time_steps,file_name_pc,intrinsic_delay)
close all



Date='2021.1.25';
nth_take=1;
% delta_t_matrix=1
T_matrix=[1];

% v_0_matrix=[11:0.5:20]
% dt=10^-2
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6

%% Running main.m
[num_transitions_matrix,theta_plus_matrix,theta_minus_matrix,R_matrix,D_omega_matrix,D_theta_matrix,theta_0_matrix,R_recip_matrix,v_0_matrix]= Scan_main(Date,nth_take,delta_t_matrix,T_matrix,theta_0_1D_matrix,dt,intrinsic_delay,Obs_time_steps,D_0);
%% Running parameters
Bifurcation_Diagram='no';

%% Transition Rates
% omega_0_matrix=v_0_matrix./R_matrix;
time_duration=Obs_time_steps.*dt+delta_t_matrix;
% theta_0_matrix=omega_0_matrix.*delta_t_matrix;
%% Comparing E_b and k_B*T
k_B=1;
save([file_name_pc,'.mat'])
end
% 


