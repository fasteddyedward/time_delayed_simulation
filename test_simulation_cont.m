clear;
close all

% Date='2021.1.8' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% % nth_take=1 % for a=5, particle 1 fixed
% % nth_take=100 % for a=5, particle 1 not fixed
% nth_take=200 % for a=0, particle 1 not fixed
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[4:0.5:10]
% dt=10^-3
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6

% Date='2021.1.8' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=300 % for a=5, particle 1 fixed
% nth_take=400 % for a=5, particle 1 not fixed
% nth_take=500 % for a=0, particle 1 not fixed
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[0.5:0.5:16]
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5

%%
a=5
% theta_0_matrix=delta_t_matrix*v_0_matrix/(2*a)
%%
close all
figure(1)
hold on
figure(2)
hold on

% load(['nth_take=',num2str(215),'.mat']);
% load(['nth_take=',num2str(469),'.mat']);
load(['nth_take=',num2str(342),'.mat']);
figure(1)
plot(theta_0_matrix,D_eff_ratio_matrix_full)
figure(2)
plot(theta_0_matrix,D_eff_ratio_matrix_approx)

% load(['nth_take=',num2str(469),'.mat']);
% figure(1)
% plot(theta_0_matrix,D_eff_ratio_matrix_full)
% figure(2)
% plot(theta_0_matrix,D_eff_ratio_matrix_approx)
close all
save('theoretical_D_eff_ratio.mat')
%% 2020.1.7  
% load(['nth_take=',num2str(16),'.mat']);
% figure(1)
% plot(theta_0_matrix,D_eff_ratio_matrix_full)
% figure(2)
% plot(theta_0_matrix,D_eff_ratio_matrix_approx)
% 
% load(['nth_take=',num2str(115),'.mat']);
% figure(1)
% plot(theta_0_matrix,D_eff_ratio_matrix_full)
% figure(2)
% plot(theta_0_matrix,D_eff_ratio_matrix_approx)
% 
% 
% load(['nth_take=',num2str(215),'.mat']);
% figure(1)
% plot(theta_0_matrix,D_eff_ratio_matrix_full)
% figure(2)
% plot(theta_0_matrix,D_eff_ratio_matrix_approx)

% legend('a=5, particle 1 fixed','a=5, particle 1 mobile','a=0, partilce 1 mobile')
% 
% ylabel('D_{eff}/(4D_0/\theta_0^2R^2)')
% xlabel('\theta_0')