clear;
close all
Date='2021.1.7' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
nth_take=1 % for a=5, particle 1 fixed
% nth_take=100 % for a=5, particle 1 not fixed
% nth_take=200 % for a=0, particle 1 not fixed
delta_t_matrix=2
T_matrix=[1]
v_0_matrix=[4:0.5:10]
dt=10^-1
intrinsic_delay=0.0 % Intrinsic delay
Obs_time_steps=10^5


%%
a=5
% theta_0_matrix=delta_t_matrix*v_0_matrix/(2*a)
%%
close all
figure(1)
hold on
load(['nth_take=',num2str(16),'.mat']);
figure(1)
plot(theta_0_matrix,D_eff_ratio_matrix)

load(['nth_take=',num2str(115),'.mat']);
figure(1)
plot(theta_0_matrix,D_eff_ratio_matrix)

load(['nth_take=',num2str(215),'.mat']);
figure(1)
plot(theta_0_matrix,D_eff_ratio_matrix)

legend('a=5, particle 1 fixed','a=5, particle 1 mobile','a=0, partilce 1 mobile')

ylabel('D_{eff}/(4D_0/\theta_0^2R^2)')
xlabel('\theta_0')