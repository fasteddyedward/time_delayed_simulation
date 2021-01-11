%% This file compares the result from test_simulation_cont.m and test_exp.m
close all;
figure(1)
xlabel('\theta_0')
ylabel('D_{eff}/(D_0/R^2)')
hold on
figure(2)
xlabel('\theta_0')
ylabel('D_{eff}/(4D_0/(\theta_0^2*R^2))')
hold on

clear;
load('experimental_D_eff_ratio.mat')
figure(1)
plot(theta_0_matrix,D_eff_ratio_matrix_full)
figure(2)
plot(theta_0_matrix,D_eff_ratio_matrix_approx)



clear;
load('theoretical_D_eff_ratio.mat')
figure(1)
plot(theta_0_matrix,D_eff_ratio_matrix_full)
legend('exp','sim')
figure(2)
plot(theta_0_matrix,D_eff_ratio_matrix_approx)
legend('exp','sim')