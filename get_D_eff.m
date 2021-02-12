function [theta,omega,R_mean,D_theta,D_omega2,theta_0,R_recip]= get_D_eff(x,y,v_x,v_y,Obs_time_steps,delta_t,dt,v_0)
%% This file analyzes the x y raw data, gives omega and theta, and analyzes the D_eff (D_theta and D_omega2) 
% modified from test_simulation

%% Calculating angle
% theta(1:Obs_time_steps)=0;
% omega(1:Obs_time_steps)=0;
% R(1:Obs_time_steps)=0;

%% Making length(v_x) == Obs_time_steps
v_x=[v_x 0]; % to keep v_x same length as x
v_y=[v_y 0]; % to keep v_y same length as y


%% 2021.1.15 Vectorized Version
%% Calculating theta
R1=[x(2,1:Obs_time_steps)-x(1,1:Obs_time_steps);y(2,1:Obs_time_steps)-y(1,1:Obs_time_steps);zeros(1,Obs_time_steps)];
R2=[x(2,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt))-x(1,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt));y(2,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt))-y(1,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt));zeros(1,Obs_time_steps)];
diff_x=[x(2,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt))-x(2,1:Obs_time_steps);
y(2,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt))-y(2,1:Obs_time_steps);
zeros(1,Obs_time_steps)];
%             R1_cross_diff_x=cross(R1,diff_x);
R1_cross_diff_x_sign=sign( R1(1,:).*diff_x(2,:)-diff_x(1,:).*R1(2,:));

theta=R1_cross_diff_x_sign.*acos(dot(R1,R2,1)./(vecnorm(R1,2,1).*vecnorm(R2,2,1)));
R=vecnorm(R1,2,1);
%% Calculating omega
omega=(R2(1,:).*v_y(1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt))-R2(2,:).*v_x(1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt)))./vecnorm(R2,2,1).^2;

%% Calculating D_eff
plot_hist_fit_T_eff='no';
%% Finding the sigma of the theta
% [theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0)
f=figure(80);clf;
switch plot_hist_fit_T_eff
    case 'no'
        f.Visible='off';
end
h_theta=histogram(diff(theta));
Values_theta=h_theta.Values;
Bins_theta=h_theta.BinEdges(1:end-1)+0.5*h_theta.BinWidth;

theta_0=v_0*delta_t/mean(R);

[fitresult,gof]=fit_Gaussian(Bins_theta,Values_theta,plot_hist_fit_T_eff);
mu_theta=fitresult.b1;
sigma_theta=fitresult.c1/sqrt(2);

%% 2021.1.12 hist(omega-\pm omega_+)
%% 2021.1.15 hist(omega-omega_0*sin(theta)
% [omega_stable,k_trans_omega,omega_plus,omega_minus,num_transitions_omega]=find_omega_plus(omega,theta_0);
f2=figure(81);clf;
switch plot_hist_fit_T_eff
    case 'no'
        f2.Visible='off';
end
omega_0=v_0/mean(R);
omega2=omega-omega_0*sin(theta);

h_omega2=histogram(omega2);
Values_omega2=h_omega2.Values;
Bins_omega2=h_omega2.BinEdges(1:end-1)+0.5*h_omega2.BinWidth;

[fitresult,gof]=fit_Gaussian(Bins_omega2,Values_omega2,plot_hist_fit_T_eff);

mu_omega2=fitresult.b1;
sigma_omega2=fitresult.c1/sqrt(2);
%% Calculating the D_eff and T_eff
D_theta=sigma_theta^2/(2*dt);
D_omega2=sigma_omega2^2*dt/(2);
%             D_eff_ratio_matrix_approx=[D_eff_ratio_matrix_approx D_theta/(4*D/(theta_0^2*mean(R)^2))];
%             D_eff_ratio_matrix_omega2=[D_eff_ratio_matrix_omega2 D_omega2/(D/mean(R)^2)];
%%
R_mean=mean(R);
R_recip=mean(1./R);
end