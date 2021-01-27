%% 2020.1.27 This file is same as Compare_1D_2D_pc_Vik.m, but modified to 
%% scan over omega_0 and delta_t.
clear 
% close all
% cd('dt=0.1')
%% Running the compare functions
dt=0.01
D_0=200
delta_t_matrix=[1:2]
v_0_matrix=[10:5:20]

%% Running the subfunctions
% compare_Viktor(D_0,dt,tV_max,theta_0_matrix,tau_matrix,runs)
% compare_Viktor(D_0,dt,500,linspace(1.1,1.6,15),tau_matrix,1)
% compare_Viktor(D_0,dt,50000,linspace(1.1,1.6,15),[1.8:0.1:4],1)
Scan_compare_Viktor(D_0,dt,50000,linspace(1.1,1.6,15),delta_t_matrix,1)


% compare_PC(D_0,dt,delta_t_matrix,v_0_matrix,Obs_time_steps)
% compare_PC(D_0,dt,tau_matrix,[11:0.5:20],10^4)
% compare_PC(D_0,dt,[2:0.2:4],5,10^6)
% Scan_compare_PC(D_0,dt,delta_t_matrix,v_0_matrix,10^3)
%%
% load(['2021.1.25_compare_PC,D_0=',num2str(D_0),'.mat'])

% load(['2021.1.25_compare_vik,D_0=',num2str(D_0),'.mat'],'AvN_plot','AratesV_plot','rateA_plot','rateAeff_plot','rateN_plot','iF')
load(['2021.1.25_compare_vik,D_0=',num2str(D_0),'.mat'])
%% Part 1, from modified_average_rate... from Viktor

load(['2021.1.25_compare_vik,D_0=',num2str(D_0),'.mat'],'AvN_plot','AratesV_plot','rateA_plot','rateAeff_plot','rateN_plot','iF')
% Plotting figures
figure(2);clf;title(['1-D, D_0=',num2str(D_0)])


hold on
plot(AvN_plot, AratesV_plot,'x','DisplayName','simulation')
plot(AvN_plot,1./(1*rateN_plot).*iF((-1 + AvN_plot)),'r','LineStyle','--','DisplayName','numerics')
plot(AvN_plot,1./(1*rateA_plot).*iF((-1 + AvN_plot)),'g','LineStyle','-.','DisplayName','Kramers')
plot(AvN_plot,1./(1*rateAeff_plot).*iF((-1 + AvN_plot)),'b','LineStyle',':','DisplayName','KramersEff')
plot(AvN_plot,1./(2*rateN_plot).*iF((-1 + AvN_plot)),'r','LineStyle','--','DisplayName','numerics2')
plot(AvN_plot,1./(2*rateA_plot).*iF((-1 + AvN_plot)),'g','LineStyle','-.','DisplayName','Kramers2')
plot(AvN_plot,1./(2*rateAeff_plot).*iF((-1 + AvN_plot)),'b','LineStyle',':','DisplayName','KramersEff2')
xlabel('$\omega\delta t$','Interpreter','latex');
ylabel('$1/k$','Interpreter','latex');



figure(1);clf;title(['1-D, D_0=',num2str(D_0)])
hold on
plot(AvN_plot, 1./AratesV_plot,'x','DisplayName','simulation')
plot(AvN_plot,(1*rateN_plot).*iF((-1 + AvN_plot)),'r','LineStyle','--','DisplayName','numerics')
plot(AvN_plot,(1*rateA_plot).*iF((-1 + AvN_plot)),'g','LineStyle','-.','DisplayName','Kramers')
plot(AvN_plot,(1*rateAeff_plot).*iF((-1 + AvN_plot)),'b','LineStyle',':','DisplayName','KramersEff')
plot(AvN_plot,(2*rateN_plot).*iF((-1 + AvN_plot)),'r','LineStyle','--','DisplayName','numerics2')
plot(AvN_plot,(2*rateA_plot).*iF((-1 + AvN_plot)),'g','LineStyle','-.','DisplayName','Kramers2')
plot(AvN_plot,(2*rateAeff_plot).*iF((-1 + AvN_plot)),'b','LineStyle',':','DisplayName','KramersEff2')
xlabel('$\omega\delta t$','Interpreter','latex');
ylabel('$k$','Interpreter','latex');


%% Part 1.2, plotting Kramer_eff and Kramer for 2-D, modified for Overleaf
figure(1);
set(gca,'YScale','log')
set(gca,'Xscale','log')
% legend('1-D simulation','2-D k_{D_{eff}}','2-D k{_D_{eff}} * 2','2-D k_{4D_0/(R^2\theta_0^2)}','2-D k_{4D_0/(R^2\theta_0^2)} * 2','2-D k_{2D_0/R^2}','2-D k_{2D_0/R^2} * 2','Location','southwest')
% legend('simulation','numerics','Kramers','KramersEff','numerics*2','Kramers*2','KramersEff*2','2-D simulation','2-D measured','2-D measured*2','Location','northwest')
legend('1-D simulation','1-D numerics', '1-D k_{4D_0/(R^2\theta_0^2)}','1-D k_{2D_0/(R^2)}','1-D numerics * 2','1-D k_{4D_0/(R^2\theta_0^2)} * 2','1-D k_{2D_0/(R^2)} * 2','Location','southwest')

figure(2)
set(gca,'YScale','log')
set(gca,'Xscale','log')
% legend('2-D simulation','2-D k_{D_{eff}}','2-D k{_D_{eff}} * 2','2-D k_{4D_0/(R^2\theta_0^2)}','2-D k_{4D_0/(R^2\theta_0^2)} * 2','2-D k_{2D_0/R^2}','2-D k_{2D_0/R^2} * 2','Location','northwest')
legend('1-D simulation','1-D numerics', '1-D k_{4D_0/(R^2\theta_0^2)}','1-D k_{2D_0/(R^2)}','1-D numerics * 2','1-D k_{4D_0/(R^2\theta_0^2)} * 2','1-D k_{2D_0/(R^2)} * 2','Location','northwest')

%% Part 2, from 2-D simulations from Pin Chuan 
load(['2021.1.25_compare_PC,D_0=',num2str(D_0),'.mat'])
trans_rate_sim_2D=num_transitions_matrix./time_duration;
figure(2)
plot(theta_0_matrix,1./trans_rate_sim_2D,'ko')
figure(1)
plot(theta_0_matrix,trans_rate_sim_2D,'ko')



%% Part 3, from plotting with D_eff measured with Gaussian
trans_rate_measured=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(D_theta_matrix.*theta_0_matrix.^2.*delta_t_matrix));
figure(2)
plot(theta_0_matrix,1./trans_rate_measured,'m--')
plot(theta_0_matrix,1./(2*trans_rate_measured),'m--')
figure(1)
plot(theta_0_matrix,trans_rate_measured,'m--')
plot(theta_0_matrix,2*trans_rate_measured,'m--')

%% Legending
figure(2)
axis([1 inf -inf inf])
set(gca,'YScale','log')
set(gca,'Xscale','log')
%     set(gca,'XScale','linear')
%     set(gca,'YScale','linear')
legend('simulation','numerics','Kramers','KramersEff','numerics*2','Kramers*2','KramersEff*2','2-D simulation','2-D measured','2-D measured*2','Location','northwest')

figure(1)
axis([1 inf -inf inf])
set(gca,'YScale','log')
set(gca,'Xscale','log')
%     set(gca,'XScale','linear')
%     set(gca,'YScale','linear')
legend('simulation','numerics','Kramers','KramersEff','numerics*2','Kramers*2','KramersEff*2','2-D simulation','2-D measured','2-D measured*2','Location','southwest')

%% Part 4, Plotting D_eff measured and 2*D_0/R^2
figure(3);clf;hold on
plot(D_theta_matrix./(2*D_0./R_matrix.^2))
% plot(D_theta_matrix./(2*D_0.*R_recip_matrix.^2))
axis([-inf inf 0 inf])
ylabel('$D_{eff}/(2D_0/R^2)$','Interpreter','latex');
legend('2D_0/R_{mean}^2','2D_0(1/R)_{mean}^2')

figure(4);clf;hold on
plot(D_theta_matrix./(D_0./R_matrix.^2))
% plot(D_theta_matrix./(2*D_0.*R_recip_matrix.^2))
axis([-inf inf 0 inf])
xlabel('$\theta_0$', 'Interpreter','latex');
ylabel('$D_{eff}/(D_0/R^2)$','Interpreter','latex');
% legend('2D_0/R_{mean}^2','2D_0(1/R)_{mean}^2')


%% Part 5, plotting Kramer_eff and Kramer for 2-D

figure(8);clf;hold on;title(['2-D, D_0=',num2str(D_0)])
plot(theta_0_matrix,1./trans_rate_sim_2D,'ko')
plot(theta_0_matrix,1./trans_rate_measured,'g--')
plot(theta_0_matrix,1./(2*trans_rate_measured),'g--')

figure(7);clf;hold on;title(['2-D, D_0=',num2str(D_0)])
plot(theta_0_matrix,trans_rate_sim_2D,'ko')
plot(theta_0_matrix,trans_rate_measured,'g--')
plot(theta_0_matrix,2*trans_rate_measured,'g--')



% Using R_mean_matrix
transition_rate_Kramer=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(4*D_0./R_matrix.^2.*delta_t_matrix));
transition_rate_Kramer_eff=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(2*D_0./R_matrix.^2.*delta_t_matrix.*theta_0_matrix.^2));
% Using 1/R_recip_matrix
% transition_rate_Kramer=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(4*D_0.*R_recip_matrix.^2.*delta_t_matrix));
% transition_rate_Kramer_eff=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(2*D_0.*R_recip_matrix.^2.*delta_t_matrix.*theta_0_matrix.^2));
figure(8)
plot(theta_0_matrix,1./transition_rate_Kramer,'r--')
plot(theta_0_matrix,1./(2*transition_rate_Kramer),'r--')
plot(theta_0_matrix,1./transition_rate_Kramer_eff,'b:')
plot(theta_0_matrix,1./(2*transition_rate_Kramer_eff),'b:')
figure(7)
plot(theta_0_matrix,transition_rate_Kramer,'r--')
plot(theta_0_matrix,2*transition_rate_Kramer,'r--')
plot(theta_0_matrix,transition_rate_Kramer_eff,'b:')
plot(theta_0_matrix,(2*transition_rate_Kramer_eff),'b:')

figure(7)
set(gca,'YScale','log')
set(gca,'Xscale','log')
legend('2-D simulation','2-D measured','2-D measured*2','2-D Kramer','2-D Kramer*2','2-D KramerEff','2-D KramerEff*2','Location','southwest')
figure(8)
set(gca,'YScale','log')
set(gca,'Xscale','log')

legend('2-D simulation','2-D measured','2-D measured*2','2-D Kramer','2-D Kramer*2','2-D KramerEff','2-D KramerEff*2','Location','northwest')

%% Part 5.2, plotting Kramer_eff and Kramer for 2-D, modified for Overleaf


figure(7)
set(gca,'YScale','log')
set(gca,'Xscale','log')
legend('2-D simulation','2-D k_{D_{eff}}','2-D k{_D_{eff}} * 2','2-D k_{4D_0/(R^2\theta_0^2)}','2-D k_{4D_0/(R^2\theta_0^2)} * 2','2-D k_{2D_0/R^2}','2-D k_{2D_0/R^2} * 2','Location','southwest')
% legend('2-D simulation','2-D measured','2-D measured*2','2-D Kramer','2-D Kramer*2','2-D KramerEff','2-D KramerEff*2','Location','southwest')
figure(8)
set(gca,'YScale','log')
set(gca,'Xscale','log')
legend('2-D simulation','2-D k_{D_{eff}}','2-D k{_D_{eff}} * 2','2-D k_{4D_0/(R^2\theta_0^2)}','2-D k_{4D_0/(R^2\theta_0^2)} * 2','2-D k_{2D_0/R^2}','2-D k_{2D_0/R^2} * 2','Location','northwest')
% legend('2-D simulation','2-D measured','2-D measured*2','2-D Kramer','2-D Kramer*2','2-D KramerEff','2-D KramerEff*2','Location','northwest')

%% Saving Data 
% save(['Compare_D_0=',num2str(D_0),'.mat'])
