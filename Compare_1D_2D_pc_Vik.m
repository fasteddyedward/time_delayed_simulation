%% Today let's put all the data together: 2-D, 1-D, numerics...etc
clear 
% close all
% cd('dt=0.1')
%% Running the compare functions
dt=0.01
D_0=2

%% Running the subfunctions
% compare_Viktor(D_0,dt,tV_max,theta_0_matrix,tau_matrix,runs)
% compare_Viktor(D_0,dt,50000,linspace(1.1,1.6,15),tau_matrix,1)
compare_Viktor(D_0,dt,50000,linspace(1.1,1.6,15),[1.8:0.1:4],1)

% compare_PC(D_0,dt,delta_t_matrix,v_0_matrix,Obs_time_steps)
% compare_PC(D_0,dt,tau_matrix,[11:0.5:20],10^4)
compare_PC(D_0,dt,[2:0.2:4],5,10^6)
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




%% Part 2, from 2-D simulations.
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


%% Saving Data 
% save(['Compare_D_0=',num2str(D_0),'.mat'])

%% Part 6, determine if Viktor and I had the same simulations for 1-D
% 
% load('2021.1.21_compare_pc_1D.mat')
% %% 2021.1.18 Plotting Transition Rates of uncorrected and corrected together
% transition_rate_approx=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(4*D_0./R_matrix.^2.*delta_t_matrix));
% % transition_rate_full=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(D_theta_full_matrix.*delta_t_matrix.*theta_0_matrix.^2));
% transition_rate_full=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(2*D_0./R_matrix.^2.*delta_t_matrix.*theta_0_matrix.^2));
% 
% 
% tranistion_rate_approx_double=2*transition_rate_approx;
% transition_rate_full_double=2*transition_rate_full;
% %%
%     figure(1);hold on;
%     flag=(theta_0_matrix>1.);
% 
% %     plot(theta_0_matrix(flag),trans_rate_matrix_full(flag),'xr')
%     plot(theta_0_matrix,transition_rate_full,'.k')
%     plot(theta_0_matrix,2*transition_rate_full,'.k')
%     plot(theta_0_matrix,transition_rate_approx,'.b')
%     plot(theta_0_matrix,2*transition_rate_approx,'.b')
%     title('\theta_0>1')
%     xlabel('\omega_0 \delta t')
%     ylabel('Transition Rates (1/s)')
% 
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
%     
%     set(gca,'XScale','linear')
%     set(gca,'YScale','linear')
%     %% Plotting inverse transition rate of uncorrected and corrected together
%     
%     figure(2);hold on;
%     flag=(theta_0_matrix>1 );
%     plot(theta_0_matrix(flag),1./(trans_rate_matrix_full(flag)),'xr')
%     plot(theta_0_matrix,1./transition_rate_full,'.k')
%     plot(theta_0_matrix,1./(2*transition_rate_full),'.k')
%     plot(theta_0_matrix,1./transition_rate_approx,'.b')
%     plot(theta_0_matrix,1./(2*transition_rate_approx),'.b')
%     title('\theta_0>1')
%     xlabel('\omega_0 \delta t')
%     ylabel('1/Transition Rates (s)')
% 
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
%     
%     set(gca,'XScale','linear')
%     set(gca,'YScale','linear')

