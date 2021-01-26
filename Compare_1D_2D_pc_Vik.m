%% Today let's put all the data together: 2-D, 1-D, numerics...etc
clear 
close all
   
%% Running the compare functions

D_0=20
%% Running the subfunctions
compare_Viktor(D_0)
compare_PC(D_0)
%% Part 1, from modified_average_rate... from Viktor

load(['2021.1.25_compare_vik,D_0=',num2str(D_0),'.mat'],'AvN_plot','AratesV_plot','rateA_plot','rateAeff_plot','rateN_plot','iF')
% Plotting figures
figure(2);clf;


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
set(gca,'YScale','log')
set(gca,'Xscale','log')


figure(1);clf;
hold on
plot(AvN_plot, 1./AratesV_plot,'x','DisplayName','simulation')
plot(AvN_plot,(1*rateN_plot).*iF((-1 + AvN_plot)),'r','LineStyle','--','DisplayName','numerics')
plot(AvN_plot,(1*rateA_plot).*iF((-1 + AvN_plot)),'g','LineStyle','-.','DisplayName','Kramers')
plot(AvN_plot,(1*rateAeff_plot).*iF((-1 + AvN_plot)),'b','LineStyle',':','DisplayName','KramersEff')
plot(AvN_plot,(2*rateN_plot).*iF((-1 + AvN_plot)),'r','LineStyle','--','DisplayName','numerics2')
plot(AvN_plot,(2*rateA_plot).*iF((-1 + AvN_plot)),'g','LineStyle','-.','DisplayName','Kramers2')
plot(AvN_plot,(2*rateAeff_plot).*iF((-1 + AvN_plot)),'b','LineStyle',':','DisplayName','KramersEff2')
xlabel('$\omega\delta t$','Interpreter','latex');
ylabel('$1/k$','Interpreter','latex');
set(gca,'YScale','log')
set(gca,'Xscale','log')



%% Part 2, from 2-D simulations.

load(['2021.1.25_compare_PC,D_0=',num2str(D_0),'.mat'])
figure(2)
plot(theta_0_matrix,1./trans_rate_matrix_full,'ro')
% legend('simulation','numerics','Kramers','KramersEff','2-D simulation','Location','northwest')

figure(1)
plot(theta_0_matrix,trans_rate_matrix_full,'ro')
% legend('simulation','numerics','Kramers','KramersEff','2-D simulation','Location','southwest')


%% Part 3, from plotting with D_eff measured with Gaussian
transition_rate_measured=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(D_theta_matrix.*theta_0_matrix.^2.*delta_t_matrix));
figure(2)
plot(theta_0_matrix,1./transition_rate_measured,'m--')
plot(theta_0_matrix,1./(2*transition_rate_measured),'m--')
legend('simulation','numerics','Kramers','KramersEff','2-D simulation','Location','northwest')

figure(1)
plot(theta_0_matrix,transition_rate_measured,'m--')
plot(theta_0_matrix,2*transition_rate_measured,'m--')

%% Legending
figure(2)
legend('simulation','numerics','Kramers','KramersEff','2-D simulation','Location','northwest')
figure(1)
legend('simulation','numerics','Kramers','KramersEff','2-D simulation','Location','southwest')
%% Saving Data 

save(['Compare_D_0=',num2str(D_0),'.mat'])
%% Part 3, determine if Viktor and I had the same simulations.
% clear all;
% load('2021.1.25_compare_pc.mat')
% %% 2021.1.18 Plotting Transition Rates of uncorrected and corrected together
%     trans_rate_matrix_full=num_transitions_full_matrix./time_duration;
%     %     trans_rate_matrix_approx=num_transitions_approx_matrix./time_duration;
%     trans_rate_matrix_approx=zeros(size(trans_rate_matrix_full));
% 
%     figure(1);hold on;
%     flag=(theta_0_matrix>1.);
%     plot(theta_0_matrix(flag),transition_rate_approx(flag),'b.')
%     plot(theta_0_matrix(flag),transition_rate_full(flag),'r.')
% %     plot(theta_0_matrix(flag),transition_rate_today(flag),'g-')
%     plot(theta_0_matrix(flag),trans_rate_matrix_full(flag),'xr')
%     plot(theta_0_matrix(flag),trans_rate_matrix_approx(flag),'ob')
%     plot(theta_0_matrix(flag),tranistion_rate_approx_double(flag),'b.')
%     plot(theta_0_matrix(flag),transition_rate_full_double(flag),'r.')
% %     plot(theta_0_matrix(flag),transition_rate_today_double(flag),'g-')
%     title('\theta_0>1')
%     xlabel('\omega_0 \delta t')
%     ylabel('Transition Rates (1/s)')
% %     legend('uncorrected','corrected','today','full 1D simulation','approximated 1D simulation')
% legend('approximated 1D (D_{eff}=4D_0/(\theta_0^2 R^2))','full 1D (D_{eff}=2D_0/R^2)','full 1D simulation','approximated 1D simulation')
% 
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
% %     axis([1.1 1.43 10^-7 10^-1])
%     %% Plotting inverse transition rate of uncorrected and corrected together
%     
%     figure(2);hold on;
%     flag=(theta_0_matrix>1 );
%     plot(theta_0_matrix(flag),1./transition_rate_approx(flag),'.b')
%     plot(theta_0_matrix(flag),1./transition_rate_full(flag),'.r')
% %     plot(theta_0_matrix(flag),1./transition_rate_today(flag),'-g')
%     plot(theta_0_matrix(flag),1./(trans_rate_matrix_full(flag)),'xr')
%     plot(theta_0_matrix(flag),1./(trans_rate_matrix_approx(flag)),'bo')
%     plot(theta_0_matrix(flag),1./tranistion_rate_approx_double(flag),'b.')
%     plot(theta_0_matrix(flag),1./transition_rate_full_double(flag),'.r')
% %     plot(theta_0_matrix(flag),1./transition_rate_today_double(flag),'g-')
%     title('\theta_0>1')
%     xlabel('\omega_0 \delta t')
%     ylabel('1/Transition Rates (s)')
% %     legend('uncorrected','corrected','today','full 1D simulation','approximated 1D simulation')
% legend('approximated 1D (D_{eff}=4D_0/(\theta_0^2 R^2))','full 1D (D_{eff}=2D_0/R^2)','full 1D simulation','approximated 1D simulation')
% 
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
%     axis([1.1 1.43 10^1 10^7])
%     






