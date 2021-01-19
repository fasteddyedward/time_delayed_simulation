%% This file runs main.m and analyzes the obtained information
close all
clear;

%%
% Date='2021.1.15' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=1 % for a=5, particle 1 fixed
% delta_t_matrix=[0.5:0.5:4.0]
% % delta_t_matrix=[10]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-2
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6

Date='2021.1.6'
nth_take=100
delta_t_matrix=2
T_matrix=[1]
% v_0_matrix=[0.5:0.5:14]
v_0_matrix=[5:0.5:8]
dt=10^-1
intrinsic_delay=0.0 % Intrinsic delay
Obs_time_steps=10^5

% Date='2021.1.18' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=1 % for a=5, particle 1 fixed
% delta_t_matrix=[0.5:0.5:4.0]
% % delta_t_matrix=[10]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-2
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6

% Date='2021.1.18'
% nth_take=1
% delta_t_matrix=2
% T_matrix=[1]
% % v_0_matrix=[0.5:0.5:14]
% v_0_matrix=[5:0.1:8]
% dt=10^-2
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^8


% Date='2021.1.18' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=1 % for a=5, particle 1 fixed
% delta_t_matrix=[0.5:0.1:4.0]
% % delta_t_matrix=[10]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-2
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^8

%% File_name
save_file_name=['Date=',Date,', delta_t=',num2str(min(delta_t_matrix)),':',num2str(max(delta_t_matrix)),', v_0=',num2str(min(v_0_matrix)),':',num2str(max(v_0_matrix)),', dt=',num2str(dt),', time=',num2str(Obs_time_steps),'.mat']
load(save_file_name)
%% Running main.m
[num_transitions_matrix,theta_plus_matrix,theta_minus_matrix,R_matrix,D_omega_matrix,D_theta_matrix,theta_0_matrix]= main(Date,nth_take,delta_t_matrix,T_matrix,v_0_matrix,dt,intrinsic_delay,Obs_time_steps);
%% Running parameters
Bifurcation_Diagram='yes'

%% The Bifurcation Diagram
switch Bifurcation_Diagram
    case 'yes'
% v_0_matrix(1)=[];
    %% Delta_t
    if length(delta_t_matrix)>1
        v_0=v_0_matrix;
        T=T_matrix;
        
        close all
        %% Bifurcation Diagram
        figure(2) ;clf
        hold on
        plot(v_0*delta_t_matrix./R_matrix,theta_plus_matrix,'.')
        plot(v_0*delta_t_matrix./R_matrix,theta_minus_matrix,'.')
        title(['Bifurcation Diagram, v_0= ',num2str(v_0),', T=',num2str(T)])
        xlabel('\omega_0 \delta t=v_0*\delta t/R')
        ylabel('\theta (rad)')
        theta_theory=@(x)real(sqrt(10-sqrt(1/42./x.^6+120./x-20)));
        theta_line=0:0.01:1.8;
        plot(theta_line,theta_theory(theta_line),'k')
        plot(theta_line,-theta_theory(theta_line),'k')
        legend('$\theta_+$','$\theta_-$','$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','','interpreter','latex','Location','northwest')
%         saveas(gcf,['Bifurcation Diagram, v_0= ',num2str(v_0),', T=',num2str(T),'.png'])
        
    end


    %% V_0
    if length(v_0_matrix)>1
        delta_t=delta_t_matrix;
        T=T_matrix;
        close all
        %% Bifurcation Diagram
        figure(2) ;clf
        hold on
        plot(v_0_matrix*delta_t./(R_matrix),theta_plus_matrix,'.')
        plot(v_0_matrix*delta_t./(R_matrix),theta_minus_matrix,'.')
        title(['Bifurcation Diagram, \delta t= ',num2str(delta_t),', T=',num2str(T)])
        xlabel('\omega_0 \delta t=v_0*\delta t/R')
        ylabel('\theta (rad)')
        theta_theory=@(x)real(sqrt(10-sqrt(1/42./x.^6+120./x-20)));
        theta_line=0:0.01:1.8;
        plot(theta_line,theta_theory(theta_line),'k')
        plot(theta_line,-theta_theory(theta_line),'k')
        legend('$\theta_+$','$\theta_-$','$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','','interpreter','latex','Location','northwest')


    end
end

%% Transition Rates
D=1;
notation='theta'
omega_0_matrix=v_0_matrix./R_matrix;
time_duration=Obs_time_steps.*dt+delta_t_matrix;

% this (uncorrected) transition rate is completely calculated with Kramer's theory
tranistion_rate_uncorrected=2*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(4*D./R_matrix.^2.*delta_t_matrix));
% this (today,2021.1.19) transition is an alternative rate suggested by Viktor for my
% 2-D simulation, which is very close to corrected rate and was inspired by
% D_theta/(D/R^2)~2
transition_rate_today=2*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(2*D./R_matrix.^2.*theta_0_matrix.^2.*delta_t_matrix));

% these corrected transition rates use the effective diffusion
switch notation
    case 'omega'
        transition_rate_theory=2*sqrt(2)./(pi*omega_0_matrix.*delta_t_matrix.^2).*(omega_0_matrix.*delta_t_matrix-1).*exp(-3*(omega_0_matrix.*delta_t_matrix-1).^2./(D_theta_matrix.*omega_0_matrix.^2.*delta_t_matrix.^5));
    case 'theta'
        theta_0_matrix=omega_0_matrix.*delta_t_matrix;
        transition_rate_theory=2*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(D_theta_matrix.*delta_t_matrix.*theta_0_matrix.^2));
end
%% 2021.1.18 Plotting Uncorrected Transition Rates
    trans_rate_matrix=num_transitions_matrix./time_duration;
    figure(25);clf;hold on;
    plot(theta_0_matrix,tranistion_rate_uncorrected,'-')
    plot(theta_0_matrix,trans_rate_matrix,'x')
    title('Full range of \theta_0')
    xlabel('\omega_0 \delta t')
    ylabel('Uncorrected Transition Rates (1/s)')


    figure(26);clf;hold on;
    flag=(theta_0_matrix>1.);
    plot(theta_0_matrix(flag),tranistion_rate_uncorrected(flag),'-')

    plot(theta_0_matrix(flag),trans_rate_matrix(flag),'x')
    title('\theta_0>1')
    xlabel('\omega_0 \delta t')
    ylabel('Uncorrected Transition Rates (1/s)')

    %% Plotting inverse transition rate
    figure(27);clf;hold on;
    plot(theta_0_matrix,1./tranistion_rate_uncorrected,'-')
    plot(theta_0_matrix,1./(trans_rate_matrix),'x')
    title('Full range of \theta_0')
    xlabel('\omega_0 \delta t')
    ylabel('Uncorrected 1/Transition Rates (s)')


    figure(28);clf;hold on;
    flag=(theta_0_matrix>1 );
    plot(theta_0_matrix(flag),1./tranistion_rate_uncorrected(flag),'-')
    plot(theta_0_matrix(flag),1./(trans_rate_matrix(flag)),'x')
    title('\theta_0>1')
    xlabel('\omega_0 \delta t')
    ylabel('Uncorrected 1/Transition Rates (s)') 
        
%% 2020.12.10 Plotting Corrected Transition Rates 
    trans_rate_matrix=num_transitions_matrix./time_duration;
    figure(7);clf;hold on;
    plot(theta_0_matrix,transition_rate_theory,'-')
    plot(theta_0_matrix,trans_rate_matrix,'x')
    title('Full range of \theta_0')
    xlabel('\omega_0 \delta t')
    ylabel('Corrected Transition Rates (1/s)')


    figure(8);clf;hold on;
    flag=(theta_0_matrix>1.);
    plot(theta_0_matrix(flag),transition_rate_theory(flag),'-')

    plot(theta_0_matrix(flag),trans_rate_matrix(flag),'x')
    title('\theta_0>1')
    xlabel('\omega_0 \delta t')
    ylabel('Corrected Transition Rates (1/s)')

    %% Plotting inverse transition rate
    figure(11);clf;hold on;
    plot(theta_0_matrix,1./transition_rate_theory,'-')
    plot(theta_0_matrix,1./(trans_rate_matrix),'x')
    title('Full range of \theta_0')
    xlabel('\omega_0 \delta t')
    ylabel('Corrected 1/Transition Rates (s)')


    figure(12);clf;hold on;
    flag=(theta_0_matrix>1 );
    plot(theta_0_matrix(flag),1./transition_rate_theory(flag),'-')
    plot(theta_0_matrix(flag),1./(trans_rate_matrix(flag)),'x')
    title('\theta_0>1')
    xlabel('\omega_0 \delta t')
    ylabel('Corrected 1/Transition Rates (s)')
%% 2021.1.18 Plotting Transition Rates of uncorrected and corrected together
    trans_rate_matrix=num_transitions_matrix./time_duration;
    figure(7);clf;hold on;
    plot(theta_0_matrix,tranistion_rate_uncorrected,'-')
    plot(theta_0_matrix,transition_rate_theory,'-')
    plot(theta_0_matrix,transition_rate_today,'-')
    plot(theta_0_matrix,trans_rate_matrix,'x')
    title('Full range of \theta_0')
    xlabel('\omega_0 \delta t')
    ylabel('Transition Rates (1/s)')
    legend('uncorrected','corrected','today','simulation')
%     set(gca,'YScale','log')

    figure(8);clf;hold on;
    flag=(theta_0_matrix>1.);
    plot(theta_0_matrix(flag),tranistion_rate_uncorrected(flag),'-')
    plot(theta_0_matrix(flag),transition_rate_theory(flag),'-')
    plot(theta_0_matrix(flag),transition_rate_today(flag),'-')
    plot(theta_0_matrix(flag),trans_rate_matrix(flag),'x')
    title('\theta_0>1')
    xlabel('\omega_0 \delta t')
    ylabel('Transition Rates (1/s)')
    legend('uncorrected','corrected','today','simulation')

    %% Plotting inverse transition rate
    figure(11);clf;hold on;
    plot(theta_0_matrix,1./(tranistion_rate_uncorrected),'-')
    plot(theta_0_matrix,1./transition_rate_theory,'-')
    plot(theta_0_matrix,1./transition_rate_today,'-')
    plot(theta_0_matrix,1./(trans_rate_matrix),'x')
    title('Full range of \theta_0')
    xlabel('\omega_0 \delta t')
    ylabel('1/Transition Rates (s)')
    legend('uncorrected','corrected','today','simulation')
    

    figure(12);clf;hold on;
    flag=(theta_0_matrix>1 );
    plot(theta_0_matrix(flag),1./tranistion_rate_uncorrected(flag),'-')
    plot(theta_0_matrix(flag),1./transition_rate_theory(flag),'-')
    plot(theta_0_matrix(flag),1./transition_rate_today(flag),'-')
    plot(theta_0_matrix(flag),1./(trans_rate_matrix(flag)),'x')
    title('\theta_0>1')
    xlabel('\omega_0 \delta t')
    ylabel('1/Transition Rates (s)')
    legend('uncorrected','corrected','today','simulation')
%     set(gca,'YScale','log')
%% Comparing between 1-D and 2-D % It's the same, so just for check actually
% figure(13);clf;hold on;
% % title('R vs 2a')
% % omega_0_matrix=V_0_matrix./(R_matrix);
% gamma_eff_matrix=omega_0_matrix.*Delta_t_matrix.^2/2;
% T_eff_1dim=gamma_eff_matrix.*D_eff_matrix/k_B;
% ratio=T_eff_matrix./T_eff_1dim;
% plot(x0,ratio)
% xlabel('\omega_0 \delta t')


%% Comparing E_b and k_B*T
k_B=1
figure(14);clf;
T_eff_matrix=D_theta_matrix;
plot(theta_0_matrix,-3/2*(omega_0_matrix.*delta_t_matrix-1).^2./(omega_0_matrix.*k_B.*T_eff_matrix.*delta_t_matrix.^3))
xlabel('\omega_0 \delta_t')
ylabel('E_b/k_B T')
%% Plotting R vs omega_0 delta_t

% figure(10);clf;hold on;
% errorbar(delta_t_matrix*v_0_matrix./R_matrix , R_matrix/(2*a), std_dR_matrix/(2*a));
% yline(1)
% xlabel('\omega_0 \delta t ')
% ylabel('R/(2a)')
% axis([-inf inf 0 inf])
% legend('(std of R)/2a','R=2a','Location','northwest')
%% Plotting std_R and std_dR
% figure(22); clf; hold on;
% title('Fluctuations on the Radial Direction')
% plot(theta_0_matrix,std_R_matrix)
% plot(theta_0_matrix,std_dR_matrix)
% xlabel('omega_0 \delta t')
% legend('r', 'dr')

%% Plotting T_eff and D_eff
figure(15);clf;
plot(theta_0_matrix,T_eff_matrix)
xlabel('\omega_0 \delta_t')
ylabel('T_{eff}')
figure(16);clf;
plot(theta_0_matrix,D_theta_matrix)
xlabel('\omega_0 \delta_t')
ylabel('D_{eff}')

%% Plotting T_eff and D_eff & R vs omega_0 * delta_t
% figure(21);clf;hold on
% % title('fixed \delta t')
% % title('fixed \omega_0')
% plot(theta_0_matrix,T_eff_matrix)
% plot(theta_0_matrix,D_theta_matrix)
% plot(theta_0_matrix,R_matrix/(2*a)/300)
% xlabel('\omega_0 \delta t')
% legend('T_{eff}','D_{eff}','R/(2a)/300','Location','northwest')

%% Plotting D_omega v.s D/R^2
figure(23);clf;hold on
title('D_{\omega} v.s. D_0')
plot(theta_0_matrix,D_omega_matrix./(D./R_matrix.^2))
legend('D_{\omega}/(D_0/R^2)','Location','northwest')
xlabel('\omega_0 \delta t')
axis([-inf inf 0 inf])


%% Plotting D_theta v.s 4D/(theta_0^2 R^2)
figure(24);clf;hold on
title('D_{\theta} v.s. D_0')
plot(theta_0_matrix,D_theta_matrix./(4*D./(theta_0_matrix.^2.*R_matrix.^2)))
legend('D_{\theta}/(4 D_0/(\theta_0^2 R^2))','Location','northwest')
xlabel('\omega_0 \delta t')
axis([-inf inf 0 inf])

%% Saving mat file 
% save_file_name=['Date=',Date,', delta_t=',num2str(min(delta_t_matrix)),':',num2str(max(delta_t_matrix)),', v_0=',num2str(min(v_0_matrix)),':',num2str(max(v_0_matrix)),', dt=',num2str(dt),', time=',num2str(Obs_time_steps),'.mat']
save(save_file_name)

