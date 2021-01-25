%% Same as test.m. Functionalized to compare the transition rate with Viktor's in Compare_1D_2D_pc_Vik.m
function compare_PC(D_0)
close all



Date='2021.1.25'
nth_take=1
delta_t_matrix=1
T_matrix=[1]
% v_0_matrix=[5:0.2:8]*2
v_0_matrix=[11:0.5:20]
dt=10^-3
intrinsic_delay=0.0 % Intrinsic delay
Obs_time_steps=10^6


%%
% D_0=1;
%% File_name
save_file_name=['Date=',Date,', delta_t=',num2str(min(delta_t_matrix)),':',num2str(max(delta_t_matrix)),', v_0=',num2str(min(v_0_matrix)),':',num2str(max(v_0_matrix)),', dt=',num2str(dt),', time=',num2str(Obs_time_steps),'.mat']
% pause
% load(save_file_name)
%% Running main.m
[num_transitions_matrix,theta_plus_matrix,theta_minus_matrix,R_matrix,D_omega_matrix,D_theta_matrix,theta_0_matrix]= main(Date,nth_take,delta_t_matrix,T_matrix,v_0_matrix,dt,intrinsic_delay,Obs_time_steps,D_0);
%% Running parameters
Bifurcation_Diagram='no'

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

notation='theta'
omega_0_matrix=v_0_matrix./R_matrix;
time_duration=Obs_time_steps.*dt+delta_t_matrix;

% this (uncorrected) transition rate is completely calculated with Kramer's theory
tranistion_rate_uncorrected=2*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(4*D_0./R_matrix.^2.*delta_t_matrix));
% this (today,2021.1.19) transition is an alternative rate suggested by Viktor for my
% 2-D simulation, which is very close to corrected rate and was inspired by
% D_theta/(D/R^2)~2
transition_rate_today=2*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(2*D_0./R_matrix.^2.*theta_0_matrix.^2.*delta_t_matrix));

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

  


%% Comparing E_b and k_B*T
k_B=1;


%% Saving mat file 
% save_file_name=['Date=',Date,', delta_t=',num2str(min(delta_t_matrix)),':',num2str(max(delta_t_matrix)),', v_0=',num2str(min(v_0_matrix)),':',num2str(max(v_0_matrix)),', dt=',num2str(dt),', time=',num2str(Obs_time_steps),'.mat']
% save(save_file_name)

%% 2021.1.25 Test again after fixing the transition rates calculation in find_theta_plus.m
%% Code below are copied from simlation_1_D_matrix_compare_with_Viktor
%% Start analyzing the data:
time_duration=Obs_time_steps.*dt+delta_t_matrix;
transition_rate_approx=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(4*D_0./R_matrix.^2.*delta_t_matrix));
% transition_rate_full=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(D_theta_full_matrix.*delta_t_matrix.*theta_0_matrix.^2));
transition_rate_full=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(2*D_0./R_matrix.^2.*delta_t_matrix.*theta_0_matrix.^2));



tranistion_rate_approx_double=2*transition_rate_approx;
transition_rate_full_double=2*transition_rate_full;
%% 2021.1.18 Plotting Transition Rates of uncorrected and corrected together
    trans_rate_matrix_full=num_transitions_matrix./time_duration;
    %     trans_rate_matrix_approx=num_transitions_approx_matrix./time_duration;
    trans_rate_matrix_approx=zeros(size(trans_rate_matrix_full));
    

save(['2021.1.25_compare_PC,D_0=',num2str(D_0),'.mat'])



end
