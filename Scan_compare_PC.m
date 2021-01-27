%% Same as test.m. Functionalized to compare the transition rate with Viktor's in Compare_1D_2D_pc_Vik.m
function Scan_compare_PC(D_0,dt,delta_t_matrix,v_0_matrix,Obs_time_steps)
close all



Date='2021.1.25'
nth_take=1
% delta_t_matrix=1
T_matrix=[1]

% v_0_matrix=[11:0.5:20]
% dt=10^-2
intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6

%% Running main.m
[num_transitions_matrix,theta_plus_matrix,theta_minus_matrix,R_matrix,D_omega_matrix,D_theta_matrix,theta_0_matrix,R_recip_matrix]= Scan_main(Date,nth_take,delta_t_matrix,T_matrix,v_0_matrix,dt,intrinsic_delay,Obs_time_steps,D_0);
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
% omega_0_matrix=v_0_matrix./R_matrix;
% time_duration=Obs_time_steps.*dt+delta_t_matrix;
% theta_0_matrix=omega_0_matrix.*delta_t_matrix;
%% Comparing E_b and k_B*T
k_B=1;
save(['2021.1.25_compare_PC,D_0=',num2str(D_0),'.mat'])
end
