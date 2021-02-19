%% Extension of simulation_1D_matrix, with all other kinds of D_eff and predicted trans_rate calculated, meant to compare with Viktor's codes at first
%% Updated on 2021.1.22 Basically the simulation are the same as Viktors; the transition rates calculated in find_theta_plus have been modified 
%% and now the transition rates are the same as Viktor's (check again if have time)

%% Updated on 2021.1.25: Now the theoretical predictions and the simulation datas are really close to Viktor's. This code here basically can be discarded.
clear;
close all



%% This setup starts from theta_0>1
nth_take=1 % for a=5, particle 1 fixed
delta_t_matrix=1
T_matrix=[1]

theta_0_matrix=linspace(1.1,1.6,15);
D_0=20
dt=10^-2
intrinsic_delay=0.0 % Intrinsic delay
Obs_time_steps=50000/dt
R_mean=10;
v_0_matrix=theta_0_matrix*R_mean/delta_t_matrix;

%%
plot_hist_fit_T_eff='yes';
%% Matrices for analysis

theta_plus_matrix=[];
theta_minus_matrix=[];
R_matrix=[];
% D_omega_matrix=[];
D_theta_full_matrix=[];
D_theta_approx_matrix=[];
theta_0_matrix=[];
num_transitions_full_matrix=[];
num_transitions_approx_matrix=[];
D_omega_full_matrix=[];
D_omega_full2_matrix=[];
%%
for delta_t_index=1:length(delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(v_0_matrix)
            if 1
                delta_t=delta_t_matrix(delta_t_index);
                v_0= v_0_matrix(v_0_index);

                
                D=D_0/R_mean^2;
                omega_0=v_0/R_mean;
%                 omega_0*delta_t
                % omega_0=1;
                % delta_t=100;
                
                %% Simulation setup
                theta_approx(1:1+Obs_time_steps)=0; %% theta=theta_delay is omega*delta_t
                phi(1:1+Obs_time_steps+round(delta_t/dt))=0; %% theta_polar is the ever-increasing polar angle
                theta_0=omega_0*delta_t;
                theta_plus_approx_theory=sqrt(6*(theta_0-1)/theta_0);
                r_theta_sim=normrnd(0,sqrt(8*D*dt/theta_0^2),[Obs_time_steps,1]);
                r_phi=normrnd(0,sqrt(2*D*dt),[Obs_time_steps,1]);
                %% For approximate result (We deal with this later)
                for k=1:Obs_time_steps
                    theta_approx(k+1)=theta_approx(k)+1/(3*delta_t)*(theta_plus_approx_theory^2-theta_approx(k)^2)*theta_approx(k)*dt+r_theta_sim(k);
                end
                figure;plot(theta_approx)
                %% For full simulation (This is what I'm more interested in, deal with this first)
                for k=1:Obs_time_steps
                    phi(1+k+round(delta_t/dt))=phi(k+round(delta_t/dt))+omega_0*sin(phi(k+round(delta_t/dt))-phi(k))*dt+r_phi(k);
                end
                
                %% Analyzing phi
                theta_full=phi(1+round(delta_t/dt):1+Obs_time_steps+round(delta_t/dt))-phi(1:1+Obs_time_steps);
                omega_full=diff(phi)/dt;
                %% Replacing theta with theta_sim
%                 clear theta;
%                 theta_approx;
                %%
                
                [V_parameters,theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta_full,theta_0);
                [theta_approx_stable,k_trans_theta_approx,theta_plus_approx,theta_minus_approx,num_transitions_theta_approx]=find_theta_plus(theta_approx,theta_0);
                %% Finding the sigma of the theta (full 1D)
                % [theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0)
                f=figure(80);clf;
%                 plot_hist_fit_T_eff='no';
                switch plot_hist_fit_T_eff
                    case 'no'
                        f.Visible='off';
                end
                h_theta_full=histogram(diff(theta_full));
                Values_theta_full=h_theta_full.Values;
                Bins_theta_full=h_theta_full.BinEdges(1:end-1)+0.5*h_theta_full.BinWidth;
                
%                 theta_0=v_0*delta_t/R_mean;
                
                [fitresult,gof]=fit_Gaussian(Bins_theta_full,Values_theta_full,plot_hist_fit_T_eff);
                mu_theta=fitresult.b1;
                sigma_theta=fitresult.c1/sqrt(2);
                %% Finding the sigma of the theta (approx 1D)
                % [theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0)
                f=figure(81);clf;
%                 plot_hist_fit_T_eff='no';
                switch plot_hist_fit_T_eff
                    case 'no'
                        f.Visible='off';
                end
                h_theta_approx=histogram(diff(theta_approx));
                Values_theta_approx=h_theta_approx.Values;
                Bins_theta_approx=h_theta_approx.BinEdges(1:end-1)+0.5*h_theta_approx.BinWidth;
                
%                 theta_0=v_0*delta_t/R_mean;
                
                [fitresult_approx,gof_approx]=fit_Gaussian(Bins_theta_approx,Values_theta_approx,plot_hist_fit_T_eff);
                mu_theta_approx=fitresult_approx.b1;
                sigma_theta_approx=fitresult_approx.c1/sqrt(2);
                %% 2021.1.21 test omega's fitting
                close all
%                 histogram(diff(theta_full)/dt)
                
                f=figure(82);clf;
%                 plot_hist_fit_T_eff='no';
                switch plot_hist_fit_T_eff
                    case 'no'
                        f.Visible='off';
                end
                
                h_omega_full=histogram(omega_full);
                Values_omega_full=h_omega_full.Values;
                Bins_omega_full=h_omega_full.BinEdges(1:end-1)+0.5*h_omega_full.BinWidth;
                
%                 theta_0=v_0*delta_t/R_mean;
                
                [fitresult_omega_full,gof_omega_full]=fit_Gaussian(Bins_omega_full,Values_omega_full,plot_hist_fit_T_eff);
                mu_omega_full=fitresult_omega_full.b1;
                sigma_omega_full=fitresult_omega_full.c1/sqrt(2);
                
                % Viktor's method 
                Values_omega_full = Values_omega_full/trapz(Bins_omega_full,Values_omega_full);
                sigma2 = trapz(Bins_omega_full,Bins_omega_full.^2.*Values_omega_full) - (trapz(Bins_omega_full,Bins_omega_full.*Values_omega_full))^2;
                sigma_omega_full=sqrt(sigma2);
                %% 2021.1.21 Using the equation for omega for testing.
                f=figure(83);clf;
%                 plot_hist_fit_T_eff='no';
                switch plot_hist_fit_T_eff
                    case 'no'
                        f.Visible='off';
                end
                
                h_omega_full2=histogram(diff(omega_full));
                Values_omega_full2=h_omega_full2.Values;
                Bins_omega_full2=h_omega_full2.BinEdges(1:end-1)+0.5*h_omega_full2.BinWidth;
                
%                 theta_0=v_0*delta_t/R_mean;
                
                [fitresult_omega_full2,gof_omega_full2]=fit_Gaussian(Bins_omega_full2,Values_omega_full2,plot_hist_fit_T_eff);
                mu_omega_full2=fitresult_omega_full2.b1;
                sigma_omega_full2=fitresult_omega_full2.c1/sqrt(2);
                

                %% 2021.1.21 Using the equation for time derivative of omega=omega_0*sin(...)+noise
                f=figure(84);clf;
%                 plot_hist_fit_T_eff='no';
                switch plot_hist_fit_T_eff
                    case 'no'
                        f.Visible='off';
                end
                omega=diff(phi)/dt;
                domega_dt=diff(omega)/dt;
                omega_full_3=domega_dt(1+delta_t/dt:end)-omega_0*cos(phi(1+delta_t/dt:end-2)-phi(1:end-delta_t/dt-2)).*sin(omega(1+delta_t/dt:end-1)-omega(1:end-delta_t/dt-1));
                
                h_omega_full3=histogram(diff(omega_full_3));
                Values_omega_full3=h_omega_full3.Values;
                Bins_omega_full3=h_omega_full3.BinEdges(1:end-1)+0.5*h_omega_full3.BinWidth;
                
%                 theta_0=v_0*delta_t/R_mean;
                
                [fitresult_omega_full3,gof_omega_full2]=fit_Gaussian(Bins_omega_full3,Values_omega_full3,plot_hist_fit_T_eff);
                mu_omega_full3=fitresult_omega_full3.b1;
                sigma_omega_full3=fitresult_omega_full3.c1/sqrt(2);
                
                %% Finding D_theta
                D_theta_full=sigma_theta^2/(2*dt);
                D_theta_approx=sigma_theta_approx^2/(2*dt);
                D_omega_full=sigma_omega_full^2*dt/2;
                D_omega_full2=sigma_omega_full2^2/(2*dt);
%                 D_omega_full3=sigma_omega_full2^3/(2*dt);
                %% Appending Matrices
%                 theta_0
%                 num_transitions_theta
%                 figure(82);plot(theta_approx);
                num_transitions_full_matrix=[num_transitions_full_matrix num_transitions_theta];
                num_transitions_approx_matrix=[num_transitions_approx_matrix num_transitions_theta_approx];
                theta_plus_matrix=[theta_plus_matrix theta_plus];
                theta_minus_matrix=[theta_minus_matrix theta_minus];
                R_matrix=[R_matrix R_mean];
                D_theta_full_matrix=[D_theta_full_matrix D_theta_full];
                D_theta_approx_matrix=[D_theta_approx_matrix D_theta_approx];
                D_omega_full_matrix=[D_omega_full_matrix D_omega_full];
                D_omega_full2_matrix=[D_omega_full2_matrix D_omega_full2];
                theta_0_matrix=[theta_0_matrix theta_0];

            end
            nth_take=nth_take+1
        end
        nth_take=nth_take+1
    end
    nth_take=nth_take+1
end
%%
close all
clear domega_dt omega omega_full omega_full_3 phi r_phi r_theta_sim theta_approx theta_approx_stable theta_full theta_stable
%% Start analyzing the data:
time_duration=Obs_time_steps.*dt+delta_t_matrix;
transition_rate_approx=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(4*D_0./R_matrix.^2.*delta_t_matrix));
% transition_rate_full=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(D_theta_full_matrix.*delta_t_matrix.*theta_0_matrix.^2));
transition_rate_full=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(2*D_0./R_matrix.^2.*delta_t_matrix.*theta_0_matrix.^2));



tranistion_rate_approx_double=2*transition_rate_approx;
transition_rate_full_double=2*transition_rate_full;
%% 2021.1.20 Something wrong with the transition rates? Delete this afterwards.
% close all;
% figure;hold on
% trans_rate_matrix_full=num_transitions_full_matrix./time_duration;
% D_eff_matrix=4*D./theta_0_matrix.^2;
% % D_eff_matrix(1:length(theta_0_matrix))=2*D;
% transition_rate_theory=2*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./((D_eff_matrix.*theta_0_matrix.^2.*delta_t_matrix)));
% 
% plot(theta_0_matrix,trans_rate_matrix_full,'xr')
% plot(theta_0_matrix,transition_rate_theory)
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
%     
% 
% figure;hold on
% plot(theta_0_matrix,1./trans_rate_matrix_full,'xr')
% plot(theta_0_matrix,1./transition_rate_theory)
%% 2021.1.18 Plotting Transition Rates of uncorrected and corrected together
    trans_rate_matrix_full=num_transitions_full_matrix./time_duration;
    %     trans_rate_matrix_approx=num_transitions_approx_matrix./time_duration;
    trans_rate_matrix_approx=zeros(size(trans_rate_matrix_full));
    figure(7);clf;hold on;
    plot(theta_0_matrix,transition_rate_approx,'-b')
    plot(theta_0_matrix,transition_rate_full,'-r')
%     plot(theta_0_matrix,transition_rate_today,'-g')
    plot(theta_0_matrix,trans_rate_matrix_full,'xr')
    plot(theta_0_matrix,trans_rate_matrix_approx,'ob')
    plot(theta_0_matrix,tranistion_rate_approx_double,'b-')
    plot(theta_0_matrix,transition_rate_full_double,'-r')
%     plot(theta_0_matrix,transition_rate_today_double,'-g')
    title('Full range of \theta_0')
    xlabel('\omega_0 \delta t')
    ylabel('Transition Rates (1/s)')
%     legend('approximated 1D (D_{eff}=4D_0/(\theta_0^2 R^2))','today','full (D_{eff}=2D_0/R^2)','full 1D simulation','approximated 1D simulation')
    legend('approximated 1D (D_{eff}=4D_0/(\theta_0^2 R^2))','full 1D (D_{eff}=2D_0/R^2)','full 1D simulation','approximated 1D simulation')

    
    figure(8);clf;hold on;
    flag=(theta_0_matrix>1.);
    plot(theta_0_matrix(flag),transition_rate_approx(flag),'b-')
    plot(theta_0_matrix(flag),transition_rate_full(flag),'r-')
%     plot(theta_0_matrix(flag),transition_rate_today(flag),'g-')
    plot(theta_0_matrix(flag),trans_rate_matrix_full(flag),'xr')
    plot(theta_0_matrix(flag),trans_rate_matrix_approx(flag),'ob')
    plot(theta_0_matrix(flag),tranistion_rate_approx_double(flag),'b-')
    plot(theta_0_matrix(flag),transition_rate_full_double(flag),'r-')
%     plot(theta_0_matrix(flag),transition_rate_today_double(flag),'g-')
    title('\theta_0>1')
    xlabel('\omega_0 \delta t')
    ylabel('Transition Rates (1/s)')
%     legend('uncorrected','corrected','today','full 1D simulation','approximated 1D simulation')
legend('approximated 1D (D_{eff}=4D_0/(\theta_0^2 R^2))','full 1D (D_{eff}=2D_0/R^2)','full 1D simulation','approximated 1D simulation')

    set(gca,'XScale','log')
    set(gca,'YScale','log')
%     axis([1.1 1.43 10^-7 10^-1])
    %% Plotting inverse transition rate of uncorrected and corrected together
    figure(11);clf;hold on;
    plot(theta_0_matrix,1./(transition_rate_approx),'b-')
    plot(theta_0_matrix,1./transition_rate_full,'r-')
%     plot(theta_0_matrix,1./transition_rate_today,'g-')
    plot(theta_0_matrix,1./(trans_rate_matrix_full),'rx')
    plot(theta_0_matrix,1./(trans_rate_matrix_approx),'bo')
    plot(theta_0_matrix,1./(tranistion_rate_approx_double),'b-')
    plot(theta_0_matrix,1./transition_rate_full_double,'r-')
%     plot(theta_0_matrix,1./transition_rate_today_double,'-g')
    title('Full range of \theta_0')
    xlabel('\omega_0 \delta t')
    ylabel('1/Transition Rates (s)')
%     legend('uncorrected','corrected','today','full 1D simulation','approximated 1D simulation')
legend('approximated 1D (D_{eff}=4D_0/(\theta_0^2 R^2))','full 1D (D_{eff}=2D_0/R^2)','full 1D simulation','approximated 1D simulation')

    %     set(gca,'XScale','log')
    %     set(gca,'YScale','log')
    
    figure(12);clf;hold on;
    flag=(theta_0_matrix>1 );
    plot(theta_0_matrix(flag),1./transition_rate_approx(flag),'-b')
    plot(theta_0_matrix(flag),1./transition_rate_full(flag),'-r')
%     plot(theta_0_matrix(flag),1./transition_rate_today(flag),'-g')
    plot(theta_0_matrix(flag),1./(trans_rate_matrix_full(flag)),'xr')
    plot(theta_0_matrix(flag),1./(trans_rate_matrix_approx(flag)),'bo')
    plot(theta_0_matrix(flag),1./tranistion_rate_approx_double(flag),'b-')
    plot(theta_0_matrix(flag),1./transition_rate_full_double(flag),'-r')
%     plot(theta_0_matrix(flag),1./transition_rate_today_double(flag),'g-')
    title('\theta_0>1')
    xlabel('\omega_0 \delta t')
    ylabel('1/Transition Rates (s)')
%     legend('uncorrected','corrected','today','full 1D simulation','approximated 1D simulation')
legend('approximated 1D (D_{eff}=4D_0/(\theta_0^2 R^2))','full 1D (D_{eff}=2D_0/R^2)','full 1D simulation','approximated 1D simulation')

    set(gca,'XScale','log')
    set(gca,'YScale','log')
    axis([1.1 1.43 10^1 10^7])
    
    %%
    close(figure(7))
    close(figure(11))
    figure(8)
    figure(12)
    
    %% Plotting D_eff with 2D/R^2 and 4D/(R^2 theta_0^2)
    figure(13);clf;hold on;
    plot(theta_0_matrix,D_theta_full_matrix,'rx')
    plot(theta_0_matrix,D_theta_approx_matrix,'bo')
    plot(theta_0_matrix,4*D./theta_0_matrix.^2,'-b')
    yline(2*D,'-r')
    axis([-inf inf 0 inf])
    xlabel('\omega_0 \delta t')
    legend('full 1D', 'approx 1D', '4D_0/(\theta_0^2 R^2))','2D_0/R^2')
  
    %% Plotting D_eff gained from omega_full
    close all
    figure;clf ;hold on
   
    plot(theta_0_matrix,D_omega_full_matrix)
    plot(theta_0_matrix,D_theta_full_matrix)
    
    figure;clf;hold on
    plot(theta_0_matrix, D_omega_full_matrix-D_theta_full_matrix)
    
    %% Plotting D_eff from omega_full
    close all
    clf;
    plot(theta_0_matrix,D_omega_full_matrix)
    
    %% Plotting num_transitions_matrix
    figure;clf;
    plot(theta_0_matrix,num_transitions_full_matrix)

    %% Plotting D_ef from omega_full2
    %% No I don't think this works because then D_eff depends on dt. This still has to blame that dot(omega) itself doesn't really go under D_eff=4D/(omega_0^2*delta_t^4)
    figure;clf;
   
    plot(theta_0_matrix,D_omega_full2_matrix)
     axis([-inf,inf,0,inf])
    transition_rate_omega2_measured=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(4*D_0./R_matrix.^2.*delta_t_matrix));
    transition_rate_omega2_theory=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(4*D_0./R_matrix.^2.*delta_t_matrix));
    
    %% Plotting Bifurcation
    figure;clf;hold on;
    plot(theta_0_matrix, theta_plus_matrix)
    plot(theta_0_matrix, theta_minus_matrix)
        
    %% Plotting on checking Viktor's D_eff to be 4D/theta_0^2R^2
    trans_rate_matrix_full=num_transitions_full_matrix./time_duration;
%     transition_rate_theta=2*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./((2*D).*theta_0_matrix.^2.*delta_t_matrix));
    transition_rate_theta=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./((4*D./theta_0_matrix.^2).*theta_0_matrix.^2.*delta_t_matrix));
    
    figure(7);clf;hold on;
    plot(theta_0_matrix,transition_rate_theta,'-r')
    plot(theta_0_matrix,trans_rate_matrix_full,'xr')
%     plot(theta_0_matrix,transition_rate_full_double,'-r')
    title('Full range of \theta_0')
    xlabel('\omega_0 \delta t')
    ylabel('Transition Rates (1/s)')
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
    
    figure(8);clf;hold on;
    flag=(theta_0_matrix>1.);
    plot(theta_0_matrix(flag),transition_rate_theta(flag),'r-')
    plot(theta_0_matrix(flag),trans_rate_matrix_full(flag),'xr')
%     plot(theta_0_matrix(flag),transition_rate_full_double(flag),'r-')
    title('\theta_0>1')
    xlabel('\omega_0 \delta t')
    ylabel('Transition Rates (1/s)')

    set(gca,'XScale','log')
    set(gca,'YScale','log')
    
    %% Plotting inverse

    figure(9);clf;hold on;
    plot(theta_0_matrix,1./transition_rate_theta,'-r')
    plot(theta_0_matrix,1./trans_rate_matrix_full,'xr')
%     plot(theta_0_matrix,transition_rate_full_double,'-r')
    title('Full range of \theta_0')
    xlabel('\omega_0 \delta t')
    ylabel('1/Transition Rates (s)')
%     set(gca,'XScale','log')
%     set(gca,'YScale','log')
    
    figure(10);clf;hold on;
    flag=(theta_0_matrix>1.);
    plot(theta_0_matrix(flag),1./transition_rate_theta(flag),'r-')
    plot(theta_0_matrix(flag),1./trans_rate_matrix_full(flag),'xr')
%     plot(theta_0_matrix(flag),transition_rate_full_double(flag),'r-')
    title('\theta_0>1')
    xlabel('\omega_0 \delta t')
    ylabel('1/Transition Rates (s)')

    set(gca,'XScale','log')
    set(gca,'YScale','log')
    %% Comparing num_transitiions_matrix with Viktor's file; They agree with each other
    close all
    figure;clf;plot(theta_0_matrix,num_transitions_full_matrix)
    figure;clf;plot(theta_0_matrix,num_transitions_full_matrix);
        set(gca,'XScale','log')
    set(gca,'YScale','log')
    %% Saving Matfile
    save('2021.1.21_compare_pc.mat')
    
