clear;
%% Parameters for importing mat files.

Temp=0.01
Delta_t_matrix=[ 0 0.3 0.6 0.9 1.05 1.2 1.5 1.8]
% T_matrix=[1]
V_0_matrix=2.1 % This value is guessed
% dt=10^-2
file_name_matrix=[
    'dt=0.000.tdms_angle';
    'dt=0.300.tdms_angle';
    'dt=0.600.tdms_angle';
    'dt=0.900.tdms_angle';
    'dt=1.050.tdms_angle';
    'dt=1.200.tdms_angle';
    'dt=1.500.tdms_angle';
    'dt=1.800.tdms_angle';
    ];
%% k_B
k_B=1.38*10^-23;
%% Execution Parameters
recalculate_theta='no' % normally just set to no
recalculate_hist='no' %
recalculate_R='no' % normally just set to no; R is calculated in hist already.
draw_hist='no'
autocorrelation='no'

recalculate_T_eff='yes'
    plot_hist_fit_T_eff='yes';
fit_Viktor_method='yes'
notation='theta'
%%
num_transitions_matrix=[];
theta_plus_matrix=[];
theta_minus_matrix=[];
R_matrix=[];

sigma_matrix=[];
mu_matrix=[];
T_eff_matrix=[];
D_eff_matrix=[];
time_duration_matrix=[];
%%
for file_name_index=1:size(file_name_matrix,1)
    file_name=file_name_matrix(file_name_index,:);
close all
%% Input File Name
% movie_name=['test3']
% movie_name=['2020.11.25,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% movie_name=[Date,',dt=10e-3 take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
movie_name=file_name
[movie_name,'.mat'];
% load([movie_name,'.mat'],'num_transitions','theta_plus','theta_minus')
load([movie_name,'.mat'])

%% Calculating theta
switch recalculate_theta
    case 'yes'
        Analyze_theta=tic;
        moving_avg=1;
        plot_rot='no';
        Theta_Analysis_Fixed_Center(movie_name,partition_movie,N,v_0,Obs_time_steps,partition_time_steps,delta_t,dt,moving_avg,plot_rot);
        time_analyze_theta=toc(Analyze_theta)
    case 'no'
        load([movie_name,'.mat'],'theta')
%         nth_take
end
%% Plotting histogram
switch recalculate_hist
    case 'yes'
        moving_avg=1 ;
        num_bins=100  ;
        bin_limit=2;
        f=figure(81),clf % Note that hist_analysis only works for one particle orbitting a fixed particle at the moment
        switch draw_hist
            case 'no'
                f.Visible='off';
        end
        [k_trans,theta_plus,theta_minus,num_transitions]=hist_analysis(movie_name,moving_avg,num_bins,bin_limit,Obs_time_steps,delta_t,dt);
        set(gca, 'YScale', 'linear')
        saveas(gcf,[movie_name,' (hist).png'])
        
    case 'no'
        load([movie_name,'.mat'],'num_transitions','theta_plus','theta_minus')
%         nth_take
end
%% Calculating R_mean. If recalculate_hist is 'yes' then recalculate_R can be set to 'no'
switch recalculate_R
    case 'yes'
        %% Calculating the histograms and stuff
        theta=0;
        Analyze_theta=tic;
        moving_avg=1 ;
        plot_rot='no';
        Theta_Analysis_Fixed_Center(movie_name,partition_movie,N,v_0,Obs_time_steps,partition_time_steps,delta_t,dt,moving_avg,plot_rot);
        time_analyze_theta=toc(Analyze_theta)
    case 'no'
        load([movie_name,'.mat'],'R_mean')
end
%% Autocorrelation
switch autocorrelation
    case 'yes'
        close all
        num_lag=20;
        [acf,lags,bounds] =autocorr(diff(theta(2,:)),'NumLags',num_lag);
        plot(0:dt:(num_lag)*dt,acf)
        xline(delta_t,'g')
        xlabel('\tau')
        ylabel('g(\tau)')
        title('Autocorrelation function of \Delta\theta')
        legend('g(\tau)',['\delta t=',num2str(delta_t)])
end
%% Calculating Effective Temperature
% switch recalculate_T_eff
%     case 'yes'
%         if Delta_t_matrix(file_name_index)==0
%             D_eff=0;T_eff=0;sigma=0;mu=0;
%         else
%         %% Finding the sigma of the omega (effective temperature)
%         f=figure(80);clf;
%         switch plot_hist_fit_T_eff
%             case 'no'
%                 f.Visible='off';
%         end
%         h=histogram(diff(theta(2,:))/Delta_t_matrix(file_name_index));
%         y=h.Values;
%         x=h.BinEdges(1:end-1)+0.5*h.BinWidth;
%         % plot(x,y);
%         [fitresult,gof]=fit_Gaussian(x,y,plot_hist_fit_T_eff);
%         % a1*exp(-((x-b1)/c1)^2)
%         mu=fitresult.b1;
%         sigma=fitresult.c1/sqrt(2);
%         
%         %% Calculating the D_eff and T_eff
%         D_eff=sigma^2/(2*dt);
%         T_eff=(V_0_matrix/(R_mean(2))*Delta_t_matrix(file_name_index)^2*sigma^2)/(4*k_B*dt);
%         end
%         save([movie_name,'.mat'],'D_eff','T_eff','sigma','mu','-append')
%     case 'no'
%         load([movie_name,'.mat'],'D_eff','T_eff','sigma','mu')
% end
%% (2021.1.7 New method)
switch recalculate_T_eff
    case 'yes'
        %% Finding the sigma of the omega (effective temperature) 
        f=figure(80);clf;
        switch plot_hist_fit_T_eff
            case 'no'
                f.Visible='off';
        end
        switch notation
            case 'omega'
                h=histogram(diff(theta(2,:))/delta_t);  % This is for measuring sigma with histogram(diff(omega))
            case 'theta'
                h=histogram(diff(theta(2,:))); % This is for measuring sigma with histogram(diff(theta))
        end
        Values=h.Values;
        Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
        % plot(x,y);
        [fitresult,gof]=fit_Gaussian(Bins,Values,plot_hist_fit_T_eff);
        % a1*exp(-((x-b1)/c1)^2)
        mu=fitresult.b1;
        sigma=fitresult.c1/sqrt(2);
        
        %% Calculating the D_eff and T_eff
        D_eff=sigma^2/(2*dt);
        T_eff=D_eff;
%                 T_eff=(v_0/R_mean(2))*delta_t^2*sigma^2/(4*k_B*dt); % This is for measuring sigma with histogram(diff(theta)/delta_t)=histogram(diff(omega))
%         D_eff=sigma^2*(v_0/R_mean(2)*delta_t)^2/(8*k_B*dt); % This is for measuring sigma with histogram(diff(theta))
%         T_eff=0;
        save([movie_name,'.mat'],'D_eff','T_eff','sigma','mu','-append')
    case 'no'
        load([movie_name,'.mat'],'D_eff','T_eff','sigma','mu')
end

%% Calculating time_duration
time_duration=(frame(end)-frame(1))*30*10^-3; % 50 ms/ frame
%% Appending the matrices
num_transitions_matrix=[num_transitions_matrix num_transitions];
theta_plus_matrix=[theta_plus_matrix, theta_plus];
theta_minus_matrix=[theta_minus_matrix, theta_minus];
R_matrix=[R_matrix R_mean(2)];

D_eff_matrix=[D_eff_matrix D_eff];
T_eff_matrix=[T_eff_matrix T_eff];
sigma_matrix=[sigma_matrix sigma];
mu_matrix=[mu_matrix mu];
time_duration_matrix=[time_duration_matrix time_duration];
end
%% Saving results
clear theta time v_omega v_x v_y x y
% save(['plotting transitions v_0_matrix=3.5 0.1 10.mat'])
% save(['plotting transitions delta_t_matrix=0.5 0.5 16.mat'])



%% Plotting and Analyzing
for rrrr=1
if 1==0
% v_0_matrix(1)=[];
    %% Delta_t
    if length(Delta_t_matrix)>1
        close all
        %% Bifurcation Diagram
            %% Bifurcation Diagram (Original)
            figure(1);clf
            hold on
            plot(Delta_t_matrix,theta_plus_matrix,'.')
            plot(Delta_t_matrix,theta_minus_matrix,'.')
            title(['Bifurcation Diagram, v_0= ',num2str(V_0_matrix)])
            xlabel('\delta t')
            ylabel('\theta (rad)')
            Bins=@(delta_t,v_0,a)(v_0*delta_t)./(2*a);
            theta_theory=@(x)sqrt(10-sqrt(1/42./x.^6+120./x-20));
            theta_pl=@(x)sqrt(6./x.*(x-1));
            x_line=0:Delta_t_matrix(end)/1000:Delta_t_matrix(end);
            plot(x_line,+theta_theory(Bins(x_line,V_0_matrix,a)),'k')
            plot(x_line,-theta_theory(Bins(x_line,V_0_matrix,a)),'k')
            plot(x_line,theta_pl(Bins(x_line,V_0_matrix,a)),'g')
            plot(x_line,-theta_pl(Bins(x_line,V_0_matrix,a)),'g')
            legend('$\theta_+$','$\theta_-$','$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','','$\theta=\sqrt{\frac{6}{\omega_0 \delta t}(\omega_0 \delta t-1)}$','interpreter','latex','Location','northwest')
            saveas(gcf,['Bifurcation Diagram, v_0= ',num2str(V_0_matrix),'.png'])
            %% Bifurcation Diagram (Normalized)
            figure(2) ;clf

            hold on
            plot(V_0_matrix*Delta_t_matrix./R_matrix,theta_plus_matrix,'.')
            plot(V_0_matrix*Delta_t_matrix./R_matrix,theta_minus_matrix,'.')
            title(['Bifurcation Diagram, v_0= ',num2str(V_0_matrix)])
            xlabel('v_0*\delta t/(2a)')
            ylabel('\theta (rad)')
            x_line=0:Delta_t_matrix(end)/1000:Delta_t_matrix(end);
            Bins=@(delta_t,v_0,a)(v_0*delta_t)./(2*a);
            theta_theory=@(x)sqrt(10-sqrt(1/42./x.^6+120./x-20));
            theta_pl=@(x)sqrt(6./x.*(x-1));
            plot(Bins(x_line,V_0_matrix,a),+theta_theory(Bins(x_line,V_0_matrix,a)),'k')
            plot(Bins(x_line,V_0_matrix,a),-theta_theory(Bins(x_line,V_0_matrix,a)),'k')
            plot(Bins(x_line,V_0_matrix,a),theta_pl(Bins(x_line,V_0_matrix,a)),'g')
            plot(Bins(x_line,V_0_matrix,a),-theta_pl(Bins(x_line,V_0_matrix,a)),'g')
            legend('$\theta_+$','$\theta_-$','$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','','$\theta=\sqrt{\frac{6}{\omega_0 \delta t}(\omega_0 \delta t-1)}$','interpreter','latex','Location','northwest')
            saveas(gcf,['Bifurcation Diagram (Dimensionless), v_0= ',num2str(V_0_matrix),'.png'])

        %% Transition Rates (Original)
            % Plotting original Data
            interest_flag=Delta_t_matrix>1
            delta_t_interest=Delta_t_matrix(interest_flag);
            num_transitions_interest=num_transitions_matrix(interest_flag);
            
            
            figure(3);clf
            hold on
            plot(delta_t_interest,num_transitions_interest,'.')
            plot(Delta_t_matrix(logical(1-interest_flag)),num_transitions_matrix(logical(1-interest_flag)),'x')
            xlabel('\delta t')
            ylabel('Number of Transitions')
            saveas(gcf,['Original, v_0= ',num2str(V_0_matrix),'.png'])
            
            % Fitting with MATLAB tool
            figure(4);clf
            [fitresult, gof] = fit_trans_delta_t(delta_t_interest, num_transitions_interest)
            xlabel('\delta t')
            ylabel('Number of Transitions')
            fitresult.b
            hold on
            plot(Delta_t_matrix(logical(1-interest_flag)),num_transitions_matrix(logical(1-interest_flag)),'x')
            saveas(gcf,['Original (fitted), v_0= ',num2str(V_0_matrix),'.png'])
            %% Transition Rates (Normalized)
            % Transition Rates with normalized variables
            figure(5);clf
            hold on
            R_interest=R_matrix(interest_flag);
            omega_0_interest=delta_t_interest./(R_matrix);
            %             omega_0_interest=delta_t_interest./(R_interest);
            plot(V_0_matrix*omega_0_interest,num_transitions_interest,'.')
            plot(V_0_matrix*Delta_t_matrix(logical(1-interest_flag))./R_matrix(logical(1-interest_flag)),num_transitions_matrix(logical(1-interest_flag)),'x')
            xlabel('v_0 \delta t/R')
            ylabel('Number of Transitions')
            saveas(gcf,['Dimensionless, v_0= ',num2str(V_0_matrix),'.png'])

            % Fitting with MATLAB tool with normalized variables
            figure(6);clf
            omega_delta_t_interest=omega_0_interest*V_0_matrix; % x variable
            num_transitions_interest; % y variable
            [fitresult_norm, gof_norm] = fit_trans_delta_t_norm(omega_delta_t_interest, num_transitions_interest)
            xlabel('v_0 \delta t/R')
            ylabel('Number of Transitions')
            hold on
            plot(V_0_matrix*Delta_t_matrix(logical(1-interest_flag))./R_matrix(logical(1-interest_flag)),num_transitions_matrix(logical(1-interest_flag)),'x')
            legend('data (valid) for theory)','fitted curve','data (invalid for theory)')
            saveas(gcf,['Dimensionless (fitted), v_0= ',num2str(V_0_matrix),'.png'])
            
    end


    %% V_0
    if length(V_0_matrix)>1
        %% Bifurcation Diagram
            % Original
            figure(1) ;clf

            hold on
            plot(V_0_matrix,theta_plus_matrix,'.')
            plot(V_0_matrix,theta_minus_matrix,'.')
            title(['Bifurcation Diagram, \delta t= ',num2str(delta_t),', T=',num2str(Temp)])
            xlabel('v_0')
            ylabel('\theta (rad)')
            Bins=@(delta_t,v_0,a)  (v_0*delta_t)./(2*a);
            theta_theory=@(x)sqrt(10-sqrt(1/42./x.^6+120./x-20));
            theta_pl=@(x)sqrt(6./x.*(x-1));
            x_line=0:V_0_matrix(end)/1000:V_0_matrix(end);
            plot(x_line,+theta_theory(Bins(delta_t,x_line,a)),'k')
            plot(x_line,-theta_theory(Bins(delta_t,x_line,a)),'k')
            plot(x_line,theta_pl(Bins(delta_t,x_line,a)),'g')
            plot(x_line,-theta_pl(Bins(delta_t,x_line,a)),'g')
            legend('$\theta_+$','$\theta_-$','$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','','$\theta=\sqrt{\frac{6}{\omega_0 \delta t}(\omega_0 \delta t-1)}$','interpreter','latex','Location','northwest')
         
            
            % Normalized 
            figure(2) ;clf

            hold on
            plot(V_0_matrix*delta_t./R_matrix,theta_plus_matrix,'.')
            plot(V_0_matrix*delta_t./R_matrix,theta_minus_matrix,'.')
            title(['Bifurcation Diagram, \delta t= ',num2str(delta_t),', T=',num2str(Temp)])
            xlabel('v_0*\delta t/(2a)')
            ylabel('\theta (rad)')
            Bins=@(delta_t,v_0,a)(v_0*delta_t)./(2*a);
            theta_theory=@(x)sqrt(10-sqrt(1/42./x.^6+120./x-20));
            theta_pl=@(x)sqrt(6./x.*(x-1))
            plot(Bins(delta_t,x_line,a),theta_theory(Bins(delta_t,x_line,a)),'k')
            plot(Bins(delta_t,x_line,a),-theta_theory(Bins(delta_t,x_line,a)),'k')
            plot(Bins(delta_t,x_line,a),theta_pl(Bins(delta_t,x_line,a)),'g')
            plot(Bins(delta_t,x_line,a),-theta_pl(Bins(delta_t,x_line,a)),'g')
            legend('$\theta_+$','$\theta_-$','$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','','$\theta=\sqrt{\frac{6}{\omega_0 \delta t}(\omega_0 \delta t-1)}$','interpreter','latex','Location','northwest')

        %% Transition Rates

        % Plotting original Data
        interest_flag=V_0_matrix>5.5
        v_0_interest=V_0_matrix(interest_flag);
        num_transitions_interest=num_transitions_matrix(interest_flag);
        
        figure(3);clf
        hold on
        plot(v_0_interest,num_transitions_interest,'.')
        plot(V_0_matrix(logical(1-interest_flag)),num_transitions_matrix(logical(1-interest_flag)),'x')
        xlabel('v_0')
        ylabel('Number of Transitions')
        
        % Fitting with MATLAB tool
        figure(4)
        [fitresult, gof] = fit_trans_v_0(v_0_interest, num_transitions_interest)
        xlabel('v_0')
        ylabel('Number of Transitions')
        fitresult.b
        
        % Transition Rates with normalized variables
        figure(5);clf
        hold on
        R_interest=R_matrix(interest_flag);
%         omega_0_interest=v_0_interest./(R_interest);
        omega_0_interest=v_0_interest./(2*a);
        plot(delta_t*omega_0_interest,num_transitions_interest,'.')
        plot(delta_t*V_0_matrix(logical(1-interest_flag))./R_matrix(logical(1-interest_flag)),num_transitions_matrix(logical(1-interest_flag)),'x')
        xlabel('v_0 \delta t/R')
        ylabel('Number of Transitions')
        
        % Fitting with MATLAB tool with normalized variables
        figure(6);clf
        omega_delta_t_interest=omega_0_interest*delta_t; % x variable
        num_transitions_interest; % y variable
        [fitresult_norm, gof_norm] = fit_trans_v_0_norm(omega_delta_t_interest, num_transitions_interest)
        xlabel('v_0 \delta t/R')
        ylabel('Number of Transitions')
        fitresult_norm.b/5
        

    end
end
end

%% Save Variables for Viktor
% save('For_Viktor.mat','a','D','delta_t','dt','k_B','num_transitions_matrix','T','R_matrix','v_0_matrix')
%% Transition Rates
% transition_rate_theory=sqrt(2)./(pi.*omega_0_matrix.*delta_t_matrix.^2).*(omega_0_matrix.*delta_t_matrix-1).*exp(-3/2*(omega_0_matrix.*delta_t_matrix-1).^2./(D_eff_matrix.*omega_0_matrix.*delta_t_matrix.^2/2.*omega_0_matrix.*delta_t_matrix.^3));
% omega_0_matrix=V_0_matrix./R_matrix;
% transition_rate_theory=2*sqrt(2)./(pi.*omega_0_matrix.*Delta_t_matrix.^2).*(omega_0_matrix.*Delta_t_matrix-1).*exp(-3/2*(omega_0_matrix.*Delta_t_matrix-1).^2./(omega_0_matrix.*k_B.*T_eff_matrix.*Delta_t_matrix.^3));
%% D_eff and T_eff

D_eff_matrix=sigma_matrix.^2/(2*dt);
T_eff_matrix=D_eff_matrix;
%% Transition Rates
omega_0_matrix=V_0_matrix./R_matrix;
switch notation
    case 'omega'
        
        transition_rate_theory=2*sqrt(2)./(pi*omega_0_matrix.*Delta_t_matrix.^2).*(omega_0_matrix.*Delta_t_matrix-1).*exp(-3*(omega_0_matrix.*Delta_t_matrix-1).^2./(D_eff_matrix.*omega_0_matrix.^2.*Delta_t_matrix.^5));
% transition_rate_theory=2*sqrt(2)./(pi.*omega_0_matrix.*Delta_t_matrix.^2).*(omega_0_matrix.*Delta_t_matrix-1).*exp(-3/2*(omega_0_matrix.*Delta_t_matrix-1).^2./(omega_0_matrix.*k_B.*T_eff_matrix.*Delta_t_matrix.^3));
    case 'theta'
        theta_0_matrix=omega_0_matrix.*Delta_t_matrix;
        transition_rate_theory=2*sqrt(2)./(pi*theta_0_matrix.*Delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(D_eff_matrix.*Delta_t_matrix.*theta_0_matrix.^2));
end
%% 2020.12.10 這邊繼續看要怎麼做
% time_duration=Obs_time_steps*dt+delta_t;
x0=omega_0_matrix.*Delta_t_matrix;
figure(7);clf;hold on;
plot(x0,transition_rate_theory,'o')
plot(x0,num_transitions_matrix./time_duration_matrix,'x')
xlabel('\omega_0 \delta t')
ylabel('Transition Rates (1/s)')
legend('theory','experiment','Location','southwest')

figure(8);clf;hold on;
flag=(x0>1.);
plot(x0(flag),transition_rate_theory(flag),'o')
plot(x0(flag),num_transitions_matrix(flag)./time_duration_matrix(flag),'x')
xlabel('\omega_0 \delta t')
ylabel('Transition Rates (1/s)')
legend('theory','experiment','Location','southeast')
%% Plotting R vs omega_0 delta_t
figure(10);clf;
plot(Delta_t_matrix,R_matrix/(2*a))
% xlabel('\omega_0 \delta t')
xlabel('\delta t')
ylabel('R/(2a)')
axis([-inf inf 0 inf])

%% Plotting inverse transition rate
% time_duration=Obs_time_steps*dt+delta_t;
% omega_0_matrix=v_0_matrix./(2*a);
% omega_0_matrix=v_0_matrix./R_matrix;
% transition_rate_theory=sqrt(2)./(pi.*omega_0_matrix.*delta_t_matrix.^2).*(omega_0_matrix.*delta_t_matrix-1).*exp(-3/2*(omega_0_matrix.*delta_t_matrix-1).^2./(omega_0_matrix.*k_B.*T_eff_matrix.*delta_t_matrix.^3));

% x0=omega_0_matrix.*delta_t_matrix;

figure(11);clf;hold on;
plot(x0,1./transition_rate_theory,'o')
plot(x0,1./(num_transitions_matrix./time_duration_matrix),'x')
xlabel('\omega_0 \delta t')
ylabel('1/Transition Rates (s)')
legend('theory','experiment','Location','southeast')


figure(12);clf;hold on;
flag=(x0>1 );
plot(x0(flag),1./transition_rate_theory(flag),'o')
plot(x0(flag),1./(num_transitions_matrix(flag)./time_duration_matrix(flag)),'x')
xlabel('\omega_0 \delta t')
ylabel('1/Transition Rates (s)')
legend('theory','experiment','Location','northeast')

%% Plotting the Coefficients

%% Plotting T_eff and D_eff
figure(15);clf;
plot(x0,T_eff_matrix)
xlabel('\omega_0 \delta_t')
ylabel('T_{eff}')
figure(16);clf;
plot(x0,D_eff_matrix)
xlabel('\omega_0 \delta_t')
ylabel('D_{eff}')

%% Plotting T_eff and D_eff & R vs omega_0 * delta_t
figure(21);clf;hold on
% title('fixed \delta t')
% title('fixed \omega_0')
plot(x0,T_eff_matrix)
plot(x0,D_eff_matrix)
plot(x0,R_matrix/(2*a)/300)
xlabel('\omega_0 \delta t')
legend('T_{eff}','D_{eff}','R/(2a)/300','Location','northwest')

%% Plotting D_eff v.s D/R^2
    %% We don't know D for sure, set to be 1 to show 
D=1
figure(23);clf;hold on
title('D_{eff} v.s. D_0')
plot(x0,D_eff_matrix./(D./R_matrix.^2))
legend('D_{eff}/(D/R^2)','Location','northwest')
xlabel('\omega_0 \delta t')
axis([-inf inf 0 inf])





