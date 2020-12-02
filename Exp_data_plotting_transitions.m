clear;
%% Parameters for importing mat files.

% Date='2020.11.25'
% nth_take=100
Temp=0.01
Delta_t_matrix=[ 0 0.3 0.6 0.9 1.05 1.2 1.5 1.8]
% T_matrix=[1]
V_0_matrix=1
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

%% Execution Parameters
recalculate_theta='no' % normally just set to no
recalculate_hist='no' %
recalculate_R='no' % normally just set to no; R is calculated in hist already.
draw_hist='no'
%%
num_transitions_matrix=[];
theta_plus_matrix=[];
theta_minus_matrix=[];
R_matrix=[];
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
%% Appending the matrices
num_transitions_matrix=[num_transitions_matrix num_transitions];
theta_plus_matrix=[theta_plus_matrix, theta_plus];
theta_minus_matrix=[theta_minus_matrix, theta_minus];
R_matrix=[R_matrix R_mean(2)];
end
%% Saving results
clear theta time v_omega v_x v_y x y
% save(['plotting transitions v_0_matrix=3.5 0.1 10.mat'])
% save(['plotting transitions delta_t_matrix=0.5 0.5 16.mat'])



%% Plotting and Analyzing
if 1
% v_0_matrix(1)=[];
    %% Delta_t
    if length(Delta_t_matrix)>1
        close all
        %% Bifurcation Diagram
%         theta_plus_matrix=theta_plus_matrix*pi/180;
%         theta_minus_matrix=theta_minus_matrix*pi/180;
            %% Bifurcation Diagram (Original)
            figure(1);clf
            hold on
            plot(Delta_t_matrix,theta_plus_matrix,'.')
            plot(Delta_t_matrix,theta_minus_matrix,'.')
            title(['Bifurcation Diagram, v_0= ',num2str(V_0_matrix)])
            xlabel('\delta t')
            ylabel('\theta (rad)')
            % figure(2)
            fraction_in_theory=@(delta_t,v_0,a)2*a./(v_0*delta_t);
            theta_theory=@(x)sqrt(10-sqrt(x.^6/42+120*x-20));
            V_0_matrix=2.1
            theta_theory(fraction_in_theory(Delta_t_matrix,V_0_matrix,a))
            x_line=0:Delta_t_matrix(end)/1000:Delta_t_matrix(end);
            plot(x_line,+theta_theory(fraction_in_theory(x_line,V_0_matrix,a)),'k')
            plot(x_line,-theta_theory(fraction_in_theory(x_line,V_0_matrix,a)),'k')
%             plot(Delta_t_matrix,+theta_theory(x(Delta_t_matrix,V_0_matrix,a)),'k')
%             plot(Delta_t_matrix,-theta_theory(x(Delta_t_matrix,V_0_matrix,a)),'k')
            saveas(gcf,['Bifurcation Diagram, v_0= ',num2str(V_0_matrix),'.png'])
            %% Bifurcation Diagram (Normalized)
            figure(2) ;clf

            hold on
            plot(V_0_matrix*Delta_t_matrix/(2*a),theta_plus_matrix,'.')
            plot(V_0_matrix*Delta_t_matrix/(2*a),theta_minus_matrix,'.')
            title(['Bifurcation Diagram, v_0= ',num2str(V_0_matrix)])
            xlabel('v_0*\delta t/(2a)')
            ylabel('\theta (rad)')
            % figure(2)
            fraction_in_theory=@(delta_t,v_0,a)2*a./(v_0*delta_t);
            theta_theory=@(x)sqrt(10-sqrt(x.^6/42+120*x-20));
            theta_theory(fraction_in_theory(Delta_t_matrix,V_0_matrix,a))
            plot(V_0_matrix*x_line/(2*a),+theta_theory(fraction_in_theory(x_line,V_0_matrix,a)),'k')
            plot(V_0_matrix*x_line/(2*a),-theta_theory(fraction_in_theory(x_line,V_0_matrix,a)),'k')
%             plot(V_0_matrix*Delta_t_matrix/(2*a),+theta_theory(x(Delta_t_matrix,V_0_matrix,a)),'k')
%             plot(V_0_matrix*Delta_t_matrix/(2*a),-theta_theory(x(Delta_t_matrix,V_0_matrix,a)),'k')
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
            omega_0_interest=delta_t_interest./(2*a);
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
            fraction_in_theory=@(delta_t,v_0,a)2*a./(v_0*delta_t);
            theta_theory=@(x)sqrt(10-sqrt(x.^6/42+120*x-20));
            theta_theory(fraction_in_theory(delta_t,V_0_matrix,a))
            plot(V_0_matrix,theta_theory(fraction_in_theory(delta_t,V_0_matrix,a)),'k')
            plot(V_0_matrix,-theta_theory(fraction_in_theory(delta_t,V_0_matrix,a)),'k')
            
            
            % Normalized 
            figure(2) ;clf

            hold on
            plot(V_0_matrix*delta_t/(2*a),theta_plus_matrix,'.')
            plot(V_0_matrix*delta_t/(2*a),theta_minus_matrix,'.')
            title(['Bifurcation Diagram, \delta t= ',num2str(delta_t),', T=',num2str(Temp)])
            xlabel('v_0*\delta t/(2a)')
            ylabel('\theta (rad)')
            fraction_in_theory=@(delta_t,v_0,a)2*a./(v_0*delta_t);
            theta_theory=@(x)sqrt(10-sqrt(x.^6/42+120*x-20));
            theta_theory(fraction_in_theory(delta_t,V_0_matrix,a))
            plot(V_0_matrix*delta_t/(2*a),theta_theory(fraction_in_theory(delta_t,V_0_matrix,a)),'k')
            plot(V_0_matrix*delta_t/(2*a),-theta_theory(fraction_in_theory(delta_t,V_0_matrix,a)),'k')

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


%% Save Variables for Viktor
% save('For_Viktor.mat','a','D','delta_t','dt','k_B','num_transitions_matrix','T','R_matrix','v_0_matrix')




