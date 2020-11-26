clear;
%% Parameters for importing mat files.

Date='2020.11.24'
nth_take=7
delta_t_matrix=2
T_matrix=[1]
v_0_matrix=[3.5:0.1:10]
dt=10^-2

%% Execution Parameters
recalculate_theta='no'
recalculate_hist='no'
draw_hist='no'


%%
num_transitions_matrix=[];
theta_plus_matrix=[];
theta_minus_matrix=[];
for delta_t_index=1:length(delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(v_0_matrix)
            %             if  nth_take~=7
            if 1
                %             if ismember(nth_take,nth_interest)
close all
%% Input File Name
% movie_name=['test3']
% movie_name=['2020.11.25,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% movie_name=[Date,',dt=10e-3 take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];

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
        nth_take
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
        nth_take
end

%% Appending the matrices
num_transitions_matrix=[num_transitions_matrix num_transitions];
theta_plus_matrix=[theta_plus_matrix, theta_plus];
theta_minus_matrix=[theta_minus_matrix, theta_minus];
nth_take
            end
            nth_take=nth_take+1;
        end
        nth_take=nth_take+1;
    end
    nth_take=nth_take+1;
end
%%
% v_0_matrix(1)=[];
%% Delta_t
if length(delta_t_matrix)>1
    close all
    %% Bifurcation Diagram
    figure(3) ;clf
    
    hold on
    plot(v_0*delta_t_matrix/(2*a),theta_plus_matrix,'.')
    plot(v_0*delta_t_matrix/(2*a),theta_minus_matrix,'.')
    title(['Bifurcation Diagram, v_0= ',num2str(v_0),', T=',num2str(T)])
    xlabel('v_0*\delta t/R')
    ylabel('\theta (rad)')
    % figure(2)
    x=@(delta_t,v_0,a)2*a./(v_0*delta_t);
    theta_theory=@(x)sqrt(10-sqrt(x.^6/42+120*x-20));
    theta_theory(x(delta_t_matrix,v_0,a))
    plot(v_0*delta_t_matrix/(2*a), theta_theory(x(delta_t_matrix,v_0,a)),'k')
    plot(v_0*delta_t_matrix/(2*a),-theta_theory(x(delta_t_matrix,v_0,a)),'k')
    saveas(gcf,['Bifurcation Diagram, v_0= ',num2str(v_0),', T=',num2str(T),'.png'])   
    %% Transition Rates
    figure(4)
    plot(delta_t_matrix,num_transitions_matrix,'-')
    title('Data points')
    set(gca, 'YScale', 'linear')
    figure(5)
    plot(delta_t_matrix,num_transitions_matrix,'-')
    title('Data points')
    set(gca, 'YScale', 'log')
    
    %% Fitting the intersted part
    delta_t_interest=delta_t_matrix(delta_t_matrix>2*a/v_0);
    num_transitions_interest=num_transitions_matrix(delta_t_matrix>2*a/v_0);
    figure(6);clf;
    hold on
    plot(delta_t_interest,num_transitions_interest)
    a1=60000
    a2=0.3
    f=@(x)(a1*exp(-a2*(x-1).^2)./x)
%     plot(delta_t_interest,f(delta_t_interest))
%     legend('simulation','theory')
    omega_0=v_0/(2*a)
%     row=1:8
%     fit_curve_trans(delta_t_interest(row), num_transitions_interest(row))
%     set(gca, 'YScale', 'linear')
end


%% V_0
if length(v_0_matrix)>1
    %% Bifurcation Diagram
    figure(3) ;clf
    
    hold on
    plot(v_0_matrix*delta_t/(2*a),theta_plus_matrix,'.')
    plot(v_0_matrix*delta_t/(2*a),theta_minus_matrix,'.')
    % plot(v_0_matrix,theta_plus_matrix,'-')
    % plot(v_0_matrix,theta_minus_matrix,'-')
    title(['Bifurcation Diagram, \delta t= ',num2str(delta_t),', T=',num2str(T)])
    xlabel('v_0*\delta t/R')
    ylabel('\theta (rad)')
    % figure(2)
    x=@(delta_t,v_0,a)2*a./(v_0*delta_t);
    theta_theory=@(x)sqrt(10-sqrt(x.^6/42+120*x-20));
    theta_theory(x(delta_t,v_0_matrix,a))
    plot(v_0_matrix*delta_t/(2*a),theta_theory(x(delta_t,v_0_matrix,a)),'k')
    plot(v_0_matrix*delta_t/(2*a),-theta_theory(x(delta_t,v_0_matrix,a)),'k')
    
    %% Transition Rates
    figure(4)
    plot(v_0_matrix*delta_t/(2*a),num_transitions_matrix,'-')
    xlabel('v_0*\delta t/R')
    title('Data points')
    %% Interested Transition Rates
    figure(6);clf
    v_0_interest=v_0_matrix(v_0_matrix>2*a/delta_t);
    num_transitions_interest=num_transitions_matrix(v_0_matrix>2*a/delta_t);
%     v_0_interest(1)=[];
%     num_transitions_interest(1)=[];
    % Getting rid of the 0 tail
%     v_0_interest_head=v_0_interest(num_transitions_interest~=0);
%     num_transitions_interest_head=num_transitions_interest(num_transitions_interest~=0);
%     v_0_interest_head=v_0_interest(1:25);
%     num_transitions_interest_head=num_transitions_interest(1:25);
%     plot(v_0_interest_head,num_transitions_interest_head,'o')
%     figure(5)
%     plot(v_0_matrix,num_transitions_matrix,'-')
%     title('Data points')
%     % hold off
    plot(v_0_interest,num_transitions_interest)
    hold on
%% Self-picked transition rates
v_0_interest=v_0_matrix(v_0_matrix>6.5);
num_trans_interest=num_transitions_matrix(v_0_matrix>6.5);
%% Try fitting the first time
%     close all
%     figure(1);clf
%     % plot(v_0_matrix,num_transitions_matrix,'.')
%     plot(v_0_matrix,num_transitions_matrix,'-')
%     title(['Transition rates, \delta t= ',num2str(delta_t),', T=',num2str(T)])
%     xlabel('v_0')
%     ylabel('Number of transitions')
%     figure(1)
%     hold on
%     set(gca, 'YScale', 'linear')
%     parameter_1=3000 % a in paper
%     parameter_2=3 % b in paper
%     omega_0=@(a,v_0) v_0./(a);
%     k=@(omega_0) parameter_1*exp(-parameter_2*3/2*(omega_0*delta_t-1).^2/(D*delta_t));
%     trans_theory=k(omega_0(a,v_0_matrix));
%
%     f=polyfit(v_0_matrix,log(trans_theory),2);
%     plot(v_0_matrix,exp(polyval(f,v_0_matrix)))
%     set(gca, 'YScale', 'log')
%     legend('simulation','theory')
%
%% Try fitting the second time
%
%     F = @(para,data)para(1)*exp(-para(2)*3/2*(data/(2*a)*delta_t-1).^2/(D*delta_t));
%     F_log= @(para,data)log(para(1)*exp(-para(2)*3/2*(data/(2*a)*delta_t-1).^2/(D*delta_t)));
%     x0 = [1 1];
%     [para,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,v_0_matrix,num_transitions_matrix);
%     [para_log,resnorm,~,exitflag,output] = lsqcurvefit(F_log,x0,v_0_matrix,log(num_transitions_matrix));
%
%     figure(4)
%     hold on
%     plot(v_0_matrix,F(para,v_0_matrix))
%     hold off
%     axis([3.5 7 0 inf])
%
%     set(gca, 'YScale', 'linear')
%
%     figure(5)
%     hold on
%     plot(v_0_matrix,exp(F_log(para_log,v_0_matrix)))
%     hold off
%     set(gca, 'YScale', 'linear')
%     axis([3.5 7 0 inf])
end






