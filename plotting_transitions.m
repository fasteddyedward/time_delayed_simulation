clear;
nth_take=100
% delta_t_matrix=[0.01 0.1 1 5 10]
% T_matrix=[0.01 0.1 1 10]
% v_0_matrix=[0.01 0.1 1 10]
% delta_t_matrix=[2]
delta_t_matrix=[0.5:0.5:16]
T_matrix=[1]
% v_0_matrix=[3.5:0.1:7]
v_0_matrix=5
dt=10^-2
intrinsic_delay=0 % Intrinsic delay
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
%% Output File Name
% movie_name=['2020.11.19,dt=10e-3 take ',num2str(nth_take)];
% movie_name=['2020.11.20,dt=10e-3 take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% movie_name=['test3']
movie_name=['2020.11.25,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];

[movie_name,'.mat'];
% load([movie_name,'.mat'],'num_transitions','theta_plus','theta_minus')
load([movie_name,'.mat'])
%%
% partition_movie='off';
% N=2;
% v_0=v_0_matrix(v_0_index);
%% Calculating the histograms and stuff
    theta=0;
    Analyze_theta=tic;
    moving_avg=1 ;
    plot_rot='no';
    Theta_Analysis_Fixed_Center(movie_name,partition_movie,N,v_0,Obs_time_steps,partition_time_steps,delta_t,dt,moving_avg,plot_rot);
    time_analyze_theta=toc(Analyze_theta)
    %% Plotting histogram
    num_bins=100  ;
    bin_limit=2;
    figure(81),clf % Note that hist_analysis only works for one particle orbitting a fixed particle at the moment
    
    [k_trans,theta_plus,theta_minus,num_transitions]=hist_analysis(movie_name,moving_avg,num_bins,bin_limit,Obs_time_steps,delta_t,dt);
    set(gca, 'YScale', 'linear')
    saveas(gcf,[movie_name,' (hist).png'])    

%% Appending the matrices
num_transitions_matrix=[num_transitions_matrix num_transitions];
theta_plus_matrix=[theta_plus_matrix, theta_plus];
theta_minus_matrix=[theta_minus_matrix, theta_minus];

            end
            nth_take=nth_take+1
        end
        nth_take=nth_take+1
    end
    nth_take=nth_take+1
end
%%
v_0_matrix(1)=[];
%% Delta_t
if length(delta_t_matrix)>1
    close all
    %%
    figure(3) ;clf
    
    hold on
    plot(v_0*delta_t_matrix/(2*a),theta_plus_matrix,'.')
    plot(v_0*delta_t_matrix/(2*a),theta_minus_matrix,'.')
    % plot(v_0_matrix,theta_plus_matrix,'-')
    % plot(v_0_matrix,theta_minus_matrix,'-')
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
end
%% V_0
if length(v_0_matrix)>1
    close all
    figure(1);clf
    % plot(v_0_matrix,num_transitions_matrix,'.')
    plot(v_0_matrix,num_transitions_matrix,'-')
    title(['Transition rates, \delta t= ',num2str(delta_t),', T=',num2str(T)])
    xlabel('v_0')
    ylabel('Number of transitions')
    figure(1)
    hold on
    set(gca, 'YScale', 'linear')
    parameter_1=3000 % a in paper
    parameter_2=3 % b in paper
    omega_0=@(a,v_0) v_0./(a);
    k=@(omega_0) parameter_1*exp(-parameter_2*3/2*(omega_0*delta_t-1).^2/(D*delta_t));
    trans_theory=k(omega_0(a,v_0_matrix));
    
    f=polyfit(v_0_matrix,log(trans_theory),2);
    plot(v_0_matrix,exp(polyval(f,v_0_matrix)))
    set(gca, 'YScale', 'log')
    legend('simulation','theory')
    
    %%
    % figure(2);clf
    % hold on
    % plot(v_0_matrix,theta_plus_matrix,'.')
    % plot(v_0_matrix,theta_minus_matrix,'.')
    % % plot(v_0_matrix,theta_plus_matrix,'-')
    % % plot(v_0_matrix,theta_minus_matrix,'-')
    % title(['Bifurcation Diagram, \delta t= ',num2str(delta_t),', T=',num2str(T)])
    % xlabel('v_0')
    % ylabel('\theta (rad)')
    % figure(2)
    % x=@(delta_t,v_0,a)2*a./(v_0*delta_t);
    % theta_theory=@(x)sqrt(10-sqrt(x.^6/42+120*x-20));
    % theta_theory(x(delta_t,v_0_matrix,a))
    % plot(v_0_matrix,theta_theory(x(delta_t,v_0_matrix,a)),'k')
    % plot(v_0_matrix,-theta_theory(x(delta_t,v_0_matrix,a)),'k')
    %%
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
    
    %%
    figure(4)
    plot(v_0_matrix,num_transitions_matrix,'-')
    title('Data points')
    figure(5)
    plot(v_0_matrix,num_transitions_matrix,'-')
    title('Data points')
    % hold off
    
    F = @(para,data)para(1)*exp(-para(2)*3/2*(data/(2*a)*delta_t-1).^2/(D*delta_t));
    F_log= @(para,data)log(para(1)*exp(-para(2)*3/2*(data/(2*a)*delta_t-1).^2/(D*delta_t)));
    x0 = [1 1];
    [para,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,v_0_matrix,num_transitions_matrix);
    [para_log,resnorm,~,exitflag,output] = lsqcurvefit(F_log,x0,v_0_matrix,log(num_transitions_matrix));
    
    figure(4)
    hold on
    plot(v_0_matrix,F(para,v_0_matrix))
    hold off
    axis([3.5 7 0 inf])
    
    set(gca, 'YScale', 'linear')
    
    figure(5)
    hold on
    plot(v_0_matrix,exp(F_log(para_log,v_0_matrix)))
    hold off
    set(gca, 'YScale', 'linear')
    axis([3.5 7 0 inf])
end
    



