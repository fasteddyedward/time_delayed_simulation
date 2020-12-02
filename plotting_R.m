clear;
Date='2020.12.1'
nth_take=100
Delta_t_matrix=[2]
T_matrix=[1]
V_0_matrix=[3.5:0.1:6.4]
dt=10^-2; % ms 
intrinsic_delay=0.1 % Intrinsic delay
recalculate_R='no'

%%
R_matrix=[];
for delta_t_index=1:length(Delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(V_0_matrix)
            %             if  nth_take~=7
            
            if 1
                %             if ismember(nth_take,nth_interest)
close all
%% Input File Name
% movie_name=['2020.11.19,dt=10e-3 take ',num2str(nth_take)];
% movie_name=['2020.11.20,dt=10e-3 take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% movie_name=['test3']
% movie_name=['2020.11.25,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(V_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];
[movie_name,'.mat'];
% load([movie_name,'.mat'],'num_transitions','theta_plus','theta_minus')
load([movie_name,'.mat'])
%%
% partition_movie='off';
% N=2;
% v_0=v_0_matrix(v_0_index);
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
    %% Plotting histogram
%     num_bins=100  ;
%     bin_limit=2;
%     figure(81),clf % Note that hist_analysis only works for one particle orbitting a fixed particle at the moment
    
%     [k_trans,theta_plus,theta_minus,num_transitions]=hist_analysis(movie_name,moving_avg,num_bins,bin_limit,Obs_time_steps,delta_t,dt);
%     set(gca, 'YScale', 'linear')
%     saveas(gcf,[movie_name,' (hist).png'])    

%%



%% Appending the matrices
% load([movie_name,'.mat'],'R_mean');
R_matrix=[R_matrix R_mean(2)];

            end
            nth_take=nth_take+1
        end
        nth_take=nth_take+1
    end
    nth_take=nth_take+1
end
%%
% v_0_matrix(1)=[];
if length(Delta_t_matrix)>1
    figure(1);clf;
    hold on
    plot(v_0*Delta_t_matrix/(2*a),R_matrix/(2*a))
    plot(v_0*Delta_t_matrix/(2*a),v_0*Delta_t_matrix* 2/ pi/(2*a))
    xlabel('v_0*\delta t/(2a)')
    ylabel('R/(2a)')
    title(['Orbit Radius v.s. \delta t, v_0= ',num2str(v_0),', T=',num2str(T)])
    axis([0 inf 0 inf])
    % set(gca, 'YScale','log')
    % set(gca, 'XScale','log')
    legend('simulation','v_0*\delta t*2/pi','Location','southeast')
    saveas(gcf,['Orbit Radius v.s. delta_t, v_0=',num2str(v_0),', T=',num2str(T),'.png'])
    %%
    figure(2);clf;
    hold on
    plot(v_0*Delta_t_matrix/(2*a),R_matrix/(2*a))
    plot(v_0*Delta_t_matrix/(2*a),v_0*Delta_t_matrix* 2/ pi/(2*a))
    xlabel('v_0*\delta t/(2a)')
    ylabel('R/(2a)')
    title(['Orbit Radius v.s. \delta t, v_0= ',num2str(v_0),', T=',num2str(T)])
    axis([0 inf 0 inf])
    set(gca, 'YScale','log')
    set(gca, 'XScale','log')
    legend('simulation','v_0*\delta t*2/pi','Location','southeast')
    saveas(gcf,['Orbit Radius v.s. delta_t, v_0=',num2str(v_0),', T=',num2str(T),' log plot.png'])
    %
end





%%
if length(V_0_matrix)>1
    figure(1);clf;
    hold on
    plot(V_0_matrix*delta_t/(2*a),R_matrix/(2*a))
    plot(V_0_matrix*delta_t/(2*a),V_0_matrix*delta_t* 2/ pi/(2*a))
    xlabel('v_0*\delta t/(2a)')
    ylabel('R/(2a)')
    title(['Orbit Radius v.s. v_0, \delta t= ',num2str(delta_t),', T=',num2str(T)])
    axis([0 inf 0 inf])
    % set(gca, 'YScale','log')
    % set(gca, 'XScale','log')
    legend('simulation','v_0*\delta t*2/pi','Location','southeast')
    saveas(gcf,['Orbit Radius v.s. v_0, delta_t=',num2str(delta_t),', T=',num2str(T),' log plot.png'])
    
    %%
    figure(2);clf;
    hold on
    plot(V_0_matrix*delta_t/(2*a),R_matrix/(2*a))
    plot(V_0_matrix*delta_t/(2*a),V_0_matrix*delta_t* 2/ pi/(2*a))
    xlabel('v_0*\delta t/(2a)')
    ylabel('R/(2a)')
    title(['Orbit Radius v.s. v_0, \delta t= ',num2str(delta_t),', T=',num2str(T)])
    axis([0 inf 0 inf])
    set(gca, 'YScale','log')
    set(gca, 'XScale','log')
    legend('simulation','v_0*\delta t*2/pi','Location','southeast')
    saveas(gcf,['Orbit Radius v.s. v_0, delta_t=',num2str(delta_t),', T=',num2str(T),' log plot.png'])
    
    %
end




