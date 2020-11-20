clear;
nth_take=7
% delta_t_matrix=[0.01 0.1 1 5 10]
% T_matrix=[0.01 0.1 1 10]
% v_0_matrix=[0.01 0.1 1 10]
delta_t_matrix=[2]
T_matrix=[1]
v_0_matrix=[3.5:0.1:7]
dt=10^-2
intrinsic_delay=0 % Intrinsic delay
num_transitions_matrix=[];
theta_plus_matrix=[];
theta_minus_matrix=[];
for delta_t_index=1:length(delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(v_0_matrix)
            if  nth_take~=7
% if 1
                %             if ismember(nth_take,nth_interest)
close all
%% Output File Name
% movie_name=['2020.11.19,dt=10e-3 take ',num2str(nth_take)];
% movie_name=['2020.11.20,dt=10e-3 take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% movie_name=['test3']
movie_name=['2020.11.20,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];

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
close all 
figure(1)
% plot(v_0_matrix(2:end),num_transitions_matrix,'.')
plot(v_0_matrix(2:end),num_transitions_matrix,'-')
figure(2)
hold on
plot(v_0_matrix(2:end),theta_plus_matrix,'.')
plot(v_0_matrix(2:end),theta_minus_matrix,'.')
% plot(v_0_matrix(2:end),theta_plus_matrix,'-')
% plot(v_0_matrix(2:end),theta_minus_matrix,'-')


