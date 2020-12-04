Date='2020.11.24'
nth_take=45
delta_t_matrix=[2]
T_matrix=[1]
v_0_matrix=[7.3:0.1:10]
dt=10^-2; % ms 
% intrinsic_delay=0.1 % Intrinsic delay
% intrinsic_delay=0.1 % Intrinsic delay
%%
% Nth_take=nth_take;
V_0_matrix=v_0_matrix;
Delta_t_matrix=delta_t_matrix;
%%
for delta_t_index=1:length(Delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(V_0_matrix)
            %             if  nth_take~=7
            if 1
%% Input File Name
% movie_name=['test3']
% movie_name=['2020.11.25,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(V_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];
% movie_name=[Date,',dt=10e-3 take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];

[movie_name,'.mat'];
% load([movie_name,'.mat'],'num_transitions','theta_plus','theta_minus')
load([movie_name,'.mat'])
%% Start Calculating Autocorrelation Functions
close all
num_lag=1000;
[acf,lags,bounds] =autocorr(diff(theta(2,:)),'NumLags',num_lag);
plot(0:dt:(num_lag)*dt,acf)
xline(delta_t,'g')
xlabel('\tau')
ylabel('g(\tau)')
title('Autocorrelation function of \Delta\theta')
legend('g(\tau)',['\delta t=',num2str(delta_t)])
            end
            nth_take=nth_take+1;
        end
        nth_take=nth_take+1;
    end
    nth_take=nth_take+1;
end
