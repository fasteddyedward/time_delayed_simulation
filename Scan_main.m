
function [num_transitions_matrix,theta_plus_matrix,theta_minus_matrix,R_matrix,D_omega_matrix,D_theta_matrix,theta_0_matrix,R_recip_matrix,v_0_matrix]= Scan_main(Date,nth_take,delta_t_matrix,T_matrix,theta_0_1D_matrix,dt,intrinsic_delay,Obs_time_steps,D_0)
%% Matrices to append
num_transitions_matrix(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
theta_plus_matrix(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
theta_minus_matrix(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
R_matrix(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
D_omega_matrix(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
D_theta_matrix(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
theta_0_matrix(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
R_recip_matrix(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
v_0_matrix(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;

%%
for delta_t_index=1:length(delta_t_matrix)
    for T_index=1:length(T_matrix)
        for theta_0_index=1:length(theta_0_1D_matrix)
            if 1
                %             if ismember(nth_take,nth_interest)
close all
%% Output File Name
% movie_name=['2020.11.19,dt=',num2str(dt),' take ',num2str(nth_take)];
% movie_name=['2020.11.26,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', theta_0_1D=',num2str(theta_0_1D_matrix(theta_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% movie_name=['test3']
warning('Have you modified the file name?')
%% Save mat file? (If to run on cluster, pick 'no'
save_file='no' ;
%% Setup for Running the program
N=2; % total number of particles in the simulation
delta_t=delta_t_matrix(delta_t_index); % ms


partition_time_steps=Obs_time_steps;

%% State if the particles are fixed, 0 for mobile, 1 for fixed
fixed_flag(1:N)=0;
fixed_flag(1)=1; % particle 1 is fixed
%% For hardcore interaction
% hard_collision='method_2' % or 'off'
hard_collision='test_no_elastic';
% hard_collision='method_3';
a=5 ;% Particle radius, typically 1 micrometer
b= 0.5 ;% for particle retreat during relaxation period
%% Coefficients and parameters
v_0= theta_0_1D_matrix(theta_0_index)*(2*a)/delta_t; % mm/ms
% v_0= theta_0_1D_matrix(theta_0_index)*(2*a)/sqrt(delta_t); % mm/ms
% v_0=pi/(2*theta_0_1D_matrix(theta_0_index)*delta_t);
T=T_matrix(T_index); % Kelvin 
gamma=1;
k_B=1;
% gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
% k_B=10^-23; % Boltzmann constant
% D_0=k_B*T/gamma; % Diffusitivity

%% Warning for Program Setup
if partition_time_steps>Obs_time_steps
    warning(['Please choose a partitioned time step larger than Obs_time_steps= ',num2str(Obs_time_steps)])
    pause
elseif partition_time_steps<delta_t/dt
    warning(['Please choose a partitioned time step larger than delta/dt= ',num2str(delta_t/dt)])
    pause
end
if (~exist('time', 'var'))
    warning(['Some variables e.g. time,x,y has been existent, delete the old folders if the simulations crash.'])
end%if

        
%% Initial positions
x_init(1:N)=0;
y_init(1:N)=0;
c=5;
for i=1:N
    x_init(i)=i*2*5*c;
    y_init(i)=i*2*5*c;
end
% y_init(3)=1*10^1;
y_init(3)=2*10^1;


    %% Warning for initial position
    for i=1:N
        for j=i+1:N
            diff_x=x_init(i)-x_init(j);
            diff_y=y_init(i)-y_init(j);
            if diff_x^2+diff_y^2 < (2*a)^2
                warning('The particles initial position are too close!')
                warning('If you decided to continue, the particles will appear to be stuck together; the hard core interaction does not have the ability to avoid this yet.')
                pause
            end
        end
    end
%% Start calculating finite element numericals for the equation of motion
time_simulation_start=tic;
% movie_x_max=0;movie_x_min=0;movie_y_max=0;movie_y_min=0; %% Boundary
% for the movie
[x,y,v_x,v_y,movie_x_min,movie_x_max,movie_y_min,movie_y_max,x_final,y_final]=simulation(movie_name,Obs_time_steps,partition_time_steps,0,0,0,0,N,delta_t,dt,v_0,T,gamma,k_B,D_0,x_init,y_init,hard_collision,a,b,intrinsic_delay,fixed_flag,save_file);
time_simulation=toc(time_simulation_start)
%% Start analyzing theta and omega for one-fixed-one-mobile-particle system
[theta,omega,R_mean,D_theta,D_omega2,theta_0,R_recip]= get_D_eff(x,y,v_x(2,:),v_y(2,:),Obs_time_steps,delta_t,dt,v_0);
[theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0);
%% Appending Matrices
figure;
% plot(theta)
theta_0
num_transitions_theta

num_transitions_matrix(delta_t_index,theta_0_index)=num_transitions_theta;
theta_plus_matrix(delta_t_index,theta_0_index)=theta_plus;
theta_minus_matrix(delta_t_index,theta_0_index)=theta_minus;
R_matrix(delta_t_index,theta_0_index)=R_mean;
D_omega_matrix(delta_t_index,theta_0_index)=D_omega2;
D_theta_matrix(delta_t_index,theta_0_index)=D_theta;
theta_0_matrix(delta_t_index,theta_0_index)=theta_0;
R_recip_matrix(delta_t_index,theta_0_index)=R_recip;
v_0_matrix(delta_t_index,theta_0_index)=v_0;

nth_take=nth_take+1

%% Clearing Unwanted Timing Variables
clear Analyze_rot combine_data_partitions_start making_movies time_simulation_start

            end
%             nth_take=nth_take+1
        end
%         nth_take=nth_take+1
    end
%     nth_take=nth_take+1
end

end

