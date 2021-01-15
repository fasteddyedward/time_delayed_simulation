function [num_transitions_matrix,theta_plus_matrix,theta_minus_matrix,R_matrix,D_omega_matrix,D_theta_matrix,theta_0_matrix]= main(Date,nth_take,delta_t_matrix,T_matrix,v_0_matrix,dt,intrinsic_delay,Obs_time_steps)
% Date,nth_take,delta_t_matrix,T_matrix,v_0_matrix,dt,intrinsic_delay,Obs_time_steps
% function main(testparameter1,testparameter2)
%% This file simulates the time-delayed interaction of active brownian particles, plots the movies, and analyzes the rotation
% clearvars -except nth_take


% Date='2021.1.15' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=1 % for a=5, particle 1 fixed
% delta_t_matrix=[1.8:0.1:4.0]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-2
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5


% Date='2021.1.6'
% nth_take=100
% delta_t_matrix=2
% T_matrix=[1]
% % v_0_matrix=[0.5:0.5:14]
% v_0_matrix=[5:0.5:8]
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5
%% Matrices to append
num_transitions_matrix=[];
theta_plus_matrix=[];
theta_minus_matrix=[];
R_matrix=[];
D_omega_matrix=[];
D_theta_matrix=[];
theta_0_matrix=[];

%%
for delta_t_index=1:length(delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(v_0_matrix)
            if 1
                %             if ismember(nth_take,nth_interest)
close all
%% Output File Name
% movie_name=['2020.11.19,dt=',num2str(dt),' take ',num2str(nth_take)];
% movie_name=['2020.11.26,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% movie_name=['test3']
warning('Have you modified the file name?')
%% Save mat file? (If to run on cluster, pick 'no'
save_file='no' ;
%% Setup for Running the program
N=2; % total number of particles in the simulation
delta_t=delta_t_matrix(delta_t_index); % ms
% Obs_time=Obs_time_steps*dt;

partition_time_steps=Obs_time_steps;
partition_movie='no';
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
v_0= v_0_matrix(v_0_index); % mm/ms
T=T_matrix(T_index); % Kelvin 
gamma=1;
k_B=1;
% gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
% k_B=10^-23; % Boltzmann constant
D=k_B*T/gamma; % Diffusitivity

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
[x,y,v_x,v_y,movie_x_min,movie_x_max,movie_y_min,movie_y_max,x_final,y_final]=simulation(movie_name,Obs_time_steps,partition_time_steps,0,0,0,0,N,delta_t,dt,v_0,T,gamma,k_B,D,x_init,y_init,hard_collision,a,b,intrinsic_delay,fixed_flag,save_file);
time_simulation=toc(time_simulation_start)

%% Putting all the partitioned files into one folder
switch save_file
    case 'yes'
        if (~exist(movie_name, 'dir')); mkdir(movie_name); end%if
        for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
            movefile([movie_name,' partition_',num2str(lth_partition),'.mat'],movie_name);
        end
    case 'no'
end
%% Loading all the recorded positions and save to [movie_name,'.mat'](! Might cause memery overflow)
switch partition_movie
    case 'no'
        %         combine_data_partitions_start=tic;
        switch save_file
            case 'yes'
                [x,y,v_x,v_y,time]=combine_partitions(movie_name,Obs_time_steps,partition_time_steps,delta_t,dt);
                % movie_x_max=max(x,[],'all');        % movie_y_max=max(y,[],'all');        % movie_x_min=min(x,[],'all');        % movie_y_min=min(y,[],'all');
                save([movie_name,'.mat']);
                clear x y v_x v_y time
            case 'no'
                time=(1:Obs_time_steps+delta_t/dt)*dt;
        end
%         time_combine_data_partitions=toc(combine_data_partitions_start)
end
%% Removing folder; these files are useless anyway.
if (exist(movie_name, 'dir')); rmdir(movie_name,'s'); end
%% Finding axis_scale  =   [movie_x_min  movie_x_max  movie_y_min  movie_y_max]

if (~exist('movie_x_max','var'))
    switch save_file
        case 'yes'
            movie_x_max=0;movie_x_min=0;movie_y_max=0;movie_y_min=0;
            cd(movie_name)
            for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
                load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y')
                movie_x_max=max([movie_x_max max(x)]);    movie_x_min=min([movie_x_min min(x)]);    movie_y_max=max([movie_y_max max(y)]);    movie_y_min=min([movie_y_min min(y)]);
            end
            cd ..
        case 'no'
    end
else
    switch save_file
        case 'yes'
        case 'no'
            movie_x_max=max(max(x)); movie_x_min=min(min(x)); movie_y_max=max(max(y)); movie_y_min=min(min(y));
    end
end
    % Making the plot look square.
    movie_x_cent=(movie_x_max+movie_x_min)/2;
    movie_y_cent=(movie_y_max+movie_y_min)/2;
    plot_width=max(movie_x_max-movie_x_min,movie_y_max-movie_y_min)+4*a;
    movie_x_min=movie_x_cent-plot_width/2;
    movie_x_max=movie_x_cent+plot_width/2;
    movie_y_min=movie_y_cent-plot_width/2;
    movie_y_max=movie_y_cent+plot_width/2;
    axis_scale=[movie_x_min movie_x_max movie_y_min movie_y_max];




    %% Making Movies
    %% Parameters for making the movies
    %     making_movies=tic;
    magnify=1000    ;
    control_animation_interval=10^4    ; % Record one frame in every ____ frame
    movie_create='off'   ;
    ghost='off'      ;
    force_tracks='off';
    axis_choice='lab'; %'cm' or 'lab'
    leave_trace='off'       ;
    close all
    %% Warnings
    if control_animation_interval>Obs_time_steps+dt*delta_t
        warning(['The number of movie frames is 0; Please modify control_animation_interval to a number smaller than Obs_time_steps= ',num2str(Obs_time_steps)])
    end
    (Obs_time_steps+delta_t/dt)/control_animation_interval      ;
    if ((Obs_time_steps+delta_t/dt)/control_animation_interval>100)
        warning(['The plotting will take about ',num2str((Obs_time_steps+delta_t/dt)/control_animation_interval*0.15),' seconds. For faster time please choose different control_animation_interval'])
    end
    % 1000 takes about 150 seconds.
    
    %% Start Making Movie
    %     if strcmp(partition_movie,'yes')
    switch partition_movie
        case 'yes'
            Movie_Vector=Make_Movie(movie_name,N,dt,Obs_time_steps,partition_time_steps,delta_t,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,partition_movie,movie_x_min, movie_x_max,movie_y_min,movie_y_max,a,force_tracks);
        case 'no'
            %             elseif strcmp(partition_movie,'no')
            switch save_file
                
                %             if strcmp(save_file,'yes')
                case 'yes'
                    Movie_Vector=Make_Movie(movie_name,N,dt,Obs_time_steps,partition_time_steps,delta_t,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,partition_movie,movie_x_min, movie_x_max,movie_y_min,movie_y_max,a,force_tracks);
                case 'no' % Then we don't need Make_Movie to load the mat files
                    %             elseif strcmp(save_file,'no')
                    [Movie_Vector]=make_movies_plots(N,delta_t,dt,partition_time_steps,x,y,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,a,force_tracks);
            end
    end
    
    %     Movie_Vector=Make_Movie(movie_name,N,dt,Obs_time_steps,partition_time_steps,delta_t,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,partition_movie,movie_x_min, movie_x_max,movie_y_min,movie_y_max,a,force_tracks,x,y,time,save_file);
    
    %% Save movie
    frame_rate=10;
    Save_created_movie(movie_name,movie_create,control_animation_interval,Obs_time_steps,delta_t,dt,frame_rate,Movie_Vector)
    %% Saving work space for Movie
    switch save_file
        case 'yes'
            save([movie_name,'.mat'],'Movie_Vector','magnify','control_animation_interval','movie_create','ghost','axis_choice','leave_trace','-append')
            clear Movie_Vector
        case 'no'
    end
    %     time_making_movies=toc(making_movies)
    %% Plotting x and y of the N particles
    %     plot_x_y(movie_name)
%% check point to make sure time is not wrong, or %% Drawing v_omega will break
switch save_file
    case 'yes'
        time=(1:Obs_time_steps+delta_t/dt)*dt;
        save([movie_name,'.mat'],'time','-append')
    case 'no'
end

%% Start analyzing theta and omega for one-fixed-one-mobile-particle system
[theta,omega,R_mean,D_theta,D_omega2,theta_0]= get_D_eff(x,y,v_x(2,:),v_y(2,:),Obs_time_steps,delta_t,dt,v_0);
[theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0);
% theta_0
% theta_plus
[omega_stable,k_trans_omega,omega_plus,omega_minus,num_transitions_omega]=find_omega_plus(omega,theta_0);

%% Appending Matrices
theta_0
num_transitions_theta
figure;plot(theta);
num_transitions_matrix=[num_transitions_matrix num_transitions_theta];
theta_plus_matrix=[theta_plus_matrix theta_plus];
theta_minus_matrix=[theta_minus_matrix theta_minus];
R_matrix=[R_matrix R_mean];
D_omega_matrix=[D_omega_matrix D_omega2];
D_theta_matrix=[D_theta_matrix D_theta];
theta_0_matrix=[theta_0_matrix theta_0];

%% Clearing Unwanted Timing Variables
clear Analyze_rot combine_data_partitions_start making_movies time_simulation_start

            end
            nth_take=nth_take+1
        end
        nth_take=nth_take+1
    end
    nth_take=nth_take+1
end

end
