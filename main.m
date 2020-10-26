%% This file runs modulized_time_delay_proto
%% 2020.10.14 to make the videos with several tries
for nth_take=49
clearvars -except nth_take
close all
%% Output File Name
% movie_name=['2020.10.26,dt=10e-3 take ',num2str(nth_take)];
movie_name=['test3']
warning('Have you modified the file name?')

%% Setup for Running the program
N=3; % number of particles in the play
delta_t=0.1; % ms
dt=10^-3; % ms 
% Obs_time=Obs_time_steps*dt;
Obs_time_steps=10^4

partition_time_steps=10^3
partition_movie='no'

%% Coefficients and parameters
v_0= 1; % mm/ms
T=0.01; % Kelvin 
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
for i=1:N
    x_init(i)=i*10^1;
    y_init(i)=i*10^1;
end
y_init(3)=2*10^1;

%% Start calculating finite element numericals for the equation of motion
time_simulation_start=tic;
% movie_x_max=0;movie_x_min=0;movie_y_max=0;movie_y_min=0; %% Boundary
% for the movie
[x,y,v_x,v_y,movie_x_min,movie_x_max,movie_y_min,movie_y_max,x_final,y_final]=simulation(movie_name,Obs_time_steps,partition_time_steps,0,0,0,0,N,delta_t,dt,v_0,T,gamma,k_B,D,x_init,y_init);
time_simulation=toc(time_simulation_start)

%% Putting all the partitioned files into one folder
if (~exist(movie_name, 'dir')); mkdir(movie_name); end%if
for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
    movefile([movie_name,' partition_',num2str(lth_partition),'.mat'],movie_name);
end
%% Loading all the recorded positions and save to [movie_name,'.mat'](! Might cause memery overflow)
switch partition_movie
    case 'no'
        combine_data_partitions_start=tic;
        [x,y,v_x,v_y,time]=combine_partitions(movie_name,Obs_time_steps,partition_time_steps,delta_t,dt);
        % movie_x_max=max(x,[],'all');        % movie_y_max=max(y,[],'all');        % movie_x_min=min(x,[],'all');        % movie_y_min=min(y,[],'all');
        save([movie_name,'.mat']);
        clear x y v_x v_y time
        time_combine_data_partitions=toc(combine_data_partitions_start)
end

%% Finding axis_scale  =   [movie_x_min  movie_x_max  movie_y_min  movie_y_max]
if (~exist('movie_x_max','var'))
    movie_x_max=0;movie_x_min=0;movie_y_max=0;movie_y_min=0;
    cd(movie_name)
    for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
        load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y')
        movie_x_max=max([movie_x_max max(x)]);    movie_x_min=min([movie_x_min min(x)]);    movie_y_max=max([movie_y_max max(y)]);    movie_y_min=min([movie_y_min min(y)]);
    end
    cd ..
end
axis_scale=[movie_x_min movie_x_max movie_y_min movie_y_max];




%% Making Movies
    %% Parameters for making the movies
    making_movies=tic;
    magnify=1000    ;
    control_animation_interval=10^4     ; % Record one frame in every ____ frame
    movie_create='on'   ;
    ghost='off'      ;
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
    Movie_Vector=Make_Movie(movie_name,N,dt,Obs_time_steps,partition_time_steps,delta_t,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,partition_movie,movie_x_min, movie_x_max,movie_y_min,movie_y_max);
    
    %% Save movie
    frame_rate=10;
    Save_created_movie(movie_name,movie_create,control_animation_interval,Obs_time_steps,delta_t,dt,frame_rate,Movie_Vector)
    %% Saving work space for Movie
    save([movie_name,'.mat'],'Movie_Vector','magnify','control_animation_interval','movie_create','ghost','axis_choice','leave_trace','-append')
    clear Movie_Vector
    time_making_movies=toc(making_movies)



%% Rotational Analysis: calculates and plots v_omega
Analyze_rot=tic;
moving_avg=1000 ;
Rotational_Analysis(movie_name,partition_movie,N,v_0,Obs_time_steps,partition_time_steps,delta_t,dt,moving_avg);
time_analyze_rot=toc(Analyze_rot)
    
%% check point to make sure time is not wrong, or %% Drawing v_omega will break
time=(1:Obs_time_steps+delta_t/dt)*dt;
save([movie_name,'.mat'],'time','-append')
%% Plotting v_omega

moving_avg=1000 ;
figure(98); clf %% Would be same as figure(99) if partition movie='no'
plot_v_omega(N,delta_t,movie_name,moving_avg)
saveas(gcf,[movie_name,'.png'])
%% Clearing Unwanted Timing Variables
clear Analyze_rot combine_data_partitions_start making_movies time_simulation_start
end

