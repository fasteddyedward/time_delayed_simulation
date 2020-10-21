%% This file runs modulized_time_delay_proto
%% 2020.10.14 to make the videos with several tries
for nth_take=4
clearvars -except nth_take
close all
%% Output File Name
movie_name=['2020.10.21,dt=10e-3 take ',num2str(nth_take)];
warning('Have you modified the file name?')
% pause
%% Output File Name
% movie_name=['2020.10.16,dt=10e-3 take ',num2str(nth_take)];
% answer = questdlg(['Is the file name correct? ',movie_name], ...
% 	'File Name Check', ...
% 	'Yes','No, I will  modify it','No, I will modify it ');
% switch answer
%     case 'Yes'
%         'Start running program...'
%     case 'No, I will  modify it'
%         error('Please modify the file name')
% end
% tic
%% Setup for Running the program
N=3; % number of particles in the play
delta_t=0.5; % ms
dt=10^-3; % ms 
% Obs_time=Obs_time_steps*dt;
Obs_time_steps=10^4   ;
partition_time_steps=10^3  ;
if partition_time_steps>Obs_time_steps
    warning(['Please choose a partitioned time step larger than Obs_time_steps= ',num2str(Obs_time_steps)])
    pause
end
if partition_time_steps<delta_t/dt
    warning(['Please choose a partitioned time step larger than delta/dt= ',num2str(delta_t/dt)])
    pause
end
%% Coefficients and parameters
v_0= 10^-6; % mm/ms
% gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
T=300; % Kelvin 
% k_B=10^-23; % Boltzmann constant
% D=k_B*T/gamma; % Diffusitivity


%% Initial positions
x_init(1:N)=0;
y_init(1:N)=0;
for i=1:N
    x_init(i)=i*10^-4;
    y_init(i)=i*10^-4;
end
y_init(3)=2*10^-4;
%% Boundary of the particles (for making the movie)
% movie_x_max=0;movie_x_min=0;movie_y_max=0;movie_y_min=0;


%% Start calculating finite element numericals for the equation of motion
time_simulation_start=tic;
    %% Making Folder for the Partitioned Files
    if (~exist(movie_name, 'dir')); mkdir(movie_name); end%if
    %% First stage: Diffusion. t=0 ~ delta_t (lth_partition=1)
    % Obs_time_Steps in this case is round(delta_t/dt)
    ['lth_partition= ',num2str(1),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
    [x,y,v_x,v_y,~]=first_stage_pure_diffusion(N,delta_t,dt,T,x_init,y_init);
%     movie_x_max=max([movie_x_max max(x)]);    movie_x_min=min([movie_x_min min(x)]);    movie_y_max=max([movie_y_max max(y)]);    movie_y_min=min([movie_y_min min(y)]);
    save([movie_name,' partition_',num2str(1),'.mat']); %% Takes ~ 0.16 sec
    x_temp=x;
    y_temp=y;
    %% Second stage: Delayed interaction starts. t=delta_t~Obs_time
    for lth_partition=2:round(Obs_time_steps/partition_time_steps)+1
        ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
        [x,y,~,~,v_x,v_y,~,~,~]=second_stage_delayed_int(N,delta_t,dt,partition_time_steps,v_0,T,x_temp,y_temp,lth_partition);
%         movie_x_max=max([movie_x_max max(x)]);        movie_x_min=min([movie_x_min min(x)]);        movie_y_max=max([movie_y_max max(y)]);        movie_y_min=min([movie_y_min min(y)]);
        save([movie_name,' partition_',num2str(lth_partition),'.mat']);
        x_temp=x(:,end-delta_t/dt:end); % This last postitions of delayed time delta_t goes to the next round for simulation
        y_temp=y(:,end-delta_t/dt:end);
        clear time
    end
clear x_temp y_temp
time_simulation=toc(time_simulation_start)



%% Loading all the recorded positions (! Might cause memery overflow)
combine_data_partitions_start=tic;
[x,y,v_x,v_y,time]=combine_partitions(movie_name,Obs_time_steps,partition_time_steps,delta_t,dt);
movie_x_max=max(x,[],'all');
movie_y_max=max(y,[],'all');
movie_x_min=min(x,[],'all');
movie_y_min=min(y,[],'all');
save([movie_name,'.mat']);
clear x y v_x v_y time
time_combine_data_partitions=toc(combine_data_partitions_start)

%% Parameters for making the movies
making_movies=tic;
magnify=1000    ;
control_animation_interval=10^4     ; % Record one frame in every ____ frame
movie_create='on'   ;
ghost='on'      ;
axis_choice='lab'; %'cm' or 'lab' 
leave_trace='off'       ;
close all
if control_animation_interval>Obs_time_steps+dt*delta_t
    warning(['The number of movie frames is 0; Please modify control_animation_interval to a number smaller than Obs_time_steps= ',num2str(Obs_time_steps)])
end
    (Obs_time_steps+delta_t/dt)/control_animation_interval      ;
if ((Obs_time_steps+delta_t/dt)/control_animation_interval>100)
    warning(['The plotting will take about ',num2str((Obs_time_steps+delta_t/dt)/control_animation_interval*0.15),' seconds. For faster time please choose different control_animation_interval'])
end
% 1000 takes about 150 seconds. 

%% Start plotting the movies and plots
load([movie_name,'.mat'],'time','x','y')
MovieVector=make_movies_plots(N,delta_t,NaN,dt,Obs_time_steps,x,y,NaN,NaN,NaN,NaN,NaN,NaN,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace);
clear time x y
    %% Save movie 
    switch movie_create
        case 'on'
            if control_animation_interval<=Obs_time_steps+dt*delta_t
                frame_rate=10    ;
                if exist('MovieVector')==0
                    load([movie_name,'.mat'],'MovieVector')
                end
                save_movie(MovieVector(2:end),movie_name,frame_rate);
            end
    end
    %% Saving work space
    save([movie_name,'.mat'],'MovieVector','magnify','control_animation_interval','movie_create','ghost','axis_choice','leave_trace','-append')
    clear MovieVector
time_making_movies=toc(making_movies)
%% Start Analyzing the Rotation Behavior
Analyze_rot=tic;
load([movie_name,'.mat'],'time','x','y','v_x','v_y')
v_omega=Rotation(N,x,y,v_0,v_x,v_y,time,delta_t,Obs_time_steps,dt);
save([movie_name,'.mat'],'v_omega','-append')
clear v_omega x y v_x v_y time
time_analyze_rot=toc(Analyze_rot)
%% Clearing Unwanted Timing Variables
clear Analyze_rot combine_data_partitions_start making_movies time_simulation_start
end

