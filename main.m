%% This file runs modulized_time_delay_proto
%% 2020.10.14 to make the videos with several tries
for nth_take=43
clearvars -except nth_take
close all
%% Output File Name
% movie_name=['2020.10.22,dt=10e-3 take ',num2str(nth_take)];
movie_name=['test3']
warning('Have you modified the file name?')
% pause

%% Setup for Running the program
N=3; % number of particles in the play
delta_t=1; % ms
dt=10^-3; % ms 
% Obs_time=Obs_time_steps*dt;
Obs_time_steps=10^5   ;
partition_time_steps=Obs_time_steps ;
partition_movie='no';
if partition_time_steps>Obs_time_steps
    warning(['Please choose a partitioned time step larger than Obs_time_steps= ',num2str(Obs_time_steps)])
    pause
end
if partition_time_steps<delta_t/dt
    warning(['Please choose a partitioned time step larger than delta/dt= ',num2str(delta_t/dt)])
    pause

end

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

%% Coefficients and parameters
v_0= 1; % mm/ms
% gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
T=0.01; % Kelvin 
% k_B=10^-23; % Boltzmann constant
% D=k_B*T/gamma; % Diffusitivity


%% Initial positions
x_init(1:N)=0;
y_init(1:N)=0;
for i=1:N
    x_init(i)=i*10^1;
    y_init(i)=i*10^1;
end
y_init(3)=2*10^1;
%% Boundary of the particles (for making the movie)
movie_x_max=0;movie_x_min=0;movie_y_max=0;movie_y_min=0;


%% Start calculating finite element numericals for the equation of motion
time_simulation_start=tic;
    %% First stage: Diffusion. t=0 ~ delta_t (lth_partition=1)
    % Obs_time_Steps in this case is round(delta_t/dt)
    ['lth_partition= ',num2str(1),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
    [x,y,v_x,v_y,~]=first_stage_pure_diffusion(N,delta_t,dt,T,x_init,y_init);
    movie_x_max=max([movie_x_max max(x)]);    movie_x_min=min([movie_x_min min(x)]);    movie_y_max=max([movie_y_max max(y)]);    movie_y_min=min([movie_y_min min(y)]);
    save([movie_name,' partition_',num2str(1),'.mat']); %% Takes ~ 0.16 sec
    x_temp=x;
    y_temp=y;
    %% Second stage: Delayed interaction starts. t=delta_t~Obs_time
    for lth_partition=2:round(Obs_time_steps/partition_time_steps)+1
        ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
        [x,y,~,~,v_x,v_y,~,~,~]=second_stage_delayed_int(N,delta_t,dt,partition_time_steps,v_0,T,x_temp,y_temp,lth_partition);
        movie_x_max=max([movie_x_max max(x)]);        movie_x_min=min([movie_x_min min(x)]);        movie_y_max=max([movie_y_max max(y)]);        movie_y_min=min([movie_y_min min(y)]);
        save([movie_name,' partition_',num2str(lth_partition),'.mat']);
        x_temp=x(:,end-delta_t/dt:end); % This last postitions of delayed time delta_t goes to the next round for simulation
        y_temp=y(:,end-delta_t/dt:end);
        clear time
    end
    x_final=x_temp; % We can use this for further simulation
    y_final=y_temp;
    clear x_temp y_temp
    time_simulation=toc(time_simulation_start)

%% Putting all the files into one folder
if (~exist(movie_name, 'dir')); mkdir(movie_name); end%if
for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
    movefile([movie_name,' partition_',num2str(lth_partition),'.mat'],movie_name);
end

switch partition_movie
    case 'no'
        %% Loading all the recorded positions (! Might cause memery overflow)
        combine_data_partitions_start=tic;
        [x,y,v_x,v_y,time]=combine_partitions(movie_name,Obs_time_steps,partition_time_steps,delta_t,dt);
        % movie_x_max=max(x,[],'all');        % movie_y_max=max(y,[],'all');        % movie_x_min=min(x,[],'all');        % movie_y_min=min(y,[],'all');
        save([movie_name,'.mat']);
        clear x y v_x v_y time
        time_combine_data_partitions=toc(combine_data_partitions_start)
end

%% Finding movie_x_min
if (~exist('movie_x_max'))
% if (exist(movie_x_max)==0)
    movie_x_max=0;movie_x_min=0;movie_y_max=0;movie_y_min=0;
    cd(movie_name)
    for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
        load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y')
        movie_x_max=max([movie_x_max max(x)]);    movie_x_min=min([movie_x_min min(x)]);    movie_y_max=max([movie_y_max max(y)]);    movie_y_min=min([movie_y_min min(y)]);
    end
    cd ..
end
axis_scale=[movie_x_min movie_x_max movie_y_min movie_y_max];

%% Making Movies and Analses
    %% Parameters for making the movies
    making_movies=tic;
    magnify=1000    ;
    control_animation_interval=10^2     ; % Record one frame in every ____ frame
    movie_create='on'   ;
    ghost='off'      ;
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

%% Start Making Movie and Analyzing v_omega
switch partition_movie
    case 'no'
        %% Start plotting the movies and plots
        load([movie_name,'.mat'],'time','x','y')
        Movie_Vector=make_movies_plots(N,delta_t,NaN,dt,Obs_time_steps+delta_t/dt,Obs_time_steps,x,y,NaN,NaN,NaN,NaN,NaN,NaN,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale);
        clear time x y
        %% Save movie
        switch movie_create
            case 'on'
                if control_animation_interval<=Obs_time_steps+dt*delta_t
                    frame_rate=10    ;
                    if exist('Movie_Vector')==0
                        load([movie_name,'.mat'],'Movie_Vector')
                    end
                    save_movie(Movie_Vector(2:end),movie_name,frame_rate);
                end
        end
        %% Saving work space for Movie
        save([movie_name,'.mat'],'Movie_Vector','magnify','control_animation_interval','movie_create','ghost','axis_choice','leave_trace','-append')
        clear Movie_Vector
        time_making_movies=toc(making_movies)
        %% Start Analyzing the Rotation Behavior
        Analyze_rot=tic;
        load([movie_name,'.mat'],'time','x','y','v_x','v_y')
        v_omega=Rotation(N,x,y,v_0,v_x,v_y,time,Obs_time_steps+delta_t/dt,delta_t,partition_movie);
        save([movie_name,'.mat'],'v_omega','-append')
        clear v_omega x y v_x v_y time
        time_analyze_rot=toc(Analyze_rot)

    case 'yes'
        %% Start plotting the movies and plots + Analyze v_omega
        close all
        Movie_Vector=[];
        v_omega=[];
        axis_scale=[movie_x_min movie_x_max movie_y_min movie_y_max];
        for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
            if lth_partition==1
                %% Start plotting movie for stage 1
                cd(movie_name)
                load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
                cd ..
                time=(1:delta_t/dt)*dt;
                Movie_Vector_Partition=make_movies_plots(N,delta_t,NaN,dt,delta_t/dt,NaN,x,y,NaN,NaN,NaN,NaN,NaN,NaN,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale);
                v_omega_partition=Rotation(N,x,y,v_0,v_x,v_y,time,delta_t/dt,delta_t,partition_movie);
                clear x y v_x v_y time
            else
                %% Start plotting movie for stage 2
                cd(movie_name)
                ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
                load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
                cd ..
                time=(1+(lth_partition-2)*partition_time_steps+delta_t/dt:(lth_partition-1)*partition_time_steps+delta_t/dt)*dt;
                Movie_Vector_Partition=make_movies_plots(N,delta_t,NaN,dt,partition_time_steps,NaN,x,y,NaN,NaN,NaN,NaN,NaN,NaN,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale);
                v_omega_partition=Rotation(N,x,y,v_0,v_x,v_y,time,partition_time_steps,delta_t,partition_movie);
                clear x y v_x v_y time
            end
            Movie_Vector=[Movie_Vector Movie_Vector_Partition];
            v_omega=[v_omega v_omega_partition];
%             v_omega=0;
        end
        %% Save movie
        switch movie_create
            case 'on'
                if control_animation_interval<=Obs_time_steps+dt*delta_t
                    frame_rate=10    ;
                    if exist('Movie_Vector')==0
                        load([movie_name,'.mat'],'Movie_Vector')
                    end
                    save_movie(Movie_Vector(2:end),movie_name,frame_rate);
                end
        end
        %% Saving work space for Movie and analzye v_omega
        save([movie_name,'.mat'],'v_omega','Movie_Vector','magnify','control_animation_interval','movie_create','ghost','axis_choice','leave_trace','-append')
        clear Movie_Vector v_omega
        time_making_movies=toc(making_movies)
        %% Save work space for analyzing v_omega
%         Analyze_rot=tic;
%         load([movie_name,'.mat'],'time','x','y','v_x','v_y')
%         v_omega=Rotation(N,x,y,v_0,v_x,v_y,time,Obs_time_steps+delta_t/dt,delta_t,partition_movie);
%         save([movie_name,'.mat'],'v_omega','-append')
%         clear v_omega x y v_x v_y time
%         time_analyze_rot=toc(Analyze_rot)
end


%% Clearing Unwanted Timing Variables
clear Analyze_rot combine_data_partitions_start making_movies time_simulation_start
end

