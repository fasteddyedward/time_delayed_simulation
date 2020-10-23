%% This file runs modulized_time_delay_proto
%% 2020.10.14 to make the videos with several tries
for nth_take=47
clearvars -except nth_take
close all
%% Output File Name
movie_name=['2020.10.23,dt=10e-3 take ',num2str(nth_take)];
% movie_name=['test3']
warning('Have you modified the file name?')

%% Setup for Running the program
N=3; % number of particles in the play
delta_t=1; % ms
dt=10^-3; % ms 
% Obs_time=Obs_time_steps*dt;
Obs_time_steps=10^5
partition_time_steps=10^4
partition_movie='no'
if partition_time_steps>Obs_time_steps
    warning(['Please choose a partitioned time step larger than Obs_time_steps= ',num2str(Obs_time_steps)])
    pause
elseif partition_time_steps<delta_t/dt
    warning(['Please choose a partitioned time step larger than delta/dt= ',num2str(delta_t/dt)])
    pause
end
%% Coefficients and parameters
v_0= 1; % mm/ms
% gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
T=0.01; % Kelvin 
% k_B=10^-23; % Boltzmann constant
gamma=1;
k_B=1;
D=k_B*T/gamma; % Diffusitivity



%% Initial positions
x_init(1:N)=0;
y_init(1:N)=0;
for i=1:N
    x_init(i)=i*10^1;
    y_init(i)=i*10^1;
end
y_init(3)=2*10^1;
%% Boundary of the particles (for making the movie)
% movie_x_max=0;movie_x_min=0;movie_y_max=0;movie_y_min=0;


%% Start calculating finite element numericals for the equation of motion
time_simulation_start=tic;
[x,y,v_x,v_y,movie_x_min,movie_x_max,movie_y_min,movie_y_max,x_final,y_final]=simulation(movie_name,Obs_time_steps,partition_time_steps,0,0,0,0,N,delta_t,dt,v_0,T,gamma,k_B,D,x_init,y_init);
time_simulation=toc(time_simulation_start)

%% Putting all the files into one folder
if (~exist(movie_name, 'dir')); mkdir(movie_name); end%if
for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
    movefile([movie_name,' partition_',num2str(lth_partition),'.mat'],movie_name);
end
%% Loading all the recorded positions (! Might cause memery overflow)
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
if (~exist('movie_x_max'))
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

    %% Start Making Movie and Analyzing v_omega
    switch partition_movie
        case 'no'
            %% Start plotting the movies and plots
            load([movie_name,'.mat'],'time','x','y')
            Movie_Vector=make_movies_plots(N,delta_t,dt,Obs_time_steps+delta_t/dt,x,y,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale);
            clear time x y
            %% Save movie
            switch movie_create
                case 'on'
                    if control_animation_interval<=Obs_time_steps+dt*delta_t
                        frame_rate=10    ;
                        if exist('Movie_Vector','var')==0
                            load([movie_name,'.mat'],'Movie_Vector')
                        end
                        save_movie(Movie_Vector(2:end),movie_name,frame_rate);
                    end
            end
            %% Saving work space for Movie
            save([movie_name,'.mat'],'Movie_Vector','magnify','control_animation_interval','movie_create','ghost','axis_choice','leave_trace','-append')
            clear Movie_Vector
            time_making_movies=toc(making_movies)

        case 'yes'
            %% Start plotting the movies and plots
            close all
            Movie_Vector=[];
%             v_omega=[];
            axis_scale=[movie_x_min movie_x_max movie_y_min movie_y_max];
            for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
                if lth_partition==1
                    %% Start plotting movie for stage 1
                    cd(movie_name)
                    load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
                    cd ..
                    time=(1:delta_t/dt)*dt;
                    Movie_Vector_Partition=make_movies_plots(N,delta_t,dt,delta_t/dt,x,y,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale);
                    v_omega_partition=Rotation(N,x,y,v_0,v_x,v_y,time,delta_t/dt,delta_t,partition_movie);
                    clear x y v_x v_y time
                else
                    %% Start plotting movie for stage 2
                    cd(movie_name)
                    ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
                    load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
                    cd ..
                    time=(1+(lth_partition-2)*partition_time_steps+delta_t/dt:(lth_partition-1)*partition_time_steps+delta_t/dt)*dt;
                    Movie_Vector_Partition=make_movies_plots(N,delta_t,dt,partition_time_steps,x,y,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale);
            clear x y v_x v_y time
                end
                Movie_Vector=[Movie_Vector Movie_Vector_Partition];
            end
            %% Save movie
            switch movie_create
                case 'on'
                    if control_animation_interval<=Obs_time_steps+dt*delta_t
                        frame_rate=10    ;
                        if exist('Movie_Vector','var')==0
                            load([movie_name,'.mat'],'Movie_Vector')
                        end
                        save_movie(Movie_Vector(2:end),movie_name,frame_rate);
                    end
            end
            %% Saving work space for Movie
            save([movie_name,'.mat'],'Movie_Vector','magnify','control_animation_interval','movie_create','ghost','axis_choice','leave_trace','-append')
            clear Movie_Vector
            time_making_movies=toc(making_movies)
    end
    
    
    
    
%% Rotational Analysis
Analyze_rot=tic;
switch partition_movie
    case 'no'
        load([movie_name,'.mat'],'time','x','y','v_x','v_y')
        v_omega=Rotation(N,x,y,v_0,v_x,v_y,time,Obs_time_steps+delta_t/dt,delta_t,partition_movie);
        %% Saving v_omega
        save([movie_name,'.mat'],'v_omega','-append')
        clear v_omega x y v_x v_y time
    case 'yes'
        v_omega=[];
        for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
            if lth_partition==1
                %% Start plotting v_omega for stage 1
                cd(movie_name)
                load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
                cd ..
                time=(1:delta_t/dt)*dt;
                v_omega_partition=Rotation(N,x,y,v_0,v_x,v_y,time,delta_t/dt,delta_t,partition_movie);
                clear x y v_x v_y time
            else
                %% Start plotting movie for stage 2
                cd(movie_name)
                ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
                load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
                cd ..
                time=(1+(lth_partition-2)*partition_time_steps+delta_t/dt:(lth_partition-1)*partition_time_steps+delta_t/dt)*dt;
                v_omega_partition=Rotation(N,x,y,v_0,v_x,v_y,time,partition_time_steps,delta_t,partition_movie);
                clear x y v_x v_y time
            end
            % Note appending v_omega like this might cause memory overflow.
            v_omega=[v_omega v_omega_partition];
            %v_omega=0;
        end
        %% Saving v_omega
        save([movie_name,'.mat'],'v_omega','-append')
        clear v_omega
end
time_analyze_rot=toc(Analyze_rot)
%% Drawing v_omega
load([movie_name,'.mat'],'v_omega','time')
figure(98) %% Would be same as figure(99) if partition movie='no'
hold on
for i=1:N
    plot(time,movmean(v_omega(i,:),100000))

end
xline(delta_t)
title('Normalized Rotation Speed')
xlabel('time (ms)')
ylabel('v_\omega(rad/ms) /v_0 ')
legend('1','2','3','\deltat')        


%% Clearing Unwanted Timing Variables
clear Analyze_rot combine_data_partitions_start making_movies time_simulation_start
end

