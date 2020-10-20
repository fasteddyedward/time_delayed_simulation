%% This file runs modulized_time_delay_proto
%% 2020.10.14 to make the videos with several tries
for nth_take=34
clearvars -except nth_take
close all
%% Output File Name
movie_name=['2020.10.20,dt=10e-4 take ',num2str(nth_take)];
warning('Have you modified the file name?')
% pause
tic
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
delta_t=50; % ms
dt=10^-3; % ms 
Obs_time_steps=10^6   ;
% Obs_time=Obs_time_steps*dt;
partition_time_steps=10^5   ;
if partition_time_steps<delta_t/dt
    warning(['Please choose a partitioned time step larger than delta/dt= ',num2str(delta_t/dt)])
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
movie_x_max=0;movie_x_min=0;movie_y_max=0;movie_y_min=0;
%% Start calculating finite element numericals for the equation of motion

%% First stage: Diffusion. t=0 ~ delta_t (lth_partition=1)
% Obs_time_Steps in this case is round(delta_t/dt)
[x,y,v_x,v_y,time]=pure_diffusion(N,delta_t,dt,T,x_init,y_init);
movie_x_max=max([movie_x_max max(x)]);    movie_x_min=min([movie_x_min min(x)]);    movie_y_max=max([movie_y_max max(y)]);    movie_y_min=min([movie_y_min min(y)]);
save([movie_name,'partition_',num2str(1),'.mat']); %% Takes ~ 0.16 sec
x_temp=x;
y_temp=y;
%% Second stage: Delayed interaction starts. t=delta_t~Obs_time
for lth_partition=2:round(Obs_time_steps/partition_time_steps)+1
    [x,y,~,~,v_x,v_y,~,~,time]=modulized_time_delay_proto(N,delta_t,dt,partition_time_steps,v_0,T,x_temp,y_temp,lth_partition);
    movie_x_max=max([movie_x_max max(x)]);        movie_x_min=min([movie_x_min min(x)]);        movie_y_max=max([movie_y_max max(y)]);        movie_y_min=min([movie_y_min min(y)]);
    save([movie_name,'partition_',num2str(lth_partition),'.mat']);
    x_temp=x(:,end-delta_t/dt:end); % This last postitions of delayed time delta_t goes to the next round for simulation
    y_temp=y(:,end-delta_t/dt:end);
    clear time
end
toc





tic
%% Parameters for making the movies
magnify=1000    ;
control_animation_interval=10^3     ; % Record one frame in every ____ frame
movie_create='on'   ;
ghost='off'      ;  % Currently I'm too lazy to fix this feature. Sorry.
axis_choice='lab'; %'cm' or 'lab' 
axis_scale=[movie_x_min movie_x_max movie_y_min movie_y_max];
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

%% Start making movie
close all
Movie_Vector=[];
for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
    if lth_partition==1
        %% Start plotting the movies and plots
        load([movie_name,'partition_',num2str(lth_partition),'.mat'],'time','x','y','v_x','v_y')
        [Movie_Vector_Partition,v_omega]=make_movies_plots(N,delta_t,v_0,dt,delta_t/dt,x,y,NaN,NaN,v_x,v_y,NaN,NaN,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale);
        clear x y v_x v_y time
    else
        %% Start plotting the movies and plots
        load([movie_name,'partition_',num2str(lth_partition),'.mat'],'time','x','y','v_x','v_y')
        [Movie_Vector_Partition,v_omega]=make_movies_plots(N,delta_t,v_0,dt,partition_time_steps,x,y,NaN,NaN,v_x,v_y,NaN,NaN,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale);
        clear x y v_x v_y time
    end
    Movie_Vector=[Movie_Vector Movie_Vector_Partition];
end
%% Save movie
switch movie_create
    case 'on'
        if control_animation_interval<=Obs_time_steps+dt*delta_t
            frame_rate=30    ;
            save_movie(Movie_Vector(2:end),movie_name,frame_rate);
        end
end
toc
%% Saving work space for Movie (can be used to resume last section Save Movie 
save([movie_name,'_movie.mat'],'Movie_Vector','control_animation_interval','Obs_time_steps','dt','delta_t','movie_create')
toc
%%
clear Movie_Vector

%% Loading all the recorded positions (! Might cause memery overflow)
tic
x_full_time=[];
y_full_time=[];
for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
    if lth_partition==1
        %% Start plotting the movies and plots
        load([movie_name,'partition_',num2str(lth_partition),'.mat'],'time','x','y','v_x','v_y')
        x_full_time=[x_full_time x];
        y_full_time=[y_full_time y];
    else
        %% Start plotting the movies and plots
        load([movie_name,'partition_',num2str(lth_partition),'.mat'],'time','x','y','v_x','v_y')
        x_full_time=[x_full_time x];
        y_full_time=[y_full_time y];
    end
end
toc
end





