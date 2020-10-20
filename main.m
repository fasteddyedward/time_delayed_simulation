%% This file runs modulized_time_delay_proto
%% 2020.10.14 to make the videos with several tries
for nth_take=30
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
delta_t=5; % ms
dt=10^-3; % ms 
Obs_time_steps=10^5   ;
% Obs_time=Obs_time_steps*dt;
partition_time_steps=10^4   ;
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
x_init(1:N)=[];
y_init(1:N)=[];
% for i=1:N
%     x_init(i)=0;
%     y_init(i)=0;
% end
for i=1:N
    x_init(i)=i*10^-4;
    y_init(i)=i*10^-4;
end
y_init(3)=2*10^-4;

%% Start calculating finite element numericals for the equation of motion

%% First stage: Diffusion. t=0 ~ delta_t (lth_partition=1)
% Obs_time_Steps in this case is round(delta_t/dt)
% [x,y,~,~,v_x,v_y,~,~,time]=pure_diffusion(N,delta_t,dt,delta_t/dt,v_0,T,x_init,y_init);
[x,y,v_x,v_y,time]=pure_diffusion(N,delta_t,dt,T,x_init,y_init);
save([movie_name,'partition_',num2str(1),'.mat']); %% Takes ~ 0.16 sec
x_temp=x;
y_temp=y;
%% Second stage: Delayed interaction starts. t=delta_t~Obs_time
for lth_partition=2:round(Obs_time_steps/partition_time_steps)+1
    [x,y,~,~,v_x,v_y,~,~,time]=modulized_time_delay_proto(N,delta_t,dt,partition_time_steps,v_0,T,x_temp,y_temp,lth_partition);
    save([movie_name,'partition_',num2str(lth_partition),'.mat']);
    x_temp=x(:,end-delta_t/dt:end); % Maybe end-delta_t/dt+1:end
    y_temp=y(:,end-delta_t/dt:end);
end
clear time
toc
tic
%% Parameters for making the movies
magnify=1000    ;
control_animation_interval=10^3     ; % Record one frame in every ____ frame
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

%%
close all
for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
    if lth_partition==1
        %% Start plotting the movies and plots
        % load([movie_name,'.mat'],'time','x','y','v_x','v_y')
        load([movie_name,'partition_',num2str(lth_partition),'.mat'],'time','x','y','v_x','v_y')
        [MovieVector,v_omega]=make_movies_plots(N,delta_t,v_0,dt,delta_t/dt,x,y,NaN,NaN,v_x,v_y,NaN,NaN,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace);
        clear x y v_x v_y time
    else
        %% Start plotting the movies and plots
        % load([movie_name,'.mat'],'time','x','y','v_x','v_y')
        load([movie_name,'partition_',num2str(lth_partition),'.mat'],'time','x','y','v_x','v_y')
        [MovieVector,v_omega]=make_movies_plots(N,delta_t,v_0,dt,partition_time_steps,x,y,NaN,NaN,v_x,v_y,NaN,NaN,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace);
        clear x y v_x v_y time
    end
    %% Save movie
    switch movie_create
        case 'on'
            if control_animation_interval<=Obs_time_steps+dt*delta_t
                frame_rate=5    ;
                save_movie(MovieVector(2:end),movie_name,frame_rate);
            end
    end
end
%% Saving work space
save([movie_name,'.mat'],'MovieVector','v_omega','magnify','control_animation_interval','movie_create','ghost','axis_choice','leave_trace','-append')
toc
%%
clear MovieVector v_omega
end
