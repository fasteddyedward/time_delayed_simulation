%% This file runs modulized_time_delay_proto
%% 2020.10.5 This is a branch for testing the valid dt for pure diffusion (v_0=0)
clear;
close all
%% Setup of the testing the pure diffusion with different dt
dt_test=10^-3
% dt_test=10^-4 % Beware of data leakage! About 4.5 GB of ram will be taken
% use of.
% dt_test=10^-3 with N=100 and obs_time=10^3 will crash the computer!
%%
for dt_index=1:size(dt_test,2)
    
%% Setup for Running the program
N=100; % number of particles in the play
delta_t=0; % ms
dt=dt_test(dt_index); % ms 
% Obs_tme

Obs_time=10^2
Obs_time_steps=Obs_time./dt
% Obs_time=Obs_time_steps*dt;

%% Coefficients and parameters
v_0=0; % mm/ms
% gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
T=300; % Kelvin 
% k_B=10^-23; % Boltzmann constant
% D=k_B*T/gamma; % Diffusitivity
%% Initial positions
for i=1:N
    x_init(i)=0;
    y_init(i)=0;
end
% for i=1:N
%     x_init(i)=i*10^-4;
%     y_init(i)=i*10^-4;
% end
% y_init(3)=2*10^-4;
%% Start calculating finite element numericals for the equation of motion
[x,y,F_x,F_y,v_x,v_y,delta_x,delta_y,time]=modulized_time_delay_proto(N,delta_t,dt,Obs_time_steps,v_0,T,x_init,y_init);
%% Saving work space of solution to equation of motion
tic
cd 'Matfiles'
save(['N=',num2str(N),', v_0=',num2str(v_0),' ,T=',num2str(T),', delta_t=',num2str(delta_t),', dt=',num2str(dt),'.mat'])
cd ..
toc
end
% pause

%% Checking x_rms vs time : x_rms ~ sqrt(6 D t)?
%% Here is the for loop to include do something repeatedly with the original codes
N_draw=N; % number of particles in the play
delta_t=0; % ms
dt_test_draw=dt_test
% dt_test_draw=[10^1 10^0 10^-1 10^-2]
Obs_time_draw=Obs_time
% Obs_time_draw=10^2
Obs_time_steps=Obs_time_draw./dt_test_draw
gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
T=300; % Kelvin 
k_B=10^-23; % Boltzmann constant
D=k_B*T/gamma; % Diffusitivity
v_0=0
for dt_index=1:size(dt_test_draw,2)
    close all
    dt=dt_test_draw(dt_index)
    cd 'Matfiles'
    load(['N=',num2str(N_draw),', v_0=',num2str(v_0),' ,T=',num2str(T),', delta_t=',num2str(delta_t),', dt=',num2str(dt),'.mat'])
    'loading success'
    cd ..
    
    %% Start calculating the rms of the x: <(x(t)-x(0))^2> vs t
    x_rms=[];
    x_rms(1:Obs_time_steps+delta_t/dt)=0;
    for k=1:Obs_time_steps+delta_t/dt
        for i=1:N
            x_rms(k)=x_rms(k)+(x(i,k)-x(i,1))^2+(y(i,k)-y(i,1))^2;
        end
    end
    x_rms=sqrt(x_rms/N);
%     time=dt*(1:Obs_time_steps(dt_index)+delta_t/dt)
    plot(time,x_rms)
    hold on
    plot(time,sqrt(6*D*time))
    pause
end

%% Parameters for making the movies
% magnify=1000
% control_animation_interval=10^3 % Record one frame in every ____ frame
% movie_create='on'
% ghost='off'
% axis_choice='lab'; %'cm' or 'lab'
% leave_trace='off'
% close all
% % (Obs_time_steps+delta_t/dt)/control_animation_interval
% ['The plotting will take about ',num2str((Obs_time_steps+delta_t/dt)/control_animation_interval*0.15),' seconds. For faster time please choose different control_animation_interval']
% if ((Obs_time_steps+delta_t/dt)/control_animation_interval>100)
%     warning(['The plotting will take about ',num2str((Obs_time_steps+delta_t/dt)/control_animation_interval*0.15),' seconds. For faster time please choose different control_animation_interval'])
% end
% % 1000 takes about 150 seconds. 
% %% Start plotting the movies and plots
% tic
% [MovieVector,v_omega]=make_movies_plots(N,delta_t,v_0,dt,Obs_time_steps,x,y,F_x,F_y,v_x,v_y,delta_x,delta_y,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace);
% toc
% %% Save movie 
% % movie_name='collection'
% movie_name=['N=',num2str(N),', v_0=',num2str(v_0),' ,T=',num2str(T),', delta_t=',num2str(delta_t),', dt=',num2str(dt)];
% % movie_name=['delta_t=',num2str(delta_t),', ',axis_choice,' frame, Obs_time_steps=',num2str(Obs_time_steps),', log(dt)=',num2str(log10(dt))]
% frame_rate=10
% save_movie(MovieVector,movie_name,frame_rate);
% %% Saving work space of movie frames + solution to equation of motion
% tic
% cd 'Matfiles'
% save(['N=',num2str(N),', v_0=',num2str(v_0),' ,T=',num2str(T),', delta_t=',num2str(delta_t),', dt=',num2str(dt),'.mat'])
% cd ..
% toc
