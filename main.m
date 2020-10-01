%% This file runs modulized_time_delay_proto
clear;
close all
%% Setup for Running the program
N=3; % number of particles in the play
delta_t=5; % ms
dt=10^-2; % ms 
Obs_time_steps=10^5
% Obs_time=Obs_time_steps*dt;

%% Coefficients and parameters
v_0= 10^-4; % mm/ms
% gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
T=300; % Kelvin 
% k_B=10^-23; % Boltzmann constant
% D=k_B*T/gamma; % Diffusitivity

%% Initial positions
for i=1:N
    x_init(i)=0;
    y_init(i)=0;
end
for i=1:N
    x_init(i)=i*10^-4;
    y_init(i)=i*10^-4;
end
y_init(3)=2*10^-4;
%% Start calculating finite element numericals for the equation of motion
[x,y,F_x,F_y,v_x,v_y,delta_x,delta_y,time]=modulized_time_delay_proto(N,delta_t,dt,Obs_time_steps,v_0,T,x_init,y_init);
%% Parameters for making the movies
magnify=1000
control_animation_interval=10^2 % Record one frame in every ____ frame
movie_create='on'
ghost='on'
axis_choice='lab'; %'cm' or 'lab'
leave_trace='off'
close all
(Obs_time_steps+delta_t/dt)/control_animation_interval
if ((Obs_time_steps+delta_t/dt)/control_animation_interval>100)
    warning(['The plotting will take about ',num2str((Obs_time_steps+delta_t/dt)/control_animation_interval*0.15),' seconds. For faster time please choose different control_animation_interval'])
end
% 1000 takes about 150 seconds. 

%% Start plotting the movies and plots
tic
[MovieVector,v_omega]=make_movies_plots(N,delta_t,v_0,dt,Obs_time_steps,x,y,F_x,F_y,v_x,v_y,delta_x,delta_y,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace);
toc
%% Save movie 
movie_name='collection'
% movie_name=['delta_t=',num2str(delta_t),', ',axis_choice,' frame, Obs_time_steps=',num2str(Obs_time_steps),', log(dt)=',num2str(log10(dt))]
frame_rate=10
save_movie(MovieVector,movie_name,frame_rate);
%% Saving work space
tic
save('collection')
toc