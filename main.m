%% This file runs to check Validity of Pure Diffusion
for nth_take=26
clearvars -except nth_take
close all
%% Output File Name
% warning('Have you modified the file name?')
% pause
% movie_name=['2020.10.16,dt=10e-3 take ',num2str(nth_take)];
movie_name=['2020.10.16, pure diffusion test.']
% answer = questdlg(['Is the file name correct? ',movie_name], ...
% 	'File Name Check', ...
% 	'Yes','No, I will  modify it','No, I will modify it ');
% switch answer
%     case 'Yes'
%         'Start running program...'
%     case 'No, I will  modify it'
%         error('Please modify the file name')
% end
tic
%% Setup for Running the program
N=20; % number of particles in the play
delta_t=0; % ms
dt=10^-5; % ms 
Obs_time_steps=10^6    ;
% Obs_time=Obs_time_steps*dt;

%% Coefficients and parameters
v_0= 0; % mm/ms
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
[x,y,~,~,v_x,v_y,~,~,time]=modulized_time_delay_proto(N,delta_t,dt,Obs_time_steps,v_0,T,x_init,y_init);
save([movie_name,'.mat']);
toc
tic
clear time
%% Parameters for making the movies
% magnify=1000    ;
% control_animation_interval=10^4     ; % Record one frame in every ____ frame
% movie_create='off'   ;
% ghost='on'      ;
% axis_choice='lab'; %'cm' or 'lab'
% leave_trace='off'       ;
% close all
% (Obs_time_steps+delta_t/dt)/control_animation_interval      ;
% if ((Obs_time_steps+delta_t/dt)/control_animation_interval>100)
%     warning(['The plotting will take about ',num2str((Obs_time_steps+delta_t/dt)/control_animation_interval*0.15),' seconds. For faster time please choose different control_animation_interval'])
% end
% % 1000 takes about 150 seconds. 

%% Start plotting the movies and plots
% load([movie_name,'.mat'],'time','x','y','v_x','v_y')
% [MovieVector,v_omega]=make_movies_plots(N,delta_t,v_0,dt,Obs_time_steps,x,y,NaN,NaN,v_x,v_y,NaN,NaN,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace);
% clear x y v_x v_y time
% %% Save movie 
% switch movie_create
% case 'on'
%     frame_rate=5    ;
%     save_movie(MovieVector(2:end),movie_name,frame_rate);
% end
%% Saving work space
% save([movie_name,'.mat'],'MovieVector','v_omega','magnify','control_animation_interval','movie_create','ghost','axis_choice','leave_trace','-append')
% toc
% %%
% clear MovieVector v_omega
end



%% Checking x_rms vs time : x_rms ~ sqrt(6 D t)?
%% Here is the for loop to include do something repeatedly with the original codes
load([movie_name,'.mat'],'x','y','time')
gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
T=300; % Kelvin
k_B=10^-23; % Boltzmann constant
D=k_B*T/gamma; % Diffusitivity
close all

%% Start calculating the rms of the x: <(x(t)-x(0))^2> vs t
% x_rms=[];
x_rms_square(1:Obs_time_steps+delta_t/dt)=0;
for k=1:Obs_time_steps+delta_t/dt
    for i=1:N
        x_rms_square(k)=x_rms_square(k)+(x(i,k)-x(i,1))^2+(y(i,k)-y(i,1))^2;
    end
end

%% Plotting the result
clf
x_rms_square=(x_rms_square/N);
x_rms=sqrt(x_rms_square);
hold on
plot(time,x_rms_square)
plot(time,(6*D*time))
title('<x_{rms}^2> vs t')

figure
hold on
plot(time,x_rms)
plot(time,sqrt(6*D*time))
title('<x_{rms}> vs t')

