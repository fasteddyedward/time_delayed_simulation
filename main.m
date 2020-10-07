%% This file runs modulized_time_delay_proto
%% 2020.10.6-7 This is a branch for testing deterministic motion of particles (T=0)
%% Also note that the v_0 has been modified to 10^-6 mm/ms
clear;
close all
%% Here are the control parameters set to test dt in this file (the original values will be replaced with variable_name_test)
dt_test=[10^0, 10^-1, 10^-2, 10^-3]
Obs_time_test=10^3
for nth_loop=1:length(dt_test)
    %% Setup for Running the program
    N=3; % number of particles in the play
    delta_t=50; % ms
    %     dt=10^-2; % ms
    dt=dt_test(nth_loop)
    %     Obs_time_steps=10^5
    Obs_time_steps=Obs_time_test./dt
    % Obs_time=Obs_time_steps*dt;
    
    %% Coefficients and parameters
    v_0= 10^-6; % mm/ms
    % gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
    T=0; % Kelvin
    % k_B=10^-23; % Boltzmann constant
    % D=k_B*T/gamma; % Diffusitivity
    
    %% Initial positions
    for i=1:N
        x_init(i)=0;
        y_init(i)=0;
    end
    for i=1:N
        x_init(i)=i*10^-6;
        y_init(i)=i*10^-6;
    end
    y_init(3)=1*10^-6;
    %% Start calculating finite element numericals for the equation of motion
    [x,y,F_x,F_y,v_x,v_y,delta_x,delta_y,time]=modulized_time_delay_proto(N,delta_t,dt,Obs_time_steps,v_0,T,x_init,y_init);
    %% Parameters for making the movies
    magnify=1000;
    control_animation_interval=10^2; % Record one frame in every ____ frame
    movie_create='off';
    ghost='on';
    axis_choice='lab'; %'cm' or 'lab'
    leave_trace='off';
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
    switch movie_create
        case 'on'
            %             movie_name='collection';
            movie_name=['2020.10.7_dt=',num2str(nth_loop),'.avi'];
            % movie_name=['delta_t=',num2str(delta_t),', ',axis_choice,' frame, Obs_time_steps=',num2str(Obs_time_steps),', log(dt)=',num2str(log10(dt))]
            frame_rate=10;
%             cd videos_dt_test_T=0
            save_movie(MovieVector,movie_name,frame_rate);
%             cd ..
    end
    %% Saving work space
    tic
    cd videos_dt_test_T=0
    %     save('collection')
    save(['2020.10.7_dt=',num2str(nth_loop),'.mat'])
    cd ..
    toc
    
end