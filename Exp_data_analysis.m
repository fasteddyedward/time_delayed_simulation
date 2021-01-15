%% This file runs analysis and makes movies for the recorded experimental data
file_name_matrix=[
    'dt=0.000.tdms_angle';
    'dt=0.300.tdms_angle';
    'dt=0.600.tdms_angle';
    'dt=0.900.tdms_angle';
    'dt=1.050.tdms_angle';
    'dt=1.200.tdms_angle';
    'dt=1.500.tdms_angle';
    'dt=1.800.tdms_angle';
    ];
for file_name_index=1:size(file_name_matrix,1)
    clearvars -except file_name_index file_name_matrix
%% Read Variables
file_name=file_name_matrix(file_name_index,:);
Table_read = readmatrix([file_name,'.csv']);
frame=Table_read(:,1);
particle=Table_read(:,2);
vel_para=Table_read(:,3);
vel_perp=Table_read(:,4);
op=Table_read(:,5);
angle=Table_read(:,6)*pi/180; % theta was recorded in degrees
cs=Table_read(:,7);
omega=Table_read(:,8); % omega was recorded in rad (confirmed by Xiangzun)
up=Table_read(:,9);
vp=Table_read(:,10);
x(2,:)=Table_read(:,11);
y(2,:)=Table_read(:,12);


%% 2021.1.11 omega?

% omega_self=(x(2,:).*vp'-y(2,:).*up')./sqrt(x(2,:).^2+y(2,:).^2);
% close all
% figure; histogram(omega);
% figure; histogram(omega_self)
%% Other parameters and variables
x(1,1:length(frame))=0;
y(1,1:length(frame))=0;
a=2.19/2 ;
time=1:length(frame);
%% 
% Date='2020.11.27'
% nth_take=1
delta_t_matrix=0
% T_matrix=[1]
% v_0_matrix=1
dt=30*10^-3; 

close all
%% Output File Name
% movie_name=['test3']
movie_name=file_name
%% Setup for Running the program
N=2; % total number of particles in the simulation
delta_t=delta_t_matrix(1); % ms
% Obs_time=Obs_time_steps*dt;
Obs_time_steps=length(frame);
partition_time_steps=Obs_time_steps;
partition_movie='no';
%% State if the particles are fixed, 0 for mobile, 1 for fixed
fixed_flag(1:N)=0;
fixed_flag(1)=1; % particle 1 is fixed

%% Finding axis_scale  =   [movie_x_min  movie_x_max  movie_y_min  movie_y_max]
movie_x_max=0;movie_x_min=0;movie_y_max=0;movie_y_min=0;
movie_x_max=max([movie_x_max max(x(2,:))]);    movie_x_min=min([movie_x_min min(x(2,:))]);    movie_y_max=max([movie_y_max max(y)]);    movie_y_min=min([movie_y_min min(y)]);
movie_x_cent=(movie_x_max+movie_x_min)/2;
movie_y_cent=(movie_y_max+movie_y_min)/2;
plot_width=max(movie_x_max-movie_x_min,movie_y_max-movie_y_min)+4*a;
movie_x_min=movie_x_cent-plot_width/2;
movie_x_max=movie_x_cent+plot_width/2;
movie_y_min=movie_y_cent-plot_width/2;
movie_y_max=movie_y_cent+plot_width/2;
axis_scale=[movie_x_min movie_x_max movie_y_min movie_y_max];
save([movie_name,'.mat'])



%% Making Movies
    %% Parameters for making the movies
%     making_movies=tic;
    magnify=1000    ;
    control_animation_interval=10^3    ; % Record one frame in every ____ frame
    movie_create='off'   ;
    ghost='off'      ;
    force_tracks='off';
    axis_choice='lab'; %'cm' or 'lab'
    leave_trace='off'       ;
    close all

    %% Start Making Movie
    Movie_Vector=Make_Movie(movie_name,N,dt,Obs_time_steps,partition_time_steps,delta_t,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,partition_movie,movie_x_min, movie_x_max,movie_y_min,movie_y_max,a,force_tracks);
    %% Save movie
    frame_rate=10;
    Save_created_movie(movie_name,movie_create,control_animation_interval,Obs_time_steps,delta_t,dt,frame_rate,Movie_Vector)
    %% Saving work space for Movie
%     save([movie_name,'.mat'],'Movie_Vector','magnify','control_animation_interval','movie_create','ghost','axis_choice','leave_trace','-append')
%     clear Movie_Vector
%     time_making_movies=toc(making_movies)
    %% Plotting x and y of the N particles
    figure(2)
    hold on
    plot(time,x(2,:))
    plot(time,y(2,:))
    xlabel('time(frame)')
    legend('x','y')

%% Rotational Analysis: calculates and plots theta
v_0=1
if N==2 && fixed_flag(1)==1
    theta(2,:)=angle;
    save([movie_name,'.mat'],'-append','theta')
    Analyze_theta=tic;
    moving_avg=1 ;
    plot_rot='no';
    Theta_Analysis_Fixed_Center(movie_name,partition_movie,N,v_0,Obs_time_steps,partition_time_steps,delta_t,dt,moving_avg,plot_rot);
    time_analyze_theta=toc(Analyze_theta)
    %% Plotting histogram
    moving_avg=1;
    theta(2,:)=angle;
    save([movie_name,'.mat'],'-append','theta')
    num_bins=100  ;
%     bin_limit=150;
    bin_limit=2;
    figure(81),clf % Note that hist_analysis only works for one particle orbitting a fixed particle at the moment
    
    [k_trans,theta_plus,theta_minus,num_transitions]=hist_analysis(movie_name,moving_avg,num_bins,bin_limit,Obs_time_steps,delta_t,dt);
    set(gca, 'YScale', 'linear')
    saveas(gcf,[movie_name,' (hist).png'])    
    %% Plotting Theta
    moving_avg=1;

    figure(80);clf;
    show_transitions='off';
    plot_theta(N,delta_t,movie_name,moving_avg,theta_plus,theta_minus,k_trans,show_transitions)
%     title(['Theta (Time Delay Angle), v_0 = ',num2str(v_0),', \delta t = ',num2str(delta_t),', T = ',num2str(T)])
    title(file_name)
    saveas(gcf,[movie_name,' (theta).png'])    
    num_transitions
    
end


%% Clearing Unwanted Timing Variables
end
          

