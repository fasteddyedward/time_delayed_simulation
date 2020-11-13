function v_omega=Rotation(N,x,y,v_0,v_x,v_y,time,partition_time_steps,delta_t,partition_movie,moving_avg,plot_rot)
%% Analyzing the angular frequency and plots v_omega
% movie_name
% v_omega(1:N,1:Obs_time_steps+delta_t/dt)=0;
v_omega(1:N,1:partition_time_steps)=0;
circle_center_x(1:3)=0;
circle_center_y(1:3)=0;
time_to_steady_state=600; %ms
%% Calculating Center of Mass
for i=1:N % Center of Circle of each particle's trajectory
    if time(end)>time_to_steady_state
        circle_center_x(i)=mean(x(i,time>time_to_steady_state));
        circle_center_y(i)=mean(y(i,time>time_to_steady_state));
    else
        if i==1
            warning('(Ignore this message first; Rotation.m.) The simulation is not long enough for the system to reach stable state.')
        end
        circle_center_x(i)=mean(mean(x,1),2);
        circle_center_y(i)=mean(mean(y,1),2);
    end
end

%% Calculating v_omega
% for k=1:Obs_time_steps+delta_t/dt % k=1:delta_t/dt does not move
for k=1:partition_time_steps % k=1:delta_t/dt does not move
%     cm_x=mean(x(:,k),1);
%     cm_y=mean(y(:,k),1);
    
    for i=1:N
        R_x=x(i,k)-mean(x(:,k),1);
        R_y=y(i,k)-mean(y(:,k),1);
        %                         R_x=x(i,k)-cm_x;
        %                         R_y=y(i,k)-cm_y;
        %         R_x=x(i,k)-circle_center_x(i);
        %         R_y=y(i,k)-circle_center_y(i);
        R=[R_x R_y 0];
        v=[v_x(i,k) v_y(i,k) 0];
        cross_R_v=cross(R,v);
        v_omega(i,k)=cross_R_v(3)/(norm(R))/v_0; %Normalized angular speed
    end
end

%% Plotting v_omega vs time
switch plot_rot
    case 'yes'
    figure(99)
    hold on
    switch partition_movie
        case 'yes'
            for i=1:N
                plot(time,movmean(v_omega(i,:),moving_avg))
            end
        case 'no'
            for i=1:N
                plot(time,movmean(v_omega(i,:),moving_avg))
            end
    end
    xline(delta_t)
    title('Normalized Rotation Speed (test plot)')
    xlabel('time (ms)')
    ylabel('v_\omega(rad/ms)')
    legend('1','2','3','\deltat')
end
end

