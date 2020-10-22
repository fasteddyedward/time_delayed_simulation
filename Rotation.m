function v_omega=Rotation(N,x,y,v_0,v_x,v_y,time,partition_time_steps,delta_t,partition_movie)
%% Analyzing the angular frequency
% movie_name
% v_omega(1:N,1:Obs_time_steps+delta_t/dt)=0;
v_omega(1:N,1:partition_time_steps)=0;
circle_center_x(1:3)=0;
circle_center_y(1:3)=0;
time_to_steady_state=600; %ms
for i=1:N % Center of Circle of each particle's trajectory
    if time(end)>time_to_steady_state
        circle_center_x(i)=mean(x(i,time>time_to_steady_state));
        circle_center_y(i)=mean(y(i,time>time_to_steady_state));
    else
        if i==1
            warning('The simulation is not long enough for the system to reach stable state.')
        end
        circle_center_x(i)=mean(mean(x,1),2);
        circle_center_y(i)=mean(mean(y,1),2);
    end
end

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
figure(99)
hold on
switch partition_movie
    case 'yes'
plot(time,v_omega(1,:))
plot(time,v_omega(2,:))
plot(time,v_omega(3,:))
    case 'no'
plot(time,movmean(v_omega(1,:),100000))
plot(time,movmean(v_omega(2,:),100000))
plot(time,movmean(v_omega(3,:),100000))
end
xline(delta_t)
title('Normalized Rotation Speed')
xlabel('time (ms)')
ylabel('v_\omega(rad/ms) /v_0 ')
legend('1','2','3','\deltat')        
end

