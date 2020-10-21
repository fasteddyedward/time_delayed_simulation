function [MovieVector]=make_movies_plots(N,delta_t,~,dt,Obs_time_steps,x,y,~,~,~,~,~,~,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace)
% function [MovieVector,v_omega]=make_movies_plots(N,delta_t,v_0,dt,Obs_time_steps,x,y,~,~,v_x,v_y,~,~,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace)
switch movie_create
    case 'off'
    MovieVector=[];
end
%% Plotting figures
    %% Plotting v_x
% figure
% hold on
% for i=1:N
%     plot(time,v_x(i,:))
% end
% title('v_x')
% legend('1','2','3')

    %% Plotting v_y
% figure
% hold on
% for i=1:N
%     plot(time,v_y(i,:))
% end
% title('v_y')
% legend('1','2','3')
%%
for i=1:N
    figure
    hold on
    plot(time,x(i,1:end-1))
    plot(time,y(i,1:end-1))
    title(['Position of particle ',num2str(i)])
    legend('x','y')
    xlabel('time (ms)')
    ylabel('position (mm)')
end
%% Making movie of the partilces
tic
switch movie_create
    case 'off'
        'No movie created.'
    case 'on'      
    animation_interval=control_animation_interval;
    % figure('visible','off'); % turning off the figure actually increases the
    % time

    movie_frame_index=1; % For movie recording, increases 1 for each recorded frame
    mean_x=mean(x(:,1+2*delta_t/dt:end),2); % We take the part of x and y after t=delta_t 
    mean_y=mean(y(:,1+2*delta_t/dt:end),2);
    std=0; % This is for setting the axis
    figure
    for i=1:N
        for j=i+1:N
            std=std+(mean_x(i)-mean_x(j))^2+(mean_y(i)-mean_y(j))^2;
        end
    end
    std=sqrt(std/(N-1));
    
    %     for k=1+delta_t/dt:animation_interval:Obs_time_steps+delta_t/dt % frames with k=1:delta_t/dt do not move
    %         MovieVector(1:((Obs_time_steps+delta_t/dt)/animation_interval))=0;
    for k=1:animation_interval:Obs_time_steps+delta_t/dt % frames with k=1:delta_t/dt only diffuses
        switch leave_trace
            case 'on'
            case 'off'
                clf
        end
        grid on
        %% Calculating center of mass and the std deviation of the particles
%         cm_x=mean(x(:,k),1);
%         cm_y=mean(y(:,k),1);

        for i=1:N
            hold on
            scatter(x(i,k),y(i,k),'filled')
            switch ghost
                case 'on'
                    if k>delta_t/dt
                        scatter(x(i,k-delta_t/dt),y(i,k-delta_t/dt),'o')
                    end
                case 'off'
            end
        end
        %         plot(cm_x,cm_y,'.')
        plot(mean(x(:,k),1),mean(y(:,k),1),'.')
        %% Don't use legend since it costs too much time
%         legend('1','2','3','CM','Location','northeastoutside')
        title(['Movie with ',axis_choice,' frame'])

        %% Axis of choice
        switch axis_choice
            case 'lab'
                axis_scale=[min(min(x)) max(max(x)) min(min(y)) max(max(y))];
                axis(axis_scale);
            case 'cm'
                %                                 axis([cm_x-magnify*std   cm_x+magnify*std   cm_y-magnify*std   cm_y+magnify*std]);
                %                 axis([cm_x-magnify*std   cm_x+magnify*std   cm_y-magnify*std   cm_y+magnify*std]);
                axis([mean(x(:,k),1)-magnify*std   mean(x(:,k),1)+magnify*std   mean(y(:,k),1)-magnify*std   mean(y(:,k),1)+magnify*std]);
                axis([mean(x(:,k),1)-magnify*std   mean(x(:,k),1)+magnify*std   mean(y(:,k),1)-magnify*std   mean(y(:,k),1)+magnify*std]);
                
        end
        
        %%% For making movie
        MovieVector(movie_frame_index)=getframe(gcf);
        movie_frame_index=movie_frame_index+1;
    %     drawnow
    end
toc
end
%% Analyzing the angular frequency
% % movie_name
% v_omega(1:N,1:Obs_time_steps+delta_t/dt)=0;
% circle_center_x(1:3)=0;
% circle_center_y(1:3)=0;
% time_to_steady_state=600; %ms
% for i=1:N % Center of Circle of each particle's trajectory
%     if time(end)>time_to_steady_state
%         circle_center_x(i)=mean(x(i,time>time_to_steady_state));
%         circle_center_y(i)=mean(y(i,time>time_to_steady_state));
%     else
%         if i==1
%             warning('The simulation is not long enough for the system to reach stable state.')
%         end
% %         pause
%         circle_center_x(i)=mean(mean(x,1),2);
%         circle_center_y(i)=mean(mean(y,1),2);
%     end
% end
% 
% for k=1:Obs_time_steps+delta_t/dt % k=1:delta_t/dt does not move
%     cm_x=mean(x(:,k),1);
%     cm_y=mean(y(:,k),1);
%     
%     for i=1:N
%                         R_x=x(i,k)-cm_x;
%                         R_y=y(i,k)-cm_y;
% %         R_x=x(i,k)-circle_center_x(i);
% %         R_y=y(i,k)-circle_center_y(i);
%         R=[R_x R_y 0];
%         v=[v_x(i,k) v_y(i,k) 0];
%         cross_R_v=cross(R,v);
%         v_omega(i,k)=cross_R_v(3)/(norm(R))/v_0; %Normalized angular speed
%     end
% end
% figure;
% hold on
% % plot(time,v_omega(1,:))
% % plot(time,v_omega(2,:))
% % plot(time,v_omega(3,:))
% plot(time,movmean(v_omega(1,:),100000))
% plot(time,movmean(v_omega(2,:),100000))
% plot(time,movmean(v_omega(3,:),100000))
% xline(delta_t)
% title('Normalized Rotation Speed')
% xlabel('time (ms)')
% ylabel('v_\omega(rad/ms) /v_0 ')
% legend('1','2','3','\deltat')        
end

