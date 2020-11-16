function [Movie_Vector]=make_movies_plots(N,delta_t,dt,partition_time_steps,x,y,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,a,force_tracks)
switch movie_create
    case 'off'
    Movie_Vector=[];
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
    figure(100)
    for i=1:N
        for j=i+1:N
            std=std+(mean_x(i)-mean_x(j))^2+(mean_y(i)-mean_y(j))^2;
        end
    end
    std=sqrt(std/(N-1));
    
    %     for k=1+delta_t/dt:animation_interval:Obs_time_steps+delta_t/dt % frames with k=1:delta_t/dt do not move
    %         MovieVector(1:((Obs_time_steps+delta_t/dt)/animation_interval))=0;
%     for k=1:animation_interval:Obs_time_steps+delta_t/dt % frames with k=1:delta_t/dt only diffuses
    for k=1:animation_interval:partition_time_steps % frames with k=1:delta_t/dt only diffuses
        switch leave_trace
            case 'on'
            case 'off'
                clf
        end
        grid on
        %% Calculating center of mass and the std deviation of the particles
        %         cm_x=mean(x(:,k),1);
        %         cm_y=mean(y(:,k),1);
        F_x(1:N)=0;
        F_y(1:N)=0;
        for i=1:N
            hold on
            h(i)=scatter(x(i,k),y(i,k),'filled');
            switch ghost
                case 'on'
                    if k>delta_t/dt
                        scatter(x(i,k-delta_t/dt),y(i,k-delta_t/dt),'o')
                    end
                case 'off'
            end
            switch force_tracks 
                case 'on'
%                 F_x(i)=x(i,k)-x(i,k-1);
%                 F_y(i)=y(i,k)-y(i,k-1);
                F_x(i)=x(i,k+1)-x(i,k);
                F_y(i)=y(i,k+1)-y(i,k);
                hold on
                magnify=5*10^3;
                quiver(x(i,k),y(i,k),F_x(i)*magnify,F_y(i)*magnify,'k');
            end
        end

        
        
        %         plot(cm_x,cm_y,'.')
%         plot(mean(x(:,k),1),mean(y(:,k),1),'.')
        %% Don't use legend since it costs too much time
%         legend('1','2','3','CM','Location','northeastoutside')
        title(['Movie with ',axis_choice,' frame'])

        %% Axis of choice
        switch axis_choice
            case 'lab'
%                 axis_scale=[min(min(x)) max(max(x)) min(min(y)) max(max(y))];
                axis(axis_scale);
            case 'cm'
                %                                 axis([cm_x-magnify*std   cm_x+magnify*std   cm_y-magnify*std   cm_y+magnify*std]);
                %                 axis([cm_x-magnify*std   cm_x+magnify*std   cm_y-magnify*std   cm_y+magnify*std]);
                axis([mean(x(:,k),1)-magnify*std   mean(x(:,k),1)+magnify*std   mean(y(:,k),1)-magnify*std   mean(y(:,k),1)+magnify*std]);
                axis([mean(x(:,k),1)-magnify*std   mean(x(:,k),1)+magnify*std   mean(y(:,k),1)-magnify*std   mean(y(:,k),1)+magnify*std]);
                
        end
        axis square
        ax=get(gca);
        marker_size=200/0.4*a/diff(ax.XLim);
        for i=1:N
        h(i).SizeData=marker_size^2;
        end
        %%% For making movie
        Movie_Vector(movie_frame_index)=getframe(gcf);
        movie_frame_index=movie_frame_index+1;
    %     drawnow
    end
toc
end
end

