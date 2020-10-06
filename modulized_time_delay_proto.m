%% This file is the modulated version of time_delay_proto. 
%% Basically copy everything below the section %% Checking comparison of Driving force and diffusion to modulized_time_delay_proto.m
%% And copy everything above to main.m. 
%% Just remember to comment the variables that were inputted to this file.

function [x,y,F_x,F_y,v_x,v_y,delta_x,delta_y,time]=modulized_time_delay_proto(N,delta_t,dt,Obs_time_steps,v_0,T,x_init,y_init)
% clear;
% close all
%% Setup for Running the program
% N=3; % number of particles in the play delta_t=5; % ms dt=10^-4; % ms
% Obs_time_steps=10^7
Obs_time=Obs_time_steps*dt;

%% Coefficients and parameters
% v_0= 10^-5; % mm/ms
gamma=6*pi*1.0016*10^-3*10^-6; % Stoke's drag, gamma=6*pi*eta*a
% T=300; % Kelvin
k_B=10^-23; % Boltzmann constant
D=k_B*T/gamma; % Diffusitivity
%% Parameters for making the movies
% magnify=1 control_animation_interval=10^4 movie_create='on' ghost='off'
% axis_choice='lab'; %'cm' or 'lab'
%% Checking comparison of Driving force and diffusion
% delta_x=F_x*dt+normrnd(0,sqrt(2*D*dt)) ~ v_0*dt+sqrt(2*D*dt)
['The ratio v_0*dt/sqrt(4*D*dt) is ',num2str(v_0*dt/(sqrt(4*D*dt)))]
%% Warning for choosing time steps 
%%% If the Obs_time is smaller than time delay delta_t, the F_x and F_y we
%%% see will be flat, since the particle hasn't started to move yet.
if Obs_time<delta_t
    warning('The observation time is less than time delay, please choose a larger Obs_time')
    pause
end


%% Start solving equation of motion
for i=1:N
    x(i,1)=x_init(i);
    y(i,1)=y_init(i);
end
%% First stage: Diffusion. t=0 ~ delta_t
F_x(1:N,1:Obs_time_steps+delta_t/dt)=0;
F_y(1:N,1:Obs_time_steps+delta_t/dt)=0;
delta_x(1:N,1:Obs_time_steps+delta_t/dt)=0;
delta_y(1:N,1:Obs_time_steps+delta_t/dt)=0;
for k=1:1+delta_t/dt
    for i=1:N
%             delta_x(i,k)=normrnd(0,sqrt(2*D*dt));
%             delta_y(i,k)=normrnd(0,sqrt(2*D*dt));
            delta_x(i,k)=normrnd(0,sqrt(4*D*dt));
            delta_y(i,k)=normrnd(0,sqrt(4*D*dt));
            x(i,k+1)=x(i,k)+delta_x(i,k);
            y(i,k+1)=y(i,k)+delta_y(i,k);
    end
end

%% Second stage: Delayed interaction starts. t=delta_t~Obs_time
for k=1:Obs_time_steps
    e_x(1:N)=0; % template value
    e_y(1:N)=0;
    diff_x(1:N,1:N)=0;
    diff_y(1:N,1:N)=0;
    %% Particle Interaction: calculating the F_x and F_y with position at delayed times
%     if v_0==0 %% Here I simplify the calculation for pure diffusion v_0=0 to save time
%         F_x(i,k+delta_t/dt)=0;
%         F_y(i,k+delta_t/dt)=0;
%     else
        for i=1:N
            for j=1:N
                if(j~=i)
                    diff_x(i,j)=x(j,k)-x(i,k);
                    diff_y(i,j)=y(j,k)-y(i,k);
                    if diff_x(i,j)^2+diff_y(i,j)^2==0 && diff_x(i,j)==0
                        e_x(i)=e_x(i)+0;
                    else
                        e_x(i)=e_x(i)+diff_x(i,j)/sqrt(diff_x(i,j)^2+diff_y(i,j)^2);
                    end
                    if diff_x(i,j)^2+diff_y(i,j)^2==0 && diff_y(i,j)==0
                        e_y(i)=e_y(i)+0 ;
                    else
                        e_y(i)=e_y(i)+diff_y(i,j)/sqrt(diff_x(i,j)^2+diff_y(i,j)^2);
                    end
                end
            end
            %% Singulartiy
            %%% If e_x or e_y is 0, F_x & F_y will not be calculated as 0
            %%% numerically, so we have to put it by hand.
            F_x(i,k+delta_t/dt)=e_x(i)/sqrt(e_x(i)^2+e_y(i)^2)*v_0;
            F_y(i,k+delta_t/dt)=e_y(i)/sqrt(e_x(i)^2+e_y(i)^2)*v_0;
            if e_x(i)^2+e_y(i)^2==0 && e_x(i)==0
                F_x(i,k+delta_t/dt)=0;        end
            if e_x(i)^2+e_y(i)^2==0 && e_y(i)==0
                F_y(i,k+delta_t/dt)=0;        end
%     end
        %% Storing the particle's displacement (including the fluctuation)
%         delta_x(i,k+delta_t/dt)=F_x(i,k+delta_t/dt)*dt+normrnd(0,sqrt(2*D*dt));
%         delta_y(i,k+delta_t/dt)=F_y(i,k+delta_t/dt)*dt+normrnd(0,sqrt(2*D*dt));
        delta_x(i,k+delta_t/dt)=F_x(i,k+delta_t/dt)*dt+normrnd(0,sqrt(4*D*dt));
        delta_y(i,k+delta_t/dt)=F_y(i,k+delta_t/dt)*dt+normrnd(0,sqrt(4*D*dt));
        %% Updating particle position with Particle Interaction and Diffusion
        x(i,1+k+delta_t/dt)=x(i,k+delta_t/dt)+delta_x(i,k+delta_t/dt);
        y(i,1+k+delta_t/dt)=y(i,k+delta_t/dt)+delta_y(i,k+delta_t/dt);
        end
    end
% end

%% Calculating the velocities at each timestep (from t=0+delta_t ~ Obs_time+delta_t)
v_x=delta_x/dt;
v_y=delta_y/dt;
time=(1:Obs_time_steps+delta_t/dt)*dt;



% 
% 
% %% Plotting figures
%     %%% Plotting v_x
% figure
% hold on
% for i=1:N
%     plot(time,v_x(i,:))
% end
% title('v_x')
% legend('1','2','3')
% 
%     %%% Plotting v_y
% figure
% hold on
% for i=1:N
%     plot(time,v_y(i,:))
% end
% title('v_y')
% legend('1','2','3')
% %%
% for i=1:N
%     figure
%     hold on
%     plot(x(i,:))
%     plot(y(i,:))
%     title(['Position of particle ',num2str(i)])
%     legend('x','y')
% end
% %% Making movie of the partilces
% tic
% switch movie_create
%     case 'off'
%         'No movie created.'
%     case 'on'      
%     animation_interval=control_animation_interval;
%     % figure('visible','off'); % turning off the figure actually increases the
%     % time
% 
%     movie_frame_index=1; % For movie recording, increases 1 for each recorded frame
%     mean_x=mean(x(:,1+delta_t/dt:end),2); % We take the part of x and y after t=delta_t 
%     mean_y=mean(y(:,1+delta_t/dt:end),2);
%     std=0; % This is for setting the axis
%     figure
%     for i=1:N
%         for j=i+1:N
%             std=std+(mean_x(i)-mean_x(j))^2+(mean_y(i)-mean_y(j))^2;
%         end
%     end
%     std=sqrt(std/(N-1));
% 
% %     for k=1+delta_t/dt:animation_interval:Obs_time_steps+delta_t/dt % frames with k=1:delta_t/dt do not move
%     for k=1:animation_interval:Obs_time_steps+delta_t/dt % frames with k=1:delta_t/dt only diffuses
% %         clf
%         grid on
%         %% Calculating center of mass and the std deviation of the particles
%         cm_x=mean(x(:,k),1);
%         cm_y=mean(y(:,k),1);
%         %% Running std for axis, but I guess it's not a good idea
%     %     std=0; for i=1:N
%     %         for j=i+1:N
%     %             std=std+(x(i,k)-x(j,k))^2+(y(i,k)-y(j,k))^2;
%     % %         std=std+(x(i,k)-cm_x)^2+(y(i,k)-cm_y)^2;
%     %         end
%     %     end std=sqrt(std/(N-1));
%         for i=1:N
%             hold on
%             scatter(x(i,k),y(i,k),'filled')
%             switch ghost
%                 case 'on'
%                     if k>delta_t/dt
%                         scatter(x(i,k-delta_t/dt),y(i,k-delta_t/dt),'o')
%                     end
%                 case 'off'
%             end
%         end
%         plot(cm_x,cm_y,'.')
%         %% Don't use legend since it costs too much time
% %         legend('1','2','3','CM','Location','northeastoutside')
%         title(['Movie with ',axis_choice,' frame'])
% 
%         %% Axis of choice
%         switch axis_choice
%             case 'lab'
%                 axis_scale=[min(min(x)) max(max(x)) min(min(y)) max(max(y))];
%                 axis(axis_scale);
%             case 'cm'
%                 axis([cm_x-magnify*std   cm_x+magnify*std   cm_y-magnify*std   cm_y+magnify*std]);
%         end
% 
%         %%% For making movie
%         MovieVector(movie_frame_index)=getframe(gcf);    
%         movie_frame_index=movie_frame_index+1;
%     %     drawnow
%     end
% toc
% 
%     %% Save video 
%     myWriter=VideoWriter('collection');
%     % myWriter=VideoWriter(['delta_t=',num2str(delta_t),', ',axis_choice,' frame, Obs_time_steps=',num2str(Obs_time_steps),', log(dt)=',num2str(log10(dt))] );
%     myWriter.FrameRate=20;
%     open(myWriter);
%     writeVideo(myWriter,MovieVector);
%     close(myWriter);
% end
% %% Analyzing the angular frequency
% v_omega(1:N,1:Obs_time_steps+delta_t/dt)=0;
% for k=1:Obs_time_steps+delta_t/dt % k=1:delta_t/dt does not move
%     cm_x=mean(x(:,k),1);
%     cm_y=mean(y(:,k),1);
%     for i=1:N
%         R_x=x(i,k)-cm_x;
%         R_y=y(i,k)-cm_y;
%         R=[R_x R_y 0];
%         v=[v_x(i,k) v_y(i,k) 0];
%         cross_R_v=cross(R,v);
%         v_omega(i,k)=cross_R_v(3)/(norm(R))/v_0; %Normalized angular speed
%     end
% end
% figure;
% hold on
% plot(time,v_omega(1,:))
% plot(time,v_omega(2,:))
% plot(time,v_omega(3,:))
% xline(delta_t)
% title('Normalized Rotation Speed')
% xlabel('time (ms)')
% ylabel('v_\omega(rad/ms) /v_0 ')
% legend('1','2','3','\deltat')        
    
end
    
    
    
    
    
