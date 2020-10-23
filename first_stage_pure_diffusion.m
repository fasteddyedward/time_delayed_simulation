%% This file is the modulated version of time_delay_proto. 
%% Basically copy everything below the section %% Checking comparison of Driving force and diffusion to modulized_time_delay_proto.m
%% And copy everything above to main.m. 
%% Just remember to comment the variables that were inputted to this file.

function [x,y,v_x,v_y,time]=first_stage_pure_diffusion(N,delta_t,dt,T,x_init,y_init,gamma,k_B,D)
%% Setup for Running the program
% Obs_time=Obs_time_steps*dt;

%% Coefficients and parameters
% v_0= 10^-5; % mm/ms
% gamma=1; % Stoke's drag, gamma=6*pi*eta*a
% T=300; % Kelvin
% k_B=1; % Boltzmann constant
% D=k_B*T/gamma; % Diffusitivity
%% Checking comparison of Driving force and diffusion
% delta_x=F_x*dt+normrnd(0,sqrt(4*D*dt)) ~ v_0*dt+sqrt(4*D*dt)
% ['The ratio v_0*dt/sqrt(4*D*dt) is ',num2str(v_0*dt/(sqrt(4*D*dt)))]
%% Warning for choosing time steps 
%%% If the Obs_time is smaller than time delay delta_t, the F_x and F_y we
%%% see will be flat, since the particle hasn't started to move yet.
% if Obs_time<delta_t
%     warning('The observation time is less than time delay, please choose a larger Obs_time')
% end


%% Start solving equation of motion
for i=1:N
    x(i,1)=x_init(i);
    y(i,1)=y_init(i);
end
%% First stage: Diffusion. t=0 ~ delta_t
F_x(1:N)=0;
F_y(1:N)=0;
% for k=1:1+delta_t/dt
for k=1:delta_t/dt    
    for i=1:N
            x(i,k+1)=x(i,k)+normrnd(0,sqrt(4*D*dt));
            y(i,k+1)=y(i,k)+normrnd(0,sqrt(4*D*dt));
    end
end

%% Second stage: Delayed interaction starts. t=delta_t~Obs_time
% delta_x(1:N)=0;
% delta_y(1:N)=0;
% for k=1:Obs_time_steps
%     e_x(1:N)=0; % template value
%     e_y(1:N)=0;
%     diff_x(1:N,1:N)=0;
%     diff_y(1:N,1:N)=0;
%     %% Particle Interaction: calculating the F_x and F_y with position at delayed times
% %     if v_0==0 %% Here I simplify the calculation for pure diffusion v_0=0 to save time
% %         F_x(i,k+delta_t/dt)=0;
% %         F_y(i,k+delta_t/dt)=0;
% %     else
%         for i=1:N
%             for j=1:N
%                 if(j~=i)
%                     diff_x(i,j)=x(j,k)-x(i,k);
%                     diff_y(i,j)=y(j,k)-y(i,k);
%                     if diff_x(i,j)^2+diff_y(i,j)^2==0 && diff_x(i,j)==0
%                         e_x(i)=e_x(i)+0;
%                     else
%                         e_x(i)=e_x(i)+diff_x(i,j)/sqrt(diff_x(i,j)^2+diff_y(i,j)^2);
%                     end
%                     if diff_x(i,j)^2+diff_y(i,j)^2==0 && diff_y(i,j)==0
%                         e_y(i)=e_y(i)+0 ;
%                     else
%                         e_y(i)=e_y(i)+diff_y(i,j)/sqrt(diff_x(i,j)^2+diff_y(i,j)^2);
%                     end
%                 end
%             end
%             %% We should test the duration of these for loops, since they might take up time
%             %% Singulartiy
%             %%% If e_x or e_y is 0, F_x & F_y will not be calculated as 0
%             %%% numerically, so we have to put it by hand.
%             F_x(i)=e_x(i)/sqrt(e_x(i)^2+e_y(i)^2)*v_0;
%             F_y(i)=e_y(i)/sqrt(e_x(i)^2+e_y(i)^2)*v_0;
%             if e_x(i)^2+e_y(i)^2==0 && e_x(i)==0
%                 F_x(i)=0;        end
%             if e_x(i)^2+e_y(i)^2==0 && e_y(i)==0
%                 F_y(i)=0;        end
%         %% Storing the particle's displacement (including the fluctuation)
%         delta_x(i)=F_x(i)*dt+normrnd(0,sqrt(4*D*dt));
%         delta_y(i)=F_y(i)*dt+normrnd(0,sqrt(4*D*dt));
%         %% Updating particle position with Particle Interaction and Diffusion
%         x(i,1+k+delta_t/dt)=x(i,k+delta_t/dt)+delta_x(i);
%         y(i,1+k+delta_t/dt)=y(i,k+delta_t/dt)+delta_y(i);
%         end
% end

%% Calculating the velocities at each timestep (from t=0+delta_t ~ Obs_time+delta_t)
v_x=diff(x,1,2)/dt;
v_y=diff(y,1,2)/dt;
time=(1:delta_t/dt)*dt;
end
    
    
    
    
    
