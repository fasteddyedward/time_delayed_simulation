%% 1st stage subfunction of simulation.m

function [x,y,v_x,v_y,time]=first_stage_pure_diffusion(N,delta_t,dt,T,x_init,y_init,gamma,k_B,D,hard_collision,a)

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
            
            %% Including hard core interaction
%             if hard_collision=='on'
switch hard_collision 
    case 'on'
                delta_x_attempt=x(i,k+1)-x(i,k);
                delta_y_attempt=y(i,k+1)-y(i,k);
                for j=1:N
                    if i>j % j has been updated to k+1, now updating i to k+1
                        diff_x=x(j,k+1)-x(i,k+1);
                        diff_y=y(j,k+1)-y(i,k+1);
                        if diff_x^2+diff_y^2 < (2*a)^2
                            %                             delta_x_attempt=x(i,k+1)-x(i,k);
                            %                             delta_y_attempt=y(i,k+1)-y(i,k);
                            x(i,k+1)=x(i,k+1)-0.5*delta_x_attempt; % the ith particle at k+1 (hitting j) goes backwards half its way
                            x(j,k+1)=x(j,k+1)+0.5*delta_x_attempt; % the jth particle at k+1 (being hitted by i) goes forward half i's way
                            y(i,k+1)=y(i,k+1)-0.5*delta_y_attempt;
                            y(j,k+1)=y(j,k+1)+0.5*delta_y_attempt;
                        end
                    elseif i<j % j has not been updated (now at k), now updating i to k+1
                        diff_x=x(j,k)-x(i,k+1);
                        diff_y=y(j,k)-y(i,k+1);
                        if diff_x^2+diff_y^2 < (2*a)^2
                            %                             delta_x_attempt=x(i,k+1)-x(i,k);
                            %                             delta_y_attempt=y(i,k+1)-y(i,k);
                            x(i,k+1)=x(i,k+1)-0.5*delta_x_attempt; % the ith particle at k+1 (hitting j) goes backwards half its way
                            x(j,k)=x(j,k)    +0.5*delta_x_attempt; % the jth particle at k (being hitted by i) goes forward half i's way
                            y(i,k+1)=y(i,k+1)-0.5*delta_y_attempt;
                            y(j,k)=y(j,k)    +0.5*delta_y_attempt;
                        end
                    end
                end
            end
    end
end

%% Calculating the velocities at each timestep (from t=0+delta_t ~ Obs_time+delta_t)
v_x=diff(x,1,2)/dt;
v_y=diff(y,1,2)/dt;
time=(1:delta_t/dt)*dt;
end
    
    
    
    
    
