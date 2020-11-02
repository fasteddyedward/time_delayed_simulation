%% 2nd stage subfunction of simulation.m

function [x,y,F_x,F_y,v_x,v_y,delta_x,delta_y,time]=second_stage_delayed_int(N,delta_t,dt,partition_time_steps,v_0,T,x_temp,y_temp,lth_partition,gamma,k_B,D,hard_collision,a)
x=x_temp;
y=y_temp;
%% Second stage: Delayed interaction starts. t=delta_t~Obs_time
delta_x(1:N)=0;
delta_y(1:N)=0;
for k=1:partition_time_steps
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
            %% We should test the duration of these for loops, since they might take up time
            %% Singulartiy
            %%% If e_x or e_y is 0, F_x & F_y will not be calculated as 0
            %%% numerically, so we have to put it by hand.
            F_x(i)=e_x(i)/sqrt(e_x(i)^2+e_y(i)^2)*v_0;
            F_y(i)=e_y(i)/sqrt(e_x(i)^2+e_y(i)^2)*v_0;
            if e_x(i)^2+e_y(i)^2==0 && e_x(i)==0
                F_x(i)=0;        end
            if e_x(i)^2+e_y(i)^2==0 && e_y(i)==0
                F_y(i)=0;        end
        %% Storing the particle's displacement (including the fluctuation)
        delta_x(i)=F_x(i)*dt+normrnd(0,sqrt(4*D*dt));
        delta_y(i)=F_y(i)*dt+normrnd(0,sqrt(4*D*dt));
        %% Updating particle position with Particle Interaction and Diffusion
        x(i,1+k+delta_t/dt)=x(i,k+delta_t/dt)+delta_x(i);
        y(i,1+k+delta_t/dt)=y(i,k+delta_t/dt)+delta_y(i);
        %% Including hard core interaction
        switch hard_collision
            case 'on'
                delta_x_attempt=x(i,k+1+delta_t/dt)-x(i,k+delta_t/dt);
                delta_y_attempt=y(i,k+1+delta_t/dt)-y(i,k+delta_t/dt);
                for j=1:N
                    if i>j % j has been updated to k+1+delta_t/dt, now updating i to k+1+delta_t/dt
                        diff_x=x(j,k+1+delta_t/dt)-x(i,k+1+delta_t/dt);
                        diff_y=y(j,k+1+delta_t/dt)-y(i,k+1+delta_t/dt);
                        if diff_x^2+diff_y^2 < (2*a)^2
                            %                             delta_x_attempt=x(i,k+1+delta_t/dt)-x(i,k+delta_t/dt);
                            %                             delta_y_attempt=y(i,k+1+delta_t/dt)-y(i,k+delta_t/dt);
                            x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-0.5*delta_x_attempt; % the ith particle at k+1+delta_t/dt (hitting j) goes backwards half its way
                            x(j,k+1+delta_t/dt)=x(j,k+1+delta_t/dt)+0.5*delta_x_attempt; % the jth particle at k+1+delta_t/dt (being hitted by i) goes forward half i's way
                            y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-0.5*delta_y_attempt;
                            y(j,k+1+delta_t/dt)=y(j,k+1+delta_t/dt)+0.5*delta_y_attempt;
                        end
                    elseif i<j % j has not been updated (now at k+delta_t/dt), now updating i to k+1+delta_t/dt
                        diff_x=x(j,k+delta_t/dt)-x(i,k+1+delta_t/dt);
                        diff_y=y(j,k+delta_t/dt)-y(i,k+1+delta_t/dt);
                        if diff_x^2+diff_y^2 < (2*a)^2
                            %                             delta_x_attempt=x(i,k+1+delta_t/dt)-x(i,k+delta_t/dt);
                            %                             delta_y_attempt=y(i,k+1+delta_t/dt)-y(i,k+delta_t/dt);
                            x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-0.5*delta_x_attempt; % the ith particle at k+1 (hitting j) goes backwards half its way
                            x(j,k+delta_t/dt)=x(j,k+delta_t/dt)    +0.5*delta_x_attempt; % the jth particle at k (being hitted by i) goes forward half i's way
                            y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-0.5*delta_y_attempt;
                            y(j,k+delta_t/dt)=y(j,k+delta_t/dt)    +0.5*delta_y_attempt;
                        end
                    end
                end
        end
        end
        %% What if the particles still overlap after the previous subsection?
        while 1==1
            %% Change I and J back to i, j
            check_relax(1:N,1:N)=0;
            b=0.51;
            for i=1:N
                for j=1:N
                    if j~=i % both i and j have been updated to k+1+delta_t/dt, now updating i to k+1+delta_t/dt
                        diff_x=x(j,k+1+delta_t/dt)-x(i,k+1+delta_t/dt);
                        diff_y=y(j,k+1+delta_t/dt)-y(i,k+1+delta_t/dt);
                        diff_r_sqr=diff_x^2+diff_y^2;
                        if diff_r_sqr < (2*a)^2
                            x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-b*diff_x/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr)); % the ith particle at k+1+delta_t/dt (hitting j) goes backwards half its way
                            x(j,k+1+delta_t/dt)=x(j,k+1+delta_t/dt)+b*diff_x/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr)); % the jth particle at k+1+delta_t/dt (being hitted by i) goes forward half i's way
                            y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-b*diff_y/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr));
                            y(j,k+1+delta_t/dt)=y(j,k+1+delta_t/dt)+b*diff_y/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr));
                        else 
                            check_relax(i,j)=1;
                        end
                    end
                end
            end

            if sum(sum(check_relax,2))==N^2-N % All particles have inter distance larger than 2a
                break
            else
%                 sum(sum(check_relax,2))
%                 k
            end
        end
end
x(:,1:size(x_temp,2)-1)=[];
y(:,1:size(y_temp,2)-1)=[];
%% Calculating the velocities at each timestep (from t=0+delta_t ~ Obs_time+delta_t)
v_x=diff(x,1,2)/dt;
v_y=diff(y,1,2)/dt;
time=(1+(lth_partition-2)*partition_time_steps+delta_t/dt:(lth_partition-1)*partition_time_steps+delta_t/dt)*dt;
end
    
    
    
    
    
