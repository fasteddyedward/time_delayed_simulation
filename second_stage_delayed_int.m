%% 2nd stage subfunction of simulation.m

function [x,y,F_x,F_y,v_x,v_y,delta_x,delta_y,time]=second_stage_delayed_int(N,delta_t,dt,partition_time_steps,v_0,T,x_temp,y_temp,lth_partition,gamma,k_B,D,hard_collision,a,b,int_delay,fixed_flag)
%% Introducing the intrinsic delay feature

frac=0.5;
d= 0.7754*a; % Distance from the center of the particle where the laser would shine on for maximum velocity. See laser_reduction.m and Martin_exp.m
%%
x=x_temp;
y=y_temp;
unit_vec = @(x,y) ([x,y]/norm([x,y]));
%% Second stage: Delayed interaction starts. t=delta_t~Obs_time
delta_x(1:N)=0;
delta_y(1:N)=0;
F_x(1:N)=0;
F_y(1:N)=0;
for k=1:partition_time_steps
    e_x(1:N)=0; % template value
    e_y(1:N)=0;
    diff_x(1:N,1:N)=0;
    diff_y(1:N,1:N)=0;
    %% Particle Interaction: calculating the F_x and F_y with position at delayed times
%         for i=N:-1:1
%             for j=N:-1:1
        for i=1:N
            for j=1:N
                if(j~=i)
                    diff_x(i,j)=x(j,k)-x(i,k);
                    diff_y(i,j)=y(j,k)-y(i,k);
                    if norm([diff_x(i,j),diff_y(i,j)])==0 && diff_x(i,j)==0
                        e_x(i)=e_x(i)+0;
                    else
                        e_x(i)=e_x(i)+diff_x(i,j)/norm([diff_x(i,j),diff_y(i,j)]);
                    end
                    if norm([diff_x(i,j),diff_y(i,j)])==0 && diff_y(i,j)==0
                        e_y(i)=e_y(i)+0 ;
                    else
                        e_y(i)=e_y(i)+diff_y(i,j)/norm([diff_x(i,j),diff_y(i,j)]);
                    end
                end
            end
                %% The Force and Singulartiy (We can actually get rid of Singularity part since F_x and F_y is nearly impossible to be 0)
                %%% If e_x or e_y is 0, F_x & F_y will not be calculated as 0
                %%% numerically, so we have to put it by hand.
                if int_delay==0 %% No intrinsic delay, then same as before.
                    e=unit_vec(e_x(i),e_y(i)); % This is the direction of the applied force
                    F_x(i)=e(1)*v_0;
                    F_y(i)=e(2)*v_0;
                    if a==0 %% The scenario below happens only when the particles are allowed to overlap each other
                        if norm(e_x(i),e_y(i))==0 && e_x(i)==0
                            F_x(i)=0;        end
                        if norm(e_x(i),e_y(i))==0 && e_y(i)==0
                            F_y(i)=0;        end
                    end
                    %% In case there is intrinsic delay
                else % Including finite intrinsic delay
                    e=unit_vec(e_x(i),e_y(i)); % This is the direction of the applied force
                    %                 r = d*e+[x(i,k+delta_t/dt), y(i,k+delta_t/dt)]-[x(i,k+(delta_t-int_delay)/dt), y(i,k+(delta_t-int_delay)/dt)]; % Error because of MATLAB numerical error
                    r = d*e+[x(i,k+delta_t/dt), y(i,k+delta_t/dt)]-[x(i,k+round((delta_t-int_delay)/dt)), y(i,k+round((delta_t-int_delay)/dt))]; % Need the round because of MATLAB numerical error
                    r_unit=unit_vec(r(1),r(2));
                    v=v_0*laser_reduction(norm(r)/a); % This is the reduction of the laser power based on Martin's model
                    F_x(i)=r_unit(1)*v;
                    F_y(i)=r_unit(2)*v;
                end
                
                %% Storing the particle's displacement (including the fluctuation)
                delta_x(i)=F_x(i)*dt+normrnd(0,sqrt(4*D*dt));
                delta_y(i)=F_y(i)*dt+normrnd(0,sqrt(4*D*dt));
                %% Updating particle position with Particle Interaction and Diffusion
                x(i,1+k+delta_t/dt)=x(i,k+delta_t/dt)+delta_x(i);
                y(i,1+k+delta_t/dt)=y(i,k+delta_t/dt)+delta_y(i);
                %% Including hard core interaction
                switch hard_collision
                    case 'method_1'
                        delta_x_attempt=x(i,k+1+delta_t/dt)-x(i,k+delta_t/dt);
                        delta_y_attempt=y(i,k+1+delta_t/dt)-y(i,k+delta_t/dt);
                        for j=1:N
                            if i>j % j has been updated to k+1+delta_t/dt, now updating i to k+1+delta_t/dt
                                Diff_x=x(j,k+1+delta_t/dt)-x(i,k+1+delta_t/dt);
                                Diff_y=y(j,k+1+delta_t/dt)-y(i,k+1+delta_t/dt);
                                if Diff_x^2+Diff_y^2 < (2*a)^2
                                    x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-0.5*delta_x_attempt; % the ith particle at k+1+delta_t/dt (hitting j) goes backwards half its way
                                    y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-0.5*delta_y_attempt;
                                    if fixed_flag(j)==0
                                    x(j,k+1+delta_t/dt)=x(j,k+1+delta_t/dt)+0.5*delta_x_attempt; % the jth particle at k+1+delta_t/dt (being hitted by i) goes forward half i's way
                                    y(j,k+1+delta_t/dt)=y(j,k+1+delta_t/dt)+0.5*delta_y_attempt;
                                    end
                                end
                            elseif i<j % j has not been updated (now at k+delta_t/dt), now updating i to k+1+delta_t/dt
                                Diff_x=x(j,k+delta_t/dt)-x(i,k+1+delta_t/dt);
                                Diff_y=y(j,k+delta_t/dt)-y(i,k+1+delta_t/dt);
                                if Diff_x^2+Diff_y^2 < (2*a)^2
                                    x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-0.5*delta_x_attempt; % the ith particle at k+1 (hitting j) goes backwards half its way
                                    y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-0.5*delta_y_attempt;
                                    if fixed_flag(j)==0
                                    x(j,k+delta_t/dt)=x(j,k+delta_t/dt)    +0.5*delta_x_attempt; % the jth particle at k (being hitted by i) goes forward half i's way
                                    y(j,k+delta_t/dt)=y(j,k+delta_t/dt)    +0.5*delta_y_attempt;
                                    end
                                end
                            end
                        end
                    case 'method_2'
                        delta_x_attempt=x(i,k+1+delta_t/dt)-x(i,k+delta_t/dt);
                        delta_y_attempt=y(i,k+1+delta_t/dt)-y(i,k+delta_t/dt);
                        for j=1:N
                            if i>j % j has been updated to k+1+delta_t/dt, now updating i to k+1+delta_t/dt
                                Diff_x=x(j,k+1+delta_t/dt)-x(i,k+1+delta_t/dt);
                                Diff_y=y(j,k+1+delta_t/dt)-y(i,k+1+delta_t/dt);
                                %                                 if Diff_x^2+Diff_y^2 < (2*a)^2                                
                                if norm([Diff_x Diff_y]) < 2*a
                                    reflect_normal_unit=unit_vec(Diff_x,Diff_y); % The normal vector of the reflective surface
                                    reflected_unit=unit_vec(delta_x_attempt,delta_y_attempt)-2*(unit_vec(delta_x_attempt,delta_y_attempt)*reflect_normal_unit')*reflect_normal_unit;
                                    
%                                     x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-0.5*delta_x_attempt; % the ith particle at k+1+delta_t/dt (hitting j) goes backwards half its way
%                                     y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-0.5*delta_y_attempt;
                                    x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-0.5*delta_x_attempt+0.5*norm([delta_x_attempt,delta_y_attempt])*reflected_unit(1); % the ith particle at k+1+delta_t/dt (hitting j) goes backwards half its way
                                    y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-0.5*delta_y_attempt+0.5*norm([delta_x_attempt,delta_y_attempt])*reflected_unit(2);
                                    if fixed_flag(j)==0
%                                     x(j,k+1+delta_t/dt)=x(j,k+1+delta_t/dt)+0.5*delta_x_attempt; % the jth particle at k+1+delta_t/dt (being hitted by i) goes forward half i's way
%                                     y(j,k+1+delta_t/dt)=y(j,k+1+delta_t/dt)+0.5*delta_y_attempt;
                                    x(j,k+1+delta_t/dt)=x(j,k+1+delta_t/dt)+0.5*norm([delta_x_attempt,delta_y_attempt])*reflect_normal_unit(1); % the jth particle at k+1+delta_t/dt (being hitted by i) goes forward half i's way
                                    y(j,k+1+delta_t/dt)=y(j,k+1+delta_t/dt)+0.5*norm([delta_x_attempt,delta_y_attempt])*reflect_normal_unit(2);
                                    end
                                end
                            elseif i<j % j has not been updated (now at k+delta_t/dt), now updating i to k+1+delta_t/dt
                                Diff_x=x(j,k+delta_t/dt)-x(i,k+1+delta_t/dt);
                                Diff_y=y(j,k+delta_t/dt)-y(i,k+1+delta_t/dt);
                                %                                 if Diff_x^2+Diff_y^2 < (2*a)^2
                                if norm([Diff_x Diff_y]) < 2*a
                                    reflect_normal_unit=unit_vec(Diff_x,Diff_y); % The normal vector of the reflective surface
                                    reflected_unit=unit_vec(delta_x_attempt,delta_y_attempt)-2*(unit_vec(delta_x_attempt,delta_y_attempt)*reflect_normal_unit')*reflect_normal_unit;
                                    
%                                     x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-0.5*delta_x_attempt; % the ith particle at k+1 (hitting j) goes backwards half its way
%                                     y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-0.5*delta_y_attempt;
                                    x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-0.5*delta_x_attempt+0.5*norm([delta_x_attempt,delta_y_attempt])*reflected_unit(1); % the ith particle at k+1 (hitting j) goes backwards half its way
                                    y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-0.5*delta_y_attempt+0.5*norm([delta_x_attempt,delta_y_attempt])*reflected_unit(2);
                                    if fixed_flag(j)==0
%                                     x(j,k+delta_t/dt)=x(j,k+delta_t/dt)    +0.5*delta_x_attempt; % the jth particle at k (being hitted by i) goes forward half i's way
%                                     y(j,k+delta_t/dt)=y(j,k+delta_t/dt)    +0.5*delta_y_attempt;
                                    x(j,k+delta_t/dt)=x(j,k+delta_t/dt)    +0.5*norm([delta_x_attempt,delta_y_attempt])*reflect_normal_unit(1); % the jth particle at k (being hitted by i) goes forward half i's way
                                    y(j,k+delta_t/dt)=y(j,k+delta_t/dt)    +0.5*norm([delta_x_attempt,delta_y_attempt])*reflect_normal_unit(2);
                                    end
                                end
                            end
                        end
                end
                if fixed_flag(i)==1 %% particle is fixed, overwrite the updated x(i,1+k+delta_t/dt) with x(i,k+delta_t/dt)
                    x(i,1+k+delta_t/dt)=x(i,k+delta_t/dt);
                    y(i,1+k+delta_t/dt)=y(i,k+delta_t/dt);
                end
        end
        
        %% What if the particles still overlap after the previous subsection?
        switch hard_collision
            case {'method_1', 'method_2'}
                while 1==1
                    
                    check_relax(1:N,1:N)=0;
                    for i=1:N
                        for j=1:N
                            if j~=i % both i and j have been updated to k+1+delta_t/dt, now updating i to k+1+delta_t/dt
                                diff_x=x(j,k+1+delta_t/dt)-x(i,k+1+delta_t/dt);
                                diff_y=y(j,k+1+delta_t/dt)-y(i,k+1+delta_t/dt);
                                diff_r_sqr=diff_x^2+diff_y^2;
                                if diff_r_sqr < (2*a)^2
                                    %% Note that here we comment out the fixed_flag(i)==0 condition or else the calculations will take forever.
                                    %                             if fixed_flag(i)==0
                                    x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-b*diff_x/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr)); % the ith particle at k+1+delta_t/dt (hitting j) goes backwards half its way
                                    y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-b*diff_y/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr));
                                    %                             elseif fixed_flag(i)==1
                                    %                                 x(i,1+k+delta_t/dt)=x(i,k+delta_t/dt);
                                    %                                 y(i,1+k+delta_t/dt)=y(i,k+delta_t/dt);
                                    %                             end
                                    %                             if fixed_flag(j)==0
                                    x(j,k+1+delta_t/dt)=x(j,k+1+delta_t/dt)+b*diff_x/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr)); % the jth particle at k+1+delta_t/dt (being hitted by i) goes forward half i's way
                                    y(j,k+1+delta_t/dt)=y(j,k+1+delta_t/dt)+b*diff_y/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr));
                                    %                             elseif fixed_flag(j)==1
                                    %                                 x(j,1+k+delta_t/dt)=x(j,k+delta_t/dt);
                                    %                                 y(j,1+k+delta_t/dt)=y(j,k+delta_t/dt);
                                    %                             end
                                else
                                    check_relax(i,j)=1;
                                end
                            end
                        end
                    end
                    
                    if sum(sum(check_relax,2))==(N^2-N) % All particles have inter distance larger than 2a
                        for i=1:N
                            if fixed_flag(i)==1 %% particle is fixed, overwrite the updated x(i,1+k+delta_t/dt) with x(i,k+delta_t/dt)
                                x(i,1+k+delta_t/dt)=x(i,k+delta_t/dt);
                                y(i,1+k+delta_t/dt)=y(i,k+delta_t/dt);
                            end
                        end
                        break
                    else
                        %                 diff_r_sqr-(2*a)^2
                        %                 sum(sum(check_relax,2))
%                         k
                    end
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
    
    
    
    
    
