%% 1st stage subfunction of simulation.m

function [x,y,v_x,v_y,time]=first_stage_pure_diffusion(N,delta_t,dt,T,x_init,y_init,gamma,k_B,D,hard_collision,a,b,fixed_flag)

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
        if fixed_flag(i)==0
            x(i,k+1)=x(i,k)+normrnd(0,sqrt(4*D*dt));
            y(i,k+1)=y(i,k)+normrnd(0,sqrt(4*D*dt));
            
            %% Including hard core interaction
            %             if hard_collision=='on'
            switch hard_collision
                case 'on'
                    delta_x_attempt=x(i,k+1)-x(i,k);
                    delta_y_attempt=y(i,k+1)-y(i,k);
                    for j=1:N
%                         if fixed_flag(j)==0 % j is mobile
                            if i>j % j has been updated to k+1, now updating i to k+1
                                Diff_x=x(j,k+1)-x(i,k+1);
                                Diff_y=y(j,k+1)-y(i,k+1);
                                if Diff_x^2+Diff_y^2 < (2*a)^2
                                    x(i,k+1)=x(i,k+1)-0.5*delta_x_attempt; % the ith particle at k+1 (hitting j) goes backwards half its way
                                    y(i,k+1)=y(i,k+1)-0.5*delta_y_attempt;
                                    if fixed_flag(j)==0 %  j is mobile; if j==1 then j doesn't move, the line under will not execute
                                    x(j,k+1)=x(j,k+1)+0.5*delta_x_attempt; % the jth particle at k+1 (being hitted by i) goes forward half i's way
                                    y(j,k+1)=y(j,k+1)+0.5*delta_y_attempt;
                                    end
                                end
                            elseif i<j % j has not been updated (now at k), now updating i to k+1
                                Diff_x=x(j,k)-x(i,k+1);
                                Diff_y=y(j,k)-y(i,k+1);
                                if Diff_x^2+Diff_y^2 < (2*a)^2
                                    x(i,k+1)=x(i,k+1)-0.5*delta_x_attempt; % the ith particle at k+1 (hitting j) goes backwards half its way
                                    y(i,k+1)=y(i,k+1)-0.5*delta_y_attempt;
                                    if fixed_flag(j)==0 %  j is mobile; if j==1 then j doesn't move, the line under will not execute
                                    x(j,k)=x(j,k)    +0.5*delta_x_attempt; % the jth particle at k (being hitted by i) goes forward half i's way
                                    y(j,k)=y(j,k)    +0.5*delta_y_attempt;
                                    end
                                end
                            end
%                         elseif j==1 %% j is fixed
%                             if i>j % j has been updated to k+1, now updating i to k+1
%                                 Diff_x=x(j,k+1)-x(i,k+1);
%                                 Diff_y=y(j,k+1)-y(i,k+1);
%                                 if Diff_x^2+Diff_y^2 < (2*a)^2
%                                     x(i,k+1)=x(i,k+1)-0.5*delta_x_attempt; % the ith particle at k+1 (hitting j) goes backwards half its way
% %                                     x(j,k+1)=x(j,k+1)+0.5*delta_x_attempt; % the jth particle at k+1 (being hitted by i) goes forward half i's way
%                                     y(i,k+1)=y(i,k+1)-0.5*delta_y_attempt;
% %                                     y(j,k+1)=y(j,k+1)+0.5*delta_y_attempt;
%                                 end
%                             elseif i<j % j has not been updated (now at k), now updating i to k+1
%                                 Diff_x=x(j,k)-x(i,k+1);
%                                 Diff_y=y(j,k)-y(i,k+1);
%                                 if Diff_x^2+Diff_y^2 < (2*a)^2
%                                     x(i,k+1)=x(i,k+1)-0.5*delta_x_attempt; % the ith particle at k+1 (hitting j) goes backwards half its way
% %                                     x(j,k)=x(j,k)    +0.5*delta_x_attempt; % the jth particle at k (being hitted by i) goes forward half i's way
%                                     y(i,k+1)=y(i,k+1)-0.5*delta_y_attempt;
% %                                     y(j,k)=y(j,k)    +0.5*delta_y_attempt;
%                                 end
%                             end
%                         end
                    end
            end
        end
    end
    %% What if the particles still overlap after the previous subsection?
    while 1==1
        check_relax(1:N,1:N)=0;
        %             b=0.6;
        for i=1:N
            for j=1:N
                if fixed_flag(i)==0 % i is mobile
                    if j~=i % both i and j have been updated to k+1+delta_t/dt, now updating i to k+1+delta_t/dt
                        Diff_x=x(j,k+1)-x(i,k+1);
                        Diff_y=y(j,k+1)-y(i,k+1);
                        diff_r_sqr=Diff_x^2+Diff_y^2;
                        if diff_r_sqr < (2*a)^2
                            x(i,k+1)=x(i,k+1)-b*Diff_x/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr)); % the ith particle at k+1+delta_t/dt (hitting j) goes backwards half its way
                            y(i,k+1)=y(i,k+1)-b*Diff_y/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr));
                            if fixed_flag(j)==0 % j is mobile; if j==1 then j doesn't move, the line under will not execute
                            x(j,k+1)=x(j,k+1)+b*Diff_x/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr)); % the jth particle at k+1+delta_t/dt (being hitted by i) goes forward half i's way
                            y(j,k+1)=y(j,k+1)+b*Diff_y/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr));
                            end
                        else
                            check_relax(i,j)=1;
                        end
                    end
                end
            end
        end
        
        if sum(sum(check_relax,2))==N^2-N % All particles have inter distance larger than 2a
            break
        else
            sum(sum(check_relax,2))
            k
        end
    end
end

%% Calculating the velocities at each timestep (from t=0+delta_t ~ Obs_time+delta_t)
v_x=diff(x,1,2)/dt;
v_y=diff(y,1,2)/dt;
time=(1:delta_t/dt)*dt;
end
    
    
    
    
    
