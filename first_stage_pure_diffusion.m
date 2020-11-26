%% 1st stage subfunction of simulation.m

function [x,y,v_x,v_y,time]=first_stage_pure_diffusion(N,delta_t,dt,T,x_init,y_init,gamma,k_B,D,hard_collision,a,b,fixed_flag)
epsilon=0.0001*a;
%% Start solving equation of motion
x(1:N,1:delta_t/dt+1)=0;
y(1:N,1:delta_t/dt+1)=0;
for i=1:N
    x(i,1)=x_init(i);
    y(i,1)=y_init(i);
end

unit_vec = @(x,y) ([x,y]/norm([x,y]));
%% First stage: Diffusion. t=0 ~ delta_t
% F_x(1:N)=0;
% F_y(1:N)=0;
% for k=1:1+delta_t/dt
for k=1:delta_t/dt    
    %% Diffusion Process
    for i=1:N
%         if fixed_flag(i)==1 %% particle is fixed
%             x(i,k+1)=x(i,k);
%             y(i,k+1)=y(i,k);
%         elseif fixed_flag(i)==0 %% particle is mobile
            x(i,k+1)=x(i,k)+normrnd(0,sqrt(4*D*dt));
            y(i,k+1)=y(i,k)+normrnd(0,sqrt(4*D*dt));           
            %% Including hard core interaction
            %             if hard_collision=='on'
            switch hard_collision
                case 'method_1'
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
                            %
                    end
                case 'method_2'
                    delta_x_attempt=x(i,k+1)-x(i,k);
                    delta_y_attempt=y(i,k+1)-y(i,k);
                    for j=1:N
                        if i>j % j has been updated to k+1+delta_t/dt, now updating i to k+1+delta_t/dt
                            Diff_x=x(j,k+1)-x(i,k+1);
                            Diff_y=y(j,k+1)-y(i,k+1);
                            %                                 if Diff_x^2+Diff_y^2 < (2*a)^2
                            if norm([Diff_x Diff_y]) < 2*a
                                reflect_normal_unit=unit_vec(Diff_x,Diff_y); % The normal vector of the reflective surface
                                reflected_unit=unit_vec(delta_x_attempt,delta_y_attempt)-2*(unit_vec(delta_x_attempt,delta_y_attempt)*reflect_normal_unit')*reflect_normal_unit;
                                
                                %                                     x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-0.5*delta_x_attempt; % the ith particle at k+1+delta_t/dt (hitting j) goes backwards half its way
                                %                                     y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-0.5*delta_y_attempt;
                                x(i,k+1)=x(i,k+1)-0.5*delta_x_attempt+0.5*norm([delta_x_attempt,delta_y_attempt])*reflected_unit(1); % the ith particle at k+1+delta_t/dt (hitting j) goes backwards half its way
                                y(i,k+1)=y(i,k+1)-0.5*delta_y_attempt+0.5*norm([delta_x_attempt,delta_y_attempt])*reflected_unit(2);
                                if fixed_flag(j)==0
                                    %                                     x(j,k+1+delta_t/dt)=x(j,k+1+delta_t/dt)+0.5*delta_x_attempt; % the jth particle at k+1+delta_t/dt (being hitted by i) goes forward half i's way
                                    %                                     y(j,k+1+delta_t/dt)=y(j,k+1+delta_t/dt)+0.5*delta_y_attempt;
                                    x(j,k+1)=x(j,k+1)+0.5*norm([delta_x_attempt,delta_y_attempt])*reflect_normal_unit(1); % the jth particle at k+1+delta_t/dt (being hitted by i) goes forward half i's way
                                    y(j,k+1)=y(j,k+1)+0.5*norm([delta_x_attempt,delta_y_attempt])*reflect_normal_unit(2);
                                end
                            end
                        elseif i<j % j has not been updated (now at k+delta_t/dt), now updating i to k+1+delta_t/dt
                            Diff_x=x(j,k)-x(i,k+1);
                            Diff_y=y(j,k)-y(i,k+1);
                            %                                 if Diff_x^2+Diff_y^2 < (2*a)^2
                            if norm([Diff_x Diff_y]) < 2*a
                                reflect_normal_unit=unit_vec(Diff_x,Diff_y); % The normal vector of the reflective surface
                                reflected_unit=unit_vec(delta_x_attempt,delta_y_attempt)-2*(unit_vec(delta_x_attempt,delta_y_attempt)*reflect_normal_unit')*reflect_normal_unit;
                                
                                %                                     x(i,k+1+delta_t/dt)=x(i,k+1+delta_t/dt)-0.5*delta_x_attempt; % the ith particle at k+1 (hitting j) goes backwards half its way
                                %                                     y(i,k+1+delta_t/dt)=y(i,k+1+delta_t/dt)-0.5*delta_y_attempt;
                                x(i,k+1)=x(i,k+1)-0.5*delta_x_attempt+0.5*norm([delta_x_attempt,delta_y_attempt])*reflected_unit(1); % the ith particle at k+1 (hitting j) goes backwards half its way
                                y(i,k+1)=y(i,k+1)-0.5*delta_y_attempt+0.5*norm([delta_x_attempt,delta_y_attempt])*reflected_unit(2);
                                if fixed_flag(j)==0
                                    %                                     x(j,k+delta_t/dt)=x(j,k+delta_t/dt)    +0.5*delta_x_attempt; % the jth particle at k (being hitted by i) goes forward half i's way
                                    %                                     y(j,k+delta_t/dt)=y(j,k+delta_t/dt)    +0.5*delta_y_attempt;
                                    x(j,k)=x(j,k)    +0.5*norm([delta_x_attempt,delta_y_attempt])*reflect_normal_unit(1); % the jth particle at k (being hitted by i) goes forward half i's way
                                    y(j,k)=y(j,k)    +0.5*norm([delta_x_attempt,delta_y_attempt])*reflect_normal_unit(2);
                                end
                            end
                        end
                    end
            end
            if fixed_flag(i)==1 %% particle is fixed
                x(i,k+1)=x(i,k);
                y(i,k+1)=y(i,k);
            end               
    end
    switch hard_collision
        %% Elastic Collision: This method requires that every particle is updated, and we calculate the reflected path with the dx_attempts(1:N)
        case 'method_3' % Elastic collision
            dx_attempt_vec(1:N,1:2)=0;
            for i=1:N
                dx_attempt_vec(i,:)=[x(i,k+1)-x(i,k),y(i,k+1)-y(i,k)];
%                 abs(dx_attempt_vec(i,1)-delta_x(i))<10^-6;
%                 abs(dx_attempt_vec(i,2)-delta_y(i))<10^-6;
                1;
            end
            
%             while 1
                check_relax_3(1:N,1:N)=0;
                
                for i=1:N
                    for j=1:i
                        if j~=i % both i and j have been updated to k+1, now updating i to k+1
                            diff_x=x(j,k+1)-x(i,k+1);
                            diff_y=y(j,k+1)-y(i,k+1);
                            diff_r_sqr=diff_x^2+diff_y^2;
                            if diff_r_sqr < (2*a-epsilon)^2 % The epsilon is an allowed error for the particles to overlap; if we don't set this, (2*a-sqrt(diff_r_sqr)) will be too increasingly small and the simulation will get stuck at some k.
                                %% Calculating the elastic collision according to Strating (3.7) and (3.8)
                                %                                 v1_unit=unit_vec(dx_attempt_vec(i,1),dx_attempt_vec(i,2));
                                %                                 v2_unit=unit_vec(dx_attempt_vec(j,1),dx_attempt_vec(j,2));
                                dist1=sqrt(diff_r_sqr);
                                dist2=2*a;
                                dist3=norm([x(j,k)-x(i,k),y(j,k)-y(i,k)]);
                                fraction=(dist2-dist1)/(dist3-dist1);
                                v_1=dx_attempt_vec(i,:)/dt;
                                v_2=dx_attempt_vec(j,:)/dt;
                                r_21= - unit_vec(diff_x,diff_y); % would be more accurate if i subtract the overlapping part fraction*dx_attempt_vec(i,:)
                                v_21=v_1-v_2;
                                v_1_reflected=v_1-(r_21*v_21')*r_21;
                                v_2_reflected=v_2+(r_21*v_21')*r_21;
                                
                                %                                     X1_reflected=[x(i,k),y(i,k)]+(1-fraction)*dx_attempt_vec(i,:)+(fraction)*norm(dx_attempt_vec(i,:))*unit_vec(v_1_reflected(1),v_1_reflected(2));
                                %                                     X2_reflected=[x(j,k),y(j,k)]+(1-fraction)*dx_attempt_vec(j,:)+(fraction)*norm(dx_attempt_vec(j,:))*unit_vec(v_2_reflected(1),v_2_reflected(2));
                                X1_reflected=[x(i,k+1),y(i,k+1)]-fraction*dx_attempt_vec(i,:)+(fraction)*norm(dx_attempt_vec(i,:))*unit_vec(v_1_reflected(1),v_1_reflected(2));
                                X2_reflected=[x(j,k+1),y(j,k+1)]-fraction*dx_attempt_vec(j,:)+(fraction)*norm(dx_attempt_vec(j,:))*unit_vec(v_2_reflected(1),v_2_reflected(2));
                                dx_attempt_vec(i,:)=X1_reflected-[x(i,k+1),y(i,k+1)];
                                dx_attempt_vec(j,:)=X2_reflected-[x(j,k+1),y(j,k+1)];
                                x(i,k+1)=X1_reflected(1);
                                y(i,k+1)=X1_reflected(2);
                                x(j,k+1)=X2_reflected(1);
                                y(j,k+1)=X2_reflected(2);
                                %                                     if isnan(x(i,k+1))
                                if ismember(dx_attempt_vec,NaN)
                                    %                                         dx_attempt_vec
                                    %                                         X1_reflected
%                                     %                                         fraction
%                                     
%                                     v_1_reflected
%                                     v_2_reflected
%                                     r_21
%                                     x(i,k)-x(j,k)
%                                     y(i,k)-y(j,k)
%                                     norm([x(i,k)-x(j,k),y(i,k)-y(j,k)])
%                                     k
%                                     
                                    %                                         pause
                                end
                                %                                     %% Note that here we comment out the fixed_flag(i)==0 condition or else the calculations will take forever.
                                %                                     if fixed_flag(i)==0
                                %                                         x(i,k+1)=x(i,k+1)-b*diff_x/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr)); % the ith particle at k+1 (hitting j) goes backwards half its way
                                %                                         y(i,k+1)=y(i,k+1)-b*diff_y/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr));
                                %                                     elseif fixed_flag(i)==1
                                %                                         x(i,1+k)=x(i,k);
                                %                                         y(i,1+k)=y(i,k);
                                %                                     end
                                %                                     if fixed_flag(j)==0
                                %                                         x(j,k+1)=x(j,k+1)+b*diff_x/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr)); % the jth particle at k+1 (being hitted by i) goes forward half i's way
                                %                                         y(j,k+1)=y(j,k+1)+b*diff_y/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr));
                                %                                     elseif fixed_flag(j)==1
                                %                                         x(j,1+k)=x(j,k);
                                %                                         y(j,1+k)=y(j,k);
                                %                                     end
                            else
                                check_relax_3(i,j)=1;
                            end
                        end
                    end
                end
                if sum(sum(check_relax_3,2))==(N^2-N)/2 % All particles have inter distance larger than 2a
                    for i=1:N
                        if fixed_flag(i)==1 %% particle is fixed, overwrite the updated x(i,1+k) with x(i,k)
                            x(i,1+k)=x(i,k);
                            y(i,1+k)=y(i,k);
                        end
                    end
%                     break
                else
                    % Debug tools.
                    
                    %                         diff_r_sqr-(2*a)^2
                    %                         sum(sum(check_relax_3,2))
                    %                         check_relax_3
                    %                         k
                    %                         figure(1);
                    %
                    %
                    %                         for i=1:N
                    %                             hold on
                    %                             plot(x(i,1+k),y(i,1+k),'o')
                    %                         end
                    %                         A=[x(1,1+k)-x(2,1+k),y(1,1+k)-y(2,1+k)];
                    %                         norm(A)
                end
                
%             end
    end
    
    %% Relaxation: What if the particles still overlap after the previous subsection?
    switch hard_collision
        case {'method_1', 'method_2','test_no_elastic','method_3'}
            while 1==1
                check_relax(1:N,1:N)=0;
                for i=1:N-1
                    for j=i+1:N
                        if j~=i % both i and j have been updated to k+1+delta_t/dt, now updating i to k+1+delta_t/dt
                            Diff_x=x(j,k+1)-x(i,k+1);
                            Diff_y=y(j,k+1)-y(i,k+1);
                            diff_r_sqr=Diff_x^2+Diff_y^2;
                            if diff_r_sqr < (2*a-epsilon)^2
                                %                         if fixed_flag(i)==0 % i is mobile
                                x(i,k+1)=x(i,k+1)-b*Diff_x/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr)); % the ith particle at k+1+delta_t/dt (hitting j) goes backwards half its way
                                y(i,k+1)=y(i,k+1)-b*Diff_y/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr));
                                %                         end
                                %                         if fixed_flag(j)==0 % j is mobile; if j==1 then j doesn't move, the line under will not execute
                                x(j,k+1)=x(j,k+1)+b*Diff_x/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr)); % the jth particle at k+1+delta_t/dt (being hitted by i) goes forward half i's way
                                y(j,k+1)=y(j,k+1)+b*Diff_y/sqrt(diff_r_sqr)*(2*a-sqrt(diff_r_sqr));
                                %                         end
                            else
                                check_relax(i,j)=1;
                            end
                        end
                    end
                end
                
                if sum(sum(check_relax,2))==(N^2-N)/2 % All particles have inter distance larger than 2a
                    if fixed_flag(i)==1 %% particle is fixed, overwrite the updated x(i,1+k) with x(i,k)
                        x(i,1+k)=x(i,k);
                        y(i,1+k)=y(i,k);
                    end
                    break
                else
                    sum(sum(check_relax,2))
                    k
                end
            end
    end
end

%% Calculating the velocities at each timestep (from t=0+delta_t ~ Obs_time+delta_t)
v_x=diff(x,1,2)/dt;
v_y=diff(y,1,2)/dt;
time=(1:delta_t/dt)*dt;
end
    
    
    
    
    
