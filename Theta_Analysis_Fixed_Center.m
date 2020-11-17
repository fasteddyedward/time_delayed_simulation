function Theta_Analysis_Fixed_Center(movie_name,partition_movie,N,v_0,Obs_time_steps,partition_time_steps,delta_t,dt,moving_avg,plot_rot)
figure(99)
clf
theta(1:Obs_time_steps)=0;
R1=0;
R2=0;
diff_x=[0,0];
theta(1:N,1:Obs_time_steps+delta_t/dt)=0;
switch partition_movie
    case 'no'
        load([movie_name,'.mat'],'x','y')
%         v_omega=Rotation(N,x,y,v_0,v_x,v_y,time,Obs_time_steps+delta_t/dt,delta_t,partition_movie,moving_avg,plot_rot);
        for k=1:Obs_time_steps
            for i=2:N
                R1=[x(i,k)-x(1,k),y(i,k)-y(1,k)];
                R2=[x(i,k+delta_t/dt)-x(1,k+delta_t/dt),y(i,k+delta_t/dt)-y(1,k+delta_t/dt)];
%                 norm(R1)
%                 norm(R2)
                diff_x=[x(i,k+delta_t/dt)-x(i,k),y(i,k+delta_t/dt)-y(i,k)];
                theta(i,k+delta_t/dt)=acos(((norm(R1)^2+norm(R2)^2-norm(diff_x)^2)/(2*norm(R1)*norm(R2))));
            end
        end

        %% Saving theta
        save([movie_name,'.mat'],'theta','-append')
        clear theta x y v_x v_y time
        
        

    case 'yes'
        warning('Sorry I have not written for the yes option.')
%         v_omega=[];
%         for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
%             if lth_partition==1
%                 %% Start plotting v_omega for stage 1
%                 cd(movie_name)
%                 ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
%                 load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
%                 cd ..
%                 time=(1:delta_t/dt)*dt;
%                 v_omega_partition=Rotation(N,x,y,v_0,v_x,v_y,time,delta_t/dt,delta_t,partition_movie,moving_avg,plot_rot);
%                 clear x y v_x v_y time
%             else
%                 %% Start plotting movie for stage 2
%                 cd(movie_name)
%                 ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
%                 load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
%                 cd ..
%                 time=(1+(lth_partition-2)*partition_time_steps+delta_t/dt:(lth_partition-1)*partition_time_steps+delta_t/dt)*dt;
%                 v_omega_partition=Rotation(N,x,y,v_0,v_x,v_y,time,partition_time_steps,delta_t,partition_movie,moving_avg,plot_rot);
%                 clear x y v_x v_y time
%             end
%             % Note appending v_omega like this might cause memory overflow.
%             v_omega=[v_omega v_omega_partition];
%             %v_omega=0;
%         end
%         %% Saving v_omega
%         save([movie_name,'.mat'],'v_omega','-append')
%         clear v_omega
end
end