function Rotational_Analysis(movie_name,partition_movie,N,v_0,Obs_time_steps,partition_time_steps,delta_t,dt,moving_avg,plot_rot)
figure(99)
clf
switch partition_movie
    case 'no'
        load([movie_name,'.mat'],'time','x','y','v_x','v_y')
        v_omega=Rotation(N,x,y,v_0,v_x,v_y,time,Obs_time_steps+delta_t/dt,delta_t,partition_movie,moving_avg,plot_rot);
        %% Saving v_omega
        save([movie_name,'.mat'],'v_omega','-append')
        clear v_omega x y v_x v_y time
    case 'yes'
        v_omega=[];
        for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
            if lth_partition==1
                %% Start plotting v_omega for stage 1
                cd(movie_name)
                ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
                load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
                cd ..
                time=(1:delta_t/dt)*dt;
                v_omega_partition=Rotation(N,x,y,v_0,v_x,v_y,time,delta_t/dt,delta_t,partition_movie,moving_avg,plot_rot);
                clear x y v_x v_y time
            else
                %% Start plotting movie for stage 2
                cd(movie_name)
                ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
                load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
                cd ..
                time=(1+(lth_partition-2)*partition_time_steps+delta_t/dt:(lth_partition-1)*partition_time_steps+delta_t/dt)*dt;
                v_omega_partition=Rotation(N,x,y,v_0,v_x,v_y,time,partition_time_steps,delta_t,partition_movie,moving_avg,plot_rot);
                clear x y v_x v_y time
            end
            % Note appending v_omega like this might cause memory overflow.
            v_omega=[v_omega v_omega_partition];
            %v_omega=0;
        end
        %% Saving v_omega
        save([movie_name,'.mat'],'v_omega','-append')
        clear v_omega
end