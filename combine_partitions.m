function [x_full_time,y_full_time,v_x_full_time,v_y_full_time,time]=combine_partitions(movie_name,Obs_time_steps,partition_time_steps,delta_t,dt)
x_full_time=[];
y_full_time=[];
v_x_full_time=[];
v_y_full_time=[];
for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
    if lth_partition==1
        %% Putting the data for the first stage: pure diffusion
        load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
        x_full_time=[x_full_time x(:,1:end)];
        y_full_time=[y_full_time y(:,1:end)];
        v_x_full_time=[v_x_full_time v_x(:,1:end)];
        v_y_full_time=[v_y_full_time v_y(:,1:end)];
        movefile([movie_name,' partition_',num2str(lth_partition),'.mat'],movie_name);
    else
        %% Putting the data for the second stage: interactions
        load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
        x_full_time=[x_full_time x(:,2:end)];
        y_full_time=[y_full_time y(:,2:end)];
        v_x_full_time=[v_x_full_time v_x(:,1:end)];
        v_y_full_time=[v_y_full_time v_y(:,1:end)];
        movefile([movie_name,' partition_',num2str(lth_partition),'.mat'],movie_name);
    end
end
time=(1:Obs_time_steps+delta_t/dt)*dt;
end