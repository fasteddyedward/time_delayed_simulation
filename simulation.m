function [x,y,v_x,v_y,movie_x_min,movie_x_max,movie_y_min,movie_y_max,x_final,y_final]=simulation(movie_name,Obs_time_steps,partition_time_steps,movie_x_min,movie_x_max,movie_y_min,movie_y_max,N,delta_t,dt,v_0,T,gamma,k_B,D,x_init,y_init,hard_collision,a,b,int_delay,fixed_flag)
%% First stage: Diffusion. t=0 ~ delta_t (lth_partition=1)
    % Obs_time_Steps in this case is round(delta_t/dt)
    ['lth_partition= ',num2str(1),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
    [x,y,v_x,v_y,~]=first_stage_pure_diffusion(N,delta_t,dt,T,x_init,y_init,gamma,k_B,D,hard_collision,a,b,fixed_flag);
    movie_x_max=max([movie_x_max max(x)]);    movie_x_min=min([movie_x_min min(x)]);    movie_y_max=max([movie_y_max max(y)]);    movie_y_min=min([movie_y_min min(y)]);
    save([movie_name,' partition_',num2str(1),'.mat']); %% Takes ~ 0.16 sec
    x_temp=x;
    y_temp=y;
    %% Second stage: Delayed interaction starts. t=delta_t~Obs_time
    for lth_partition=2:round(Obs_time_steps/partition_time_steps)+1
        ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
        [x,y,~,~,v_x,v_y,~,~,~]=second_stage_delayed_int(N,delta_t,dt,partition_time_steps,v_0,T,x_temp,y_temp,lth_partition,gamma,k_B,D,hard_collision,a,b,int_delay,fixed_flag);
        movie_x_max=max([movie_x_max max(x)]);        movie_x_min=min([movie_x_min min(x)]);        movie_y_max=max([movie_y_max max(y)]);        movie_y_min=min([movie_y_min min(y)]);
        save([movie_name,' partition_',num2str(lth_partition),'.mat']);
        x_temp=x(:,end-delta_t/dt:end); % This last postitions of delayed time delta_t goes to the next round for simulation
        y_temp=y(:,end-delta_t/dt:end);
        clear time
    end
    x_final=x_temp; % We can use this for further simulation
    y_final=y_temp;
    clear x_temp y_temp
end