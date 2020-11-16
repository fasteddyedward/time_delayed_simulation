function Movie_Vector=Make_Movie(movie_name,N,dt,Obs_time_steps,partition_time_steps,delta_t,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,partition_movie,movie_x_min, movie_x_max,movie_y_min,movie_y_max,a,force_tracks)
switch partition_movie
    case 'no'
        %% Start plotting the movies and plots
        load([movie_name,'.mat'],'time','x','y')
        Movie_Vector=make_movies_plots(N,delta_t,dt,Obs_time_steps+delta_t/dt,x,y,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,a,force_tracks);
        clear time x y
    case 'yes'
        %% Start plotting the movies and plots
        close all
        Movie_Vector=[];
        %             v_omega=[];
        axis_scale=[movie_x_min movie_x_max movie_y_min movie_y_max];
        for lth_partition=1:round(Obs_time_steps/partition_time_steps)+1
            if lth_partition==1
                %% Start plotting movie for stage 1
                cd(movie_name)
                load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
                cd ..
                time=(1:delta_t/dt)*dt;
                Movie_Vector_Partition=make_movies_plots(N,delta_t,dt,delta_t/dt,x,y,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,a,force_tracks);
                clear x y v_x v_y time
            else
                %% Start plotting movie for stage 2
                cd(movie_name)
                ['lth_partition= ',num2str(lth_partition),' out of ',num2str(round(Obs_time_steps/partition_time_steps)+1)]
                load([movie_name,' partition_',num2str(lth_partition),'.mat'],'x','y','v_x','v_y')
                cd ..
                time=(1+(lth_partition-2)*partition_time_steps+delta_t/dt:(lth_partition-1)*partition_time_steps+delta_t/dt)*dt;
                Movie_Vector_Partition=make_movies_plots(N,delta_t,dt,partition_time_steps,x,y,time,magnify,control_animation_interval,movie_create,ghost,axis_choice,leave_trace,axis_scale,a,force_tracks);
                clear x y v_x v_y time
            end
            Movie_Vector=[Movie_Vector Movie_Vector_Partition];
        end
        
        
end
