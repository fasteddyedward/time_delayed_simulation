function Save_created_movie(movie_name,movie_create,control_animation_interval,Obs_time_steps,delta_t,dt,frame_rate,Movie_Vector)
switch movie_create
    case 'on'
        if control_animation_interval<=Obs_time_steps+dt*delta_t
%             frame_rate=10    ;
            if exist('Movie_Vector','var')==0
                load([movie_name,'.mat'],'Movie_Vector')
            end
            save_movie(Movie_Vector(1:end),movie_name,frame_rate);
        end
end
end