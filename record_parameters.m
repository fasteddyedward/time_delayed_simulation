function record_parameters(movie_name,delta_t,v_0,T,N,Obs_time_steps,a,b)
fid=fopen([movie_name,'_parameters.txt']);
fprint('delta_t = ',num2str(delta_t),'\n')
fprint('v_0 = ',num2str(v_0),'\n')
fprint('T = ',num2str(T),'\n')
fprint('\n')
fprint('N = ',num2str(N),'\n')
fprint('Obs_time_steps = ',num2str(Obs_time_steps),'\n')
fprint('a = ',num2str(a),'\n')
fprint('b = ',num2str(b),'\n')


fclose(fid);
end