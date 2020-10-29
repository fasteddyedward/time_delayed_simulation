% lower=78
% for nth_take=lower:154
    nth_take=123
    moving_avg=100000
    movie_name=['2020.10.27,dt=10e-3 take ',num2str(nth_take)];
%     if (exist([movie_name,'.mat'],'file')==0)
%         continue
%     end

    load([movie_name,'.mat'],'delta_t','T','v_0');
    
    N=3;
    figure(98); clf %% Would be same as figure(99) if partition movie='no'
    plot_v_omega(N,delta_t,movie_name,moving_avg)
    title(['Normalized Rotation Speed, v_0 = ',num2str(v_0),', \delta t = ',num2str(delta_t),', T = ',num2str(T),', moving_{avg}= ',num2str(moving_avg)])
    % subtitle(['v_0=',num2str(v_0),', delta_t=',num2str(delta_t),', T=',num2str(T)])
    saveas(gcf,[movie_name,'.png'])
% end1