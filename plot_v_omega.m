function plot_v_omega(N,delta_t,movie_name,moving_avg)

load([movie_name,'.mat'],'v_omega','time')
hold on
for i=1:N
    plot(time,movmean(v_omega(i,:),moving_avg))
end
xline(delta_t)
title('Normalized Rotation Speed')
xlabel('time (ms)')
ylabel('v_\omega(rad/ms) ')
legend('1','2','3','\deltat')        