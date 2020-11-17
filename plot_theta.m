function plot_theta(N,delta_t,movie_name,moving_avg)

load([movie_name,'.mat'],'theta','time')
hold on
for i=2:N
    plot(time,movmean(theta(i,:),moving_avg))
end
xline(delta_t)
title('\theta')
xlabel('time (ms)')
ylabel('\theta (rad/ms) ')
legend('1','2','3','\deltat')        