function plot_theta(N,delta_t,movie_name,moving_avg,theta_plus,theta_minus,k_trans,show_transitions)

load([movie_name,'.mat'],'theta','time')
hold on
for i=2:N
    plot(time,movmean(theta(i,1:length(time)),moving_avg))
end
xline(delta_t)
yline(theta_plus)
yline(theta_minus)
title('\theta')
xlabel('time (ms)')
ylabel('\theta (rad/ms) ')
legend('particle','\deltat',['\theta_+=',num2str(theta_plus)],['\theta_-=',num2str(theta_minus)])    
switch show_transitions
    case 'on'
        for k=1:length(k_trans)
            xline(time(k_trans(k)),'r')
        end
end
clear theta time