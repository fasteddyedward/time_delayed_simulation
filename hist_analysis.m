function [k_trans,theta_plus,theta_minus,num_transitions]=hist_analysis(movie_name,moving_avg,num_bins,bin_limit,Obs_time_steps,delta_t,dt)

moving_avg=1;
load([movie_name,'.mat'],'theta','time')

bin_interval=2*bin_limit/num_bins;
bin_loc=-bin_limit:bin_interval:bin_limit;
%     h=histogram(movmean(theta(2,:),moving_avg),1000)
h=histogram(movmean(theta(2,:),moving_avg),bin_loc);
theta_plus_index=find(h.Values==max(h.Values(end/2:end)));
theta_minus_index=find(h.Values==max(h.Values(1:end/2)));
if theta_plus_index==num_bins-find(h.Values==max(h.Values(1:end/2)))
    'The histogram is symmetric.'
end
theta_plus=bin_loc(theta_plus_index)+0.5*bin_interval;
theta_minus=bin_loc(theta_minus_index)+0.5*bin_interval;
hold on
xline(theta_plus)
xline(theta_minus)
legend('Entries',['\theta_+=',num2str(theta_plus)],['\theta_-=',num2str(theta_minus)])
title('Histogram of \theta')
xlabel('\theta (rad)')
ylabel('Entries')

%% Start Calculating Transition Rates
theta_sing=theta(2,:);


sign_old=0; % Initial 'order parameter' for the orbit. +1 for stable orbit with theta_plus, -1 for stable orbit with theta_miunus
sign_current=0;
num_transitions=0;
k_trans=[];
for k=1:Obs_time_steps+delta_t/dt
    if theta_sing(k)>theta_plus
        sign_current=1;
    elseif theta_sing(k)<theta_minus
        sign_current=-1;
    end
    % Updating sign: sign_old -> sign_current
    if sign_current*sign_old==-1
        num_transitions=num_transitions+1;
        k_trans=[k_trans k];
    end
    sign_old=sign_current;
    
end
clear theta

%% Fitting with polynomial (not used at the moment)
% bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
% f=polyfit(bins,h.Values,6);
% fitted_hist=polyval(f,bins);
% figure(82);clf;
% plot(bins,h.Values)
% hold on
% plot(bins,fitted_hist)
% %     f=polyfit(delta_t(flag_not_NaN),log(transition_rate(flag_not_NaN)),2);
% h.Values(end/2+1)
end