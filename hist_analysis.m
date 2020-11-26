function [k_trans,theta_plus,theta_minus,num_transitions]=hist_analysis(movie_name,moving_avg,num_bins,bin_limit,Obs_time_steps,delta_t,dt)

% moving_avg=1;
load([movie_name,'.mat'],'theta','time')

bin_interval=2*bin_limit/num_bins;
bin_loc=-bin_limit:bin_interval:bin_limit;
%     h=histogram(movmean(theta(2,:),moving_avg),1000)
h=histogram(movmean(theta(2,:),moving_avg),bin_loc);

%%% Note that before reaching steady state, theta=0, and that is not the
%%% part we are interested in.
% theta_plus_index=find(h.Values==max(h.Values(end/2+2:end))); % end/2+2 instead of end/2 to avoid theta=0;
% theta_minus_index=find(h.Values==max(h.Values(1:end/2-2)));  % end/2-2 instead of end/2 to avoid theta=0;
theta_plus_index=find(h.Values==max(h.Values(end/2:end))); % end/2+2 instead of end/2 to avoid theta=0;
theta_minus_index=find(h.Values==max(h.Values(1:end/2)));  % end/2-2 instead of end/2 to avoid theta=0;
if theta_plus_index==num_bins-find(h.Values==max(h.Values(1:end/2)))
    'The histogram is symmetric.'
end

% if length(theta_plus_index)==1 
%     theta_plus=bin_loc(theta_plus_index)+0.5*bin_interval;
% else 
%     theta_plus=0;
%     warning('There is more than one maximum for theta_plus')
% end
% if  length(theta_minus_index)==1
%     theta_minus=bin_loc(theta_minus_index)+0.5*bin_interval;
% else 
%     theta_minus=0;
%     warning('There is more than one maximum for theta_minus')
% end

if length(theta_plus_index)==1 
    theta_plus=bin_loc(theta_plus_index)+0.5*bin_interval;
elseif isempty(theta_plus_index)
    theta_plus=0;
else
    theta_plus=mean(bin_loc(theta_plus_index))+0.5*bin_interval;
    warning('There is more than one maximum for theta_plus')
end
if  length(theta_minus_index)==1
    theta_minus=bin_loc(theta_minus_index)+0.5*bin_interval;
elseif isempty(theta_minus_index)
    theta_minus=0;
else
    theta_minus=mean(bin_loc(theta_minus_index))+0.5*bin_interval;
    warning('There is more than one maximum for theta_minus')
end

%%
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
save([movie_name,'.mat'],'theta_plus','theta_minus','num_transitions','-append')
clear theta
end