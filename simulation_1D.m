clear;
close all
N=10^5;
theta(1:N+1)=0;
dt=10^-1;
D=0.5;
omega_0=10;
delta_t=1;
theta_0=omega_0*delta_t;
theta_plus=sqrt(6*(theta_0-1)/theta_0);
r=normrnd(0,sqrt(2*D*dt),[N,1]);
%% For approximate result
for k=1:N
    theta(k+1)=theta(k)+1/(3*delta_t)*(theta_plus^2-theta(k)^2)*theta(k)*dt+r(k);
end
%% For full simulation
for k=delta_t/dt:N
    theta(k+1)=theta(k)-
% plot(theta)
% histogram(theta,100)

%% Histogram: determining theta+ and theta-
num_bins=100;
bin_limit=1
bin_interval=2*bin_limit/num_bins;
bin_loc=-bin_limit:bin_interval:bin_limit;
figure(1)
h=histogram(theta/theta_0,bin_loc)
theta_plus_index=find(h.Values==max(h.Values(end/2:end))); % end/2+2 instead of end/2 to avoid theta=0;
theta_minus_index=find(h.Values==max(h.Values(1:end/2)));  % end/2-2 instead of end/2 to avoid theta=0;
if theta_plus_index==num_bins-find(h.Values==max(h.Values(1:end/2)))
    'The histogram is symmetric.'
end
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
%% Calculating the transitions
sign_old=0; % Initial 'order parameter' for the orbit. +1 for stable orbit with theta_plus, -1 for stable orbit with theta_miunus
sign_current=0;
num_transitions=0;
k_trans=[];
for k=1:N
    if theta(k)>theta_plus
        sign_current=1;
    elseif theta(k)<theta_minus
        sign_current=-1;
    end
    % Updating sign: sign_old -> sign_current
    if sign_current*sign_old==-1
        num_transitions=num_transitions+1;
        k_trans=[k_trans k];
    end
    sign_old=sign_current;
end
%%
close all
figure(2)
plot(theta)
figure(3)
plot(1:N+1,theta)
hold on
plot(k_trans,r(k_trans-1),'o')

std(r)
std(r(k_trans-1))



