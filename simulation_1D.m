clear;
close all
N=10^7;

dt=10^-3
D=2;
omega_0=1;
delta_t=100;
theta(1:1+N)=0; %% theta=theta_delay is omega*delta_t
phi(1:1+N+delta_t/dt)=0; %% theta_polar is the ever-increasing polar angle
theta_0=omega_0*delta_t;
theta_plus=sqrt(6*(theta_0-1)/theta_0);
r_theta=normrnd(0,sqrt(8*D*dt/theta_0^2),[N,1]);
r_phi=normrnd(0,sqrt(2*D*dt),[N,1]);
%% For approximate result
for k=1:N
    theta(k+1)=theta(k)+1/(3*delta_t)*(theta_plus^2-theta(k)^2)*theta(k)*dt+r_theta(k);
end
%% For full simulation
for k=1:N
    phi(1+k+delta_t/dt)=phi(k+delta_t/dt)+omega_0*sin(phi(k+delta_t/dt)-phi(k))*dt+r_phi(k);
end


%% Calculating D_eff for theta
recalculate_T_eff='yes'
plot_hist_fit_T_eff='no'
switch recalculate_T_eff
    case 'yes'

        %% Finding the sigma of the omega (effective temperature)
        f=figure(80);clf;
        switch plot_hist_fit_T_eff
            case 'no'
                f.Visible='off';
        end
        %         h=histogram(diff(theta)/delta_t);
        h=histogram(diff(theta));
        y=h.Values;
        x=h.BinEdges(1:end-1)+0.5*h.BinWidth;
        % plot(x,y);
        [fitresult,gof]=fit_Gaussian(x,y,plot_hist_fit_T_eff);
        % a1*exp(-((x-b1)/c1)^2)
        mu=fitresult.b1;
        sigma=fitresult.c1/sqrt(2);
        
        %% Calculating the D_eff and T_eff
        k_B=1
        D_eff_theta=sigma^2*theta_0^2/(8*dt);
%         T_eff_theta=(omega_0*delta_t^2*sigma^2)/(4*k_B*dt);
        
%         save([movie_name,'.mat'],'D_eff','T_eff','sigma','mu','-append')
    case 'no'
%         load([movie_name,'.mat'],'D_eff','T_eff','sigma','mu')
end


%% Calculating D_eff for phi
switch recalculate_T_eff
    case 'yes'
        %% Finding the sigma of the omega (effective temperature)
        f=figure(80);clf;
        switch plot_hist_fit_T_eff
            case 'no'
                f.Visible='off';
        end
        %         h=histogram(diff(theta)/delta_t);
        h2=histogram(diff(phi));
        y2=h2.Values;
        x2=h2.BinEdges(1:end-1)+0.5*h2.BinWidth;
        % plot(x,y);
        [fitresult,gof]=fit_Gaussian(x2,y2,plot_hist_fit_T_eff);
        % a1*exp(-((x-b1)/c1)^2)
        mu2=fitresult.b1;
        sigma2=fitresult.c1/sqrt(2);
        
        %% Calculating the D_eff and T_eff
        k_B=1
        D_eff_phi=sigma2^2/(2*dt);
%         T_eff_phi=(omega_0*delta_t^2*sigma2^2)/(4*k_B*dt);
        
%         save([movie_name,'.mat'],'D_eff','T_eff','sigma','mu','-append')
    case 'no'
%         load([movie_name,'.mat'],'D_eff','T_eff','sigma','mu')
end

%%
D_eff_theta
D_eff_phi
%%
% 
%% Histogram: determining theta+ and theta-
num_bins=100;
bin_limit=1;
bin_interval=2*bin_limit/num_bins;
bin_loc=-bin_limit:bin_interval:bin_limit;
figure(1)
%% Here the histogram is normalized by theta_0; the obtained theta_plus and theta_minus will be restored of this factor in the end of this section
h=histogram(theta/theta_0,bin_loc);
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
theta_plus=theta_plus*theta_0;
theta_minus=theta_minus*theta_0;
%% Calculating the transitions
sign_old=0; % Initial 'order parameter' for the orbit. +1 for stable orbit with theta_plus, -1 for stable orbit with theta_miunus
sign_current=0;
% sign_current_matrix(1:N+1)=0;
num_transitions=0;
k_trans=[];
for k=1:N+1
    if theta(k)>theta_plus
        sign_current=1;
%         sign_current_matrix(k)=1;
    elseif theta(k)<theta_minus
        sign_current=-1;
%         sign_current_matrix(k)=-1;
    end
    % Updating sign: sign_old -> sign_current
    if sign_current*sign_old== -1
        num_transitions=num_transitions+1;
        k_trans=[k_trans k];
    end
    sign_old=sign_current;
end
%%
close all
% figure(2)
% plot(theta_delay)
figure(3)
plot(1:N+1,theta,'.')
hold on
plot(k_trans,r(k_trans),'o')
% plot(1:N+1,sign_current_matrix,'x')
% plot(k_trans,0,'o')

std(r)
std(r(k_trans-1))
theta_plus
theta_minus
num_transitions





