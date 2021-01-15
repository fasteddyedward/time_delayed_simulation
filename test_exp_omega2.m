clear;
close all
load('dt=1.500.tdms_angle.mat')
omega_0=2.1/R_mean(2);
theta_0=1.5*omega_0;
    D=0.065; % Measured from Exp_Free_Diffusion.m from Xiangzun's old diffusive particle
% omega2=omega'-omega_0*sin(theta(2,:));
% plot(omega2)
% figure;
% plot(omega*1.5)
% hold on
% plot(theta(2,:))

% figure;
% plot(theta(2,:)./omega')

%%
theta=theta(2,:)';

%% 2021.1.12 hist(omega-\pm omega_+)

[omega_stable,k_trans_omega,omega_plus,omega_minus,num_transitions_omega]=find_omega_plus(omega,theta_0);

%         theta_plus=sqrt(10-sqrt(1/(42*theta_0)+120/theta_0-20));
%         theta_minus=-sqrt(10-sqrt(1/(42*theta_0)+120/theta_0-20));

% omega3=omega-omega_stable;
%         omega_0=v_0/mean(R);
% omega_0=2.1/R_mean(2);
omega2=omega-omega_0*sin(theta);

h_omega2=histogram(omega2);
Values_omega2=h_omega2.Values;
Bins_omega2=h_omega2.BinEdges(1:end-1)+0.5*h_omega2.BinWidth;

%         Values_omega2(abs(Bins_omega2)<0.5)=[];
%         Bins_omega2(abs(Bins_omega2)<0.5)=[];

% [fitresult,gof]=fit_Gaussian(Bins_omega2,Values_omega2,plot_hist_fit_T_eff);
        [fitresult,gof]=fit_Gaussian(Bins_omega2,Values_omega2,'yes');
title('hist omega')
% a1*exp(-((x-b1)/c1)^2)
%         histogram(omega2)
% histogram(omega2)
gof
theta_0
mu_omega2=fitresult.b1;
sigma_omega2=fitresult.c1/sqrt(2);
D_omega2=sigma_omega2^2*dt/(2);
D_ratio=D_omega2/(D/R_mean(2)^2)
        %% Calculating omega_+ and k_trans with hist(omega)
function [omega_stable,k_trans_omega,omega_plus,omega_minus,num_transitions]=find_omega_plus(omega,theta_0)
% omega_stable is omega_+-
omega_stable(1:length(omega))=0;
if theta_0<0.9
    omega_stable(1:length(omega))=0;
    k_trans_omega=[];
    omega_plus=0;
    omega_minus=0;
    num_transitions=0;
else
    h=histogram(omega);
    Values=h.Values;
    Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    omega_plus_index=find(Values==max(Values(end/2:end))); % end/2+2 instead of end/2 to avoid theta=0;
    omega_minus_index=find(Values==max(Values(1:end/2)));  % end/2-2 instead of end/2 to avoid theta=0;
    if length(omega_plus_index)==1
        omega_plus=Bins(omega_plus_index);
    elseif isempty(omega_plus_index)
        omega_plus=0;
    else
        omega_plus=mean(Bins(omega_plus_index));
        warning('There is more than one maximum for theta_plus')
    end
    
    if  length(omega_minus_index)==1
        omega_minus=Bins(omega_minus_index);
    elseif isempty(omega_minus_index)
        omega_minus=0;
    else
        omega_minus=mean(Bins(omega_minus_index));
        warning('There is more than one maximum for theta_minus')
    end   
    %% Start Calculating Transition Rates
    sign_old=0; % Initial 'order parameter' for the orbit. +1 for stable orbit with theta_plus, -1 for stable orbit with theta_miunus
    sign_current=0;
    num_transitions=0;
    k_trans_omega=[];
    for k=1:length(omega)
        if omega(k)>omega_plus
            sign_current=1;
%             theta_stable(k)=theta_plus;
        elseif omega(k)<omega_minus
            sign_current=-1;
%             theta_stable(k)=theta_minus;
        end
        % Updating sign: sign_old -> sign_current
        if sign_current*sign_old==-1
            num_transitions=num_transitions+1;
            k_trans_omega=[k_trans_omega k];
        end
        
        if sign_current==1
            omega_stable(k)=omega_plus;
        elseif sign_current==-1
            omega_stable(k)=omega_minus;
        end
        
        sign_old=sign_current;
    end
    
end
end