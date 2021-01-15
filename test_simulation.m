%% This file tests if the pure diffusion case gives the correct diffusion constant
%% For 2 particles purely diffusing with D=1, the diffusive constant measured is correct and gave 1.
clear
close all

% Date='2021.1.8' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=1 % for a=5, particle 1 fixed
% nth_take=100 % for a=5, particle 1 not fixed
% nth_take=200 % for a=0, particle 1 not fixed
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[4:0.5:10]
% dt=10^-3
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6
% 
% Date='2021.1.8' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=300; % for a=5, particle 1 fixed
% % nth_take=400 % for a=5, particle 1 not fixed
% % nth_take=500 % for a=0, particle 1 not fixed
% delta_t_matrix=2;
% T_matrix=[1];
% v_0_matrix=[0.5:0.5:20];
% dt=10^-1;
% intrinsic_delay=0.0; % Intrinsic delay
% Obs_time_steps=10^5;

Date='2021.1.8' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
nth_take=400 % for a=5, particle 1 fixed
delta_t_matrix=[1.8:0.1:4.0]
T_matrix=[1]
v_0_matrix=5
dt=10^-1
intrinsic_delay=0.0 % Intrinsic delay
Obs_time_steps=10^4

% Date='2021.1.8' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=500 % for a=5, particle 1 fixed
% delta_t_matrix=[1.8:0.1:4.0]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-2
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date='2021.1.6'
% nth_take=100
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[0.5:0.5:14]
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5

% Date='2021.1.12' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=400 % for a=5, particle 1 fixed
% delta_t_matrix=[0.5:0.5:4]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-3
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6


D_eff_ratio_matrix_approx=[];
D_eff_ratio_matrix_full=[];
D_eff_ratio_matrix_theta2=[];
D_eff_ratio_matrix_omega2=[];
theta_0_matrix=[];
for delta_t_index=1:length(delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(v_0_matrix)
%% Input File Name
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% [movie_name,'.mat'];
load([movie_name,'.mat'])

%% Calculating angle
theta(1:Obs_time_steps)=0;
omega(1:Obs_time_steps)=0;
R(1:Obs_time_steps)=0;

%%
v_x=diff(x(2,:))/dt;
v_y=diff(y(2,:))/dt;
v_x=[v_x 0]; % to keep v_x same length as x
v_y=[v_y 0]; % to keep v_y same length as y
% v=[v_x; v_y ;zeros(1,length(v_x))];
%% Old method (correct but slower)
% for k=1:Obs_time_steps
% % for k=1+round(delta_t/dt):Obs_time_steps+(delta_t/dt)
%     %% Method one
% %     R1=[x(2,k)-x(1,k),y(2,k)-y(1,k),0];
% %     R2=[x(2,k+round(delta_t/dt))-x(1,k+round(delta_t/dt)),y(2,k+round(delta_t/dt))-y(1,k+round(delta_t/dt)),0];
% %     R1_cross_R2=cross(R1,R2);
% %     theta(k)=asin(R1_cross_R2(3)/(norm(R1)*norm(R2)));
%     %% Method two
%             R1=[x(2,k)-x(1,k),y(2,k)-y(1,k),0];
%             R2=[x(2,k+round(delta_t/dt))-x(1,k+round(delta_t/dt)),y(2,k+round(delta_t/dt))-y(1,k+round(delta_t/dt)),0];
%             diff_x=[x(2,k+round(delta_t/dt))-x(2,k),y(2,k+round(delta_t/dt))-y(2,k),0];
%             %% Determining Sign of theta
%             R1_cross_diff_x=cross(R1,diff_x);
% %             theta(k)=sign(R1_cross_diff_x(3))*acos(((norm(R1)^2+norm(R2)^2-norm(diff_x)^2)/(2*norm(R1)*norm(R2))));
%             theta(k)=sign(R1_cross_diff_x(3))*acos(dot(R1,R2)/(norm(R1)*norm(R2)));
%  %%   
%     if fixed_flag(1)==1
%         R(k)=norm(R1);
%     elseif fixed_flag(1)==0
%         R(k)=norm(R1)/2;
%     end
%     %% Calculating omega
%     omega(k)=(R2(1)*v_y(k+round(delta_t/dt))-R2(2)*v_x(k+round(delta_t/dt)))/norm(R2)^2;
% end

%% 2021.1.15 Vectorized Version
R1=[x(2,1:Obs_time_steps)-x(1,1:Obs_time_steps);y(2,1:Obs_time_steps)-y(1,1:Obs_time_steps);zeros(1,Obs_time_steps)];
R2=[x(2,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt))-x(1,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt));y(2,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt))-y(1,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt));zeros(1,Obs_time_steps)];
diff_x=[x(2,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt))-x(2,1:Obs_time_steps);
    y(2,1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt))-y(2,1:Obs_time_steps);
    zeros(1,Obs_time_steps)];
%             R1_cross_diff_x=cross(R1,diff_x);
R1_cross_diff_x_sign=sign( R1(1,:).*diff_x(2,:)-diff_x(1,:).*R1(2,:));
theta=R1_cross_diff_x_sign.*acos(dot(R1,R2,1)./(vecnorm(R1,2,1).*vecnorm(R2,2,1)));
if fixed_flag(1)==1
    R=vecnorm(R1,2,1);
elseif fixed_flag(1)==0
    R=vecnorm(R1,2,1)/2;
end

    %% Calculating omega
omega=(R2(1,:).*v_y(1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt))-R2(2,:).*v_x(1+round(delta_t/dt):Obs_time_steps+round(delta_t/dt)))./vecnorm(R2,2,1).^2;
% figure(1);histogram(theta);
% figure(2);histogram(omega);
% figure(3);histogram(diff(omega));
%%
% figure
% plot(omega./theta)
% figure(19);clf;hold on; plot(omega*delta_t);plot(theta)
% figure(20);clf;plot(omega);
% figure(21);clf;plot(theta);
%                             figure(1);clf;histogram(omega*delta_t);
%                             figure(2);clf;histogram(theta);
%                             figure(3);clf;plot(omega*delta_t)
%                             figure(4);clf;plot(theta);
%% Calculating D_eff
% close all
% recalculate_T_eff='yes';
    plot_hist_fit_T_eff='no';
% notation='theta';
% switch recalculate_T_eff
%     case 'yes'
        %% Finding the sigma of the omega (effective temperature)
        f=figure(80);clf;
        switch plot_hist_fit_T_eff
            case 'no'
                f.Visible='off';
        end
        %         switch notation
        %             case 'omega' % This is for measuring sigma with histogram(diff(omega))
        %                 h=histogram(diff(theta(2,:))/delta_t);
        %             case 'theta' % This is for measuring sigma with histogram(diff(theta))
        
        h_theta=histogram(diff(theta));
        Values_theta=h_theta.Values;
        Bins_theta=h_theta.BinEdges(1:end-1)+0.5*h_theta.BinWidth;
        
        
        theta_0=v_0*delta_t/mean(R);
%         if theta_0<=1
%             h_phi=histogram(abs(theta-theta_0*sin(theta)));
%         else
%             h_phi=histogram(abs(theta)-theta_0*sin(abs(theta)));
%         end

%         Values_phi=h_phi.Values;
%         Bins_phi=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
%         Values_phi(Bins_phi<0.08)=[];
%         Bins_phi(Bins_phi<0.08)=[];
%         Bins_phi=[-flip(Bins_phi) Bins_phi];
%         Values_phi=[flip(Values_phi) Values_phi];

%         [fitresult,gof]=fit_Gaussian(Bins_phi,Values_phi,plot_hist_fit_T_eff);
%         mu_phi=fitresult.b1;
%         sigma_phi=fitresult.c1/sqrt(2);
        
        
        [fitresult,gof]=fit_Gaussian(Bins_theta,Values_theta,plot_hist_fit_T_eff);
        mu_theta=fitresult.b1;
        sigma_theta=fitresult.c1/sqrt(2);
        %% 2021.1.12 hist(theta-\pm theta_+)
        
%         [theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0);
%         theta2(1:length(theta))=theta-theta_stable;
%         h_theta2=histogram(theta2);
%         Values_theta2=h_theta2.Values;
%         Bins_theta2=h_theta2.BinEdges(1:end-1)+0.5*h_theta2.BinWidth;
%         [fitresult,gof]=fit_Gaussian(Bins_theta2,Values_theta2,plot_hist_fit_T_eff);
%         % a1*exp(-((x-b1)/c1)^2)
%         mu_theta2=fitresult.b1;
%         sigma_theta2=fitresult.c1/sqrt(2);
        %% 2021.1.12 hist(omega-\pm omega_+)
        
        [omega_stable,k_trans_omega,omega_plus,omega_minus,num_transitions_omega]=find_omega_plus(omega,theta_0);

        omega_0=v_0/mean(R);
        omega2=omega-omega_0*sin(theta);
        
        h_omega2=histogram(omega2);
        Values_omega2=h_omega2.Values;
        Bins_omega2=h_omega2.BinEdges(1:end-1)+0.5*h_omega2.BinWidth;
        
        [fitresult,gof]=fit_Gaussian(Bins_omega2,Values_omega2,plot_hist_fit_T_eff);
        % a1*exp(-((x-b1)/c1)^2)
        mu_omega2=fitresult.b1;
        sigma_omega2=fitresult.c1/sqrt(2);
        %%

        %% Calculating the D_eff and T_eff
%         D_phi=sigma_phi^2/(2);
%         T_phi=D_phi;
        
        D_theta=sigma_theta^2/(2*dt);
%         T_theta=D_theta;

%         D_theta2=sigma_theta2^2/(2);
        
        D_omega2=sigma_omega2^2*dt/(2);
%     case 'no'
%         load([movie_name,'.mat'],'D_eff','T_eff','sigma','mu')
% end

theta_0_matrix=[theta_0_matrix theta_0];
% D_eff_ratio_matrix_full=[D_eff_ratio_matrix_full D_phi/(D*delta_t^2/(mean(R)^2))];
D_eff_ratio_matrix_approx=[D_eff_ratio_matrix_approx D_theta/(4*D/(theta_0^2*mean(R)^2))];
% D_eff_ratio_matrix_theta2=[D_eff_ratio_matrix_theta2 D_theta2/(D*delta_t^2/(mean(R)^2))];
D_eff_ratio_matrix_omega2=[D_eff_ratio_matrix_omega2 D_omega2/(D/mean(R)^2)];
%             D_phi
            nth_take=nth_take+1;
        end
        nth_take=nth_take+1;
    end
    nth_take=nth_take+1;
end
%%

% close all
% save(['nth_take=',num2str(nth_take),'.mat']);
% nth_take

%% 2021.1.12 combining test_simulation_cont.m

close all


figure(1);clf;
% plot(theta_0_matrix,D_eff_ratio_matrix_full)
figure(2);clf;
plot(theta_0_matrix,D_eff_ratio_matrix_approx)
figure(3);clf;
% plot(theta_0_matrix,D_eff_ratio_matrix_theta2)
figure(4);clf;
plot(theta_0_matrix,D_eff_ratio_matrix_omega2)
axis([-inf inf 0 inf])
% close all
save('theoretical_D_eff_ratio.mat')


%% Calculating theta_+ and k_trans with hist(theta)
% function [theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0)
% % theta_stable is theta_+-
% theta_stable(1:length(theta))=0;
% if theta_0<0.9
%     theta_stable=0;
%     k_trans_theta=[];
%     theta_plus=0;
%     theta_minus=0;
%     num_transitions_theta=0;
% else
%     h=histogram(theta);
%     Values=h.Values;
%     Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
%     theta_plus_index=find(Values==max(Values(end/2:end))); % end/2+2 instead of end/2 to avoid theta=0;
%     theta_minus_index=find(Values==max(Values(1:end/2)));  % end/2-2 instead of end/2 to avoid theta=0;
%     if length(theta_plus_index)==1
%         theta_plus=Bins(theta_plus_index);
%     elseif isempty(theta_plus_index)
%         theta_plus=0;
%     else
%         theta_plus=mean(Bins(theta_plus_index));
%         warning('There is more than one maximum for theta_plus')
%     end
%     
%     if  length(theta_minus_index)==1
%         theta_minus=Bins(theta_minus_index);
%     elseif isempty(theta_minus_index)
%         theta_minus=0;
%     else
%         theta_minus=mean(Bins(theta_minus_index));
%         warning('There is more than one maximum for theta_minus')
%     end   
%     %% Start Calculating Transition Rates
%     sign_old=0; % Initial 'order parameter' for the orbit. +1 for stable orbit with theta_plus, -1 for stable orbit with theta_miunus
%     sign_current=0;
%     num_transitions_theta=0;
%     k_trans_theta=[];
%     for k=1:length(theta)
%         if theta(k)>theta_plus
%             sign_current=1;
% %             theta_stable(k)=theta_plus;
%         elseif theta(k)<theta_minus
%             sign_current=-1;
% %             theta_stable(k)=theta_minus;
%         end
%         % Updating sign: sign_old -> sign_current
%         if sign_current*sign_old==-1
%             num_transitions_theta=num_transitions_theta+1;
%             k_trans_theta=[k_trans_theta k];
%         end
%         
%         if sign_current==1
%             theta_stable(k)=theta_plus;
%         elseif sign_current==-1
%             theta_stable(k)=theta_minus;
%         end
%         
%         sign_old=sign_current;
%     end
%     
% end
% end
%% Calculating omega_+ and k_trans with hist(omega)
% function [omega_stable,k_trans_omega,omega_plus,omega_minus,num_transitions]=find_omega_plus(omega,theta_0)
% % omega_stable is omega_+-
% omega_stable(1:length(omega))=0;
% if theta_0<0.9
%     omega_stable(1:length(omega))=0;
%     k_trans_omega=[];
%     omega_plus=0;
%     omega_minus=0;
%     num_transitions=0;
% else
%     h=histogram(omega);
%     Values=h.Values;
%     Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
%     omega_plus_index=find(Values==max(Values(end/2:end))); % end/2+2 instead of end/2 to avoid theta=0;
%     omega_minus_index=find(Values==max(Values(1:end/2)));  % end/2-2 instead of end/2 to avoid theta=0;
%     if length(omega_plus_index)==1
%         omega_plus=Bins(omega_plus_index);
%     elseif isempty(omega_plus_index)
%         omega_plus=0;
%     else
%         omega_plus=mean(Bins(omega_plus_index));
%         warning('There is more than one maximum for theta_plus')
%     end
%     
%     if  length(omega_minus_index)==1
%         omega_minus=Bins(omega_minus_index);
%     elseif isempty(omega_minus_index)
%         omega_minus=0;
%     else
%         omega_minus=mean(Bins(omega_minus_index));
%         warning('There is more than one maximum for theta_minus')
%     end   
%     %% Start Calculating Transition Rates
%     sign_old=0; % Initial 'order parameter' for the orbit. +1 for stable orbit with theta_plus, -1 for stable orbit with theta_miunus
%     sign_current=0;
%     num_transitions=0;
%     k_trans_omega=[];
%     for k=1:length(omega)
%         if omega(k)>omega_plus
%             sign_current=1;
% %             theta_stable(k)=theta_plus;
%         elseif omega(k)<omega_minus
%             sign_current=-1;
% %             theta_stable(k)=theta_minus;
%         end
%         % Updating sign: sign_old -> sign_current
%         if sign_current*sign_old==-1
%             num_transitions=num_transitions+1;
%             k_trans_omega=[k_trans_omega k];
%         end
%         
%         if sign_current==1
%             omega_stable(k)=omega_plus;
%         elseif sign_current==-1
%             omega_stable(k)=omega_minus;
%         end
%         
%         sign_old=sign_current;
%     end
%     
% end
% end
%% Calculating diffusion coefficients
% if 1==0
%     %% Calculating sigma_x for particle 1
% %     figure(1)
%     h_phi=histogram(diff(x(1,:)),100);
%     Values_phi=h_phi.Values;
%     Bins_phi=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
%     % plot(x,y);
%     [fitresult,gof]=fit_Gaussian(Bins_phi,Values_phi,plot_hist_fit_T_eff);
%     % a1*exp(-((x-b1)/c1)^2)
%     mu_x=fitresult.b1;
%     sigma_x=fitresult.c1/sqrt(2);
%     D_x=sigma_x^2/(2*dt);
%     %% Calculating sigma_y for particle 1
% %     figure(2)
%     h_phi=histogram(diff(y(1,:)),100);
%     Values_phi=h_phi.Values;
%     Bins_phi=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
%     % plot(x,y);
%     [fitresult,gof]=fit_Gaussian(Bins_phi,Values_phi,plot_hist_fit_T_eff);
%     % a1*exp(-((x-b1)/c1)^2)
%     mu_y=fitresult.b1;
%     sigma_y=fitresult.c1/sqrt(2);
%     D_y=sigma_y^2/(2*dt);
%     %% Calculating sigma_xy (I guess this doesn't work for circular orbits) for particle 1
%     h_phi=histogram(diff(x(1,:)).^2+diff(y(1,:)).^2)
%     perp=h_phi.Values;
%     hori=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
%     % plot(x,y);
%     [fitresult,gof]=fit_Gaussian(hori,perp,plot_hist_fit_T_eff);
%     % a1*exp(-((x-b1)/c1)^2)
%     mu_xy=fitresult.b1;
%     sigma_xy=fitresult.c1/sqrt(2);
%     D_xy=sigma_xy^2/(2*dt);
% 
%     
%     
% 
%     %% Calculating sigma_x for particle 2
% %     figure(1)
%     h_phi=histogram(diff(x(2,:)),100);
%     Values_phi=h_phi.Values;
%     Bins_phi=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
%     % plot(x,y);
%     [fitresult,gof]=fit_Gaussian(Bins_phi,Values_phi,plot_hist_fit_T_eff);
%     % a1*exp(-((x-b1)/c1)^2)
%     mu_x=fitresult.b1;
%     sigma_x=fitresult.c1/sqrt(2);
%     D_x=sigma_x^2/(2*dt);
%     %% Calculating sigma_y for particle 2
% %     figure(2)
%     h_phi=histogram(diff(y(2,:)),100);
%     Values_phi=h_phi.Values;
%     Bins_phi=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
%     % plot(x,y);
%     [fitresult,gof]=fit_Gaussian(Bins_phi,Values_phi,plot_hist_fit_T_eff);
%     % a1*exp(-((x-b1)/c1)^2)
%     mu_y=fitresult.b1;
%     sigma_y=fitresult.c1/sqrt(2);
%     D_y=sigma_y^2/(2*dt);
%     %% Calculating sigma_xy (I guess this doesn't work for circular orbits) for particle 2
%     h_phi=histogram(diff(x(2,:)).^2+diff(y(2,:)).^2)
%     perp=h_phi.Values;
%     hori=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
%     % plot(x,y);
%     [fitresult,gof]=fit_Gaussian(hori,perp,plot_hist_fit_T_eff);
%     % a1*exp(-((x-b1)/c1)^2)
%     mu_xy=fitresult.b1;
%     sigma_xy=fitresult.c1/sqrt(2);
%     D_xy=sigma_xy^2/(2*dt);
% end
