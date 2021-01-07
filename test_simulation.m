%% This file tests if the pure diffusion case gives the correct diffusion constant
%% For 2 particles purely diffusing with D=1, the diffusive constant measured is correct and gave 1.
clear
close all

Date='2021.1.7' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=1 % for a=5, particle 1 fixed
% nth_take=100 % for a=5, particle 1 not fixed
nth_take=200 % for a=0, particle 1 not fixed
delta_t_matrix=2
T_matrix=[1]
v_0_matrix=[4:0.5:10]
dt=10^-1
intrinsic_delay=0.0 % Intrinsic delay
Obs_time_steps=10^5

D_eff_ratio_matrix=[];
for delta_t_index=1:length(delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(v_0_matrix)
%% Input File Name
% movie_name=['test3']
% movie_name=['2020.11.25,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% movie_name=[Date,',dt=10e-3 take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];

[movie_name,'.mat'];
load([movie_name,'.mat'])

%% Calculating angle
theta=[];
R=[]
for k=1:Obs_time_steps
    %% Method one
%     R1=[x(2,k)-x(1,k),y(2,k)-y(1,k),0];
%     R2=[x(2,k+round(delta_t/dt))-x(1,k+round(delta_t/dt)),y(2,k+round(delta_t/dt))-y(1,k+round(delta_t/dt)),0];
%     R1_cross_R2=cross(R1,R2);
%     theta(k)=asin(R1_cross_R2(3)/(norm(R1)*norm(R2)));
%     
    
    
    %% Method two
            R1=[x(2,k)-x(1,k),y(2,k)-y(1,k),0];
            R2=[x(2,k+round(delta_t/dt))-x(1,k+round(delta_t/dt)),y(2,k+round(delta_t/dt))-y(1,k+round(delta_t/dt)),0];
            diff_x=[x(2,k+round(delta_t/dt))-x(2,k),y(2,k+round(delta_t/dt))-y(2,k),0];
            %% Determining Sign of theta
            R1_cross_diff_x=cross(R1,diff_x);
%             theta(k)=sign(R1_cross_diff_x(3))*acos(((norm(R1)^2+norm(R2)^2-norm(diff_x)^2)/(2*norm(R1)*norm(R2))));
            theta(k)=sign(R1_cross_diff_x(3))*acos(dot(R1,R2)/(norm(R1)*norm(R2)));
 %%   
    
    
    if fixed_flag(1)==1
        R(k)=norm(R1);
    elseif fixed_flag(1)==0
        R(k)=norm(R1)/2;
    end
end

plot(theta)
pause
%% Calculating D_eff
recalculate_T_eff='yes'
    plot_hist_fit_T_eff='yes'
notation='theta'
% D_eff=0
switch recalculate_T_eff
    case 'yes'
        %% Finding the sigma of the omega (effective temperature)
        f=figure(80);clf;
        switch plot_hist_fit_T_eff
            case 'no'
                f.Visible='off';
        end
        switch notation
%             case 'omega'
%                 h=histogram(diff(theta(2,:))/delta_t);  % This is for measuring sigma with histogram(diff(omega))
            case 'theta'
                h=histogram(diff(theta)); % This is for measuring sigma with histogram(diff(theta))
        end
        Values=h.Values;
        Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
        % plot(x,y);
        [fitresult,gof]=fit_Gaussian(Bins,Values,plot_hist_fit_T_eff);
        % a1*exp(-((x-b1)/c1)^2)
        mu=fitresult.b1;
        sigma=fitresult.c1/sqrt(2);
        
        %% Calculating the D_eff and T_eff
        D_eff=sigma^2/(2*dt);
        T_eff=D_eff;
%                 T_eff=(v_0/R_mean(2))*delta_t^2*sigma^2/(4*k_B*dt); % This is for measuring sigma with histogram(diff(theta)/delta_t)=histogram(diff(omega))
%         D_eff=sigma^2*(v_0/R_mean(2)*delta_t)^2/(8*k_B*dt); % This is for measuring sigma with histogram(diff(theta))
%         T_eff=0;
%         save([movie_name,'.mat'],'D_eff','T_eff','sigma','mu','-append')
    case 'no'
%         load([movie_name,'.mat'],'D_eff','T_eff','sigma','mu')
end
% D_eff
% D_eff/(D/mean(R)^2)
D_eff_ratio_matrix=[D_eff_ratio_matrix D_eff/(D/mean(R)^2)];
            
            nth_take=nth_take+1;
        end
        nth_take=nth_take+1;
    end
    nth_take=nth_take+1;
end
%%
close all
save(['nth_take=',num2str(nth_take),'.mat']);

%% Calculating diffusion coefficients
if 1==0
    %% Calculating sigma_x for particle 1
%     figure(1)
    h=histogram(diff(x(1,:)),100);
    Values=h.Values;
    Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(Bins,Values,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_x=fitresult.b1;
    sigma_x=fitresult.c1/sqrt(2);
    D_x=sigma_x^2/(2*dt);
    %% Calculating sigma_y for particle 1
%     figure(2)
    h=histogram(diff(y(1,:)),100);
    Values=h.Values;
    Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(Bins,Values,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_y=fitresult.b1;
    sigma_y=fitresult.c1/sqrt(2);
    D_y=sigma_y^2/(2*dt);
    %% Calculating sigma_xy (I guess this doesn't work for circular orbits) for particle 1
    h=histogram(diff(x(1,:)).^2+diff(y(1,:)).^2)
    perp=h.Values;
    hori=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(hori,perp,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_xy=fitresult.b1;
    sigma_xy=fitresult.c1/sqrt(2);
    D_xy=sigma_xy^2/(2*dt);

    
    

    %% Calculating sigma_x for particle 2
%     figure(1)
    h=histogram(diff(x(2,:)),100);
    Values=h.Values;
    Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(Bins,Values,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_x=fitresult.b1;
    sigma_x=fitresult.c1/sqrt(2);
    D_x=sigma_x^2/(2*dt);
    %% Calculating sigma_y for particle 2
%     figure(2)
    h=histogram(diff(y(2,:)),100);
    Values=h.Values;
    Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(Bins,Values,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_y=fitresult.b1;
    sigma_y=fitresult.c1/sqrt(2);
    D_y=sigma_y^2/(2*dt);
    %% Calculating sigma_xy (I guess this doesn't work for circular orbits) for particle 2
    h=histogram(diff(x(2,:)).^2+diff(y(2,:)).^2)
    perp=h.Values;
    hori=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(hori,perp,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_xy=fitresult.b1;
    sigma_xy=fitresult.c1/sqrt(2);
    D_xy=sigma_xy^2/(2*dt);
end
    
    
    
    
    
    
    