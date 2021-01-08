%% This file calculates the D_x D_y D_eff for the experimental data, and see if D_eff=D_x/R^2
%% But still there is no pure diffusion data from the experiments, so we don't know.
clear;
close all
load('dt=1.500.tdms_angle.mat')
%%
dt=0.03 % 30 ms as Xiangzun said
notation='theta'
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
        switch notation
            case 'omega'
                h=histogram(diff(theta(2,:))/delta_t);  % This is for measuring sigma with histogram(diff(omega))
            case 'theta'
                h=histogram(diff(theta(2,:))); % This is for measuring sigma with histogram(diff(theta))
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
        save([movie_name,'.mat'],'D_eff','T_eff','sigma','mu','-append')
    case 'no'
        load([movie_name,'.mat'],'D_eff','T_eff','sigma','mu')
end

%% Calculating sigma_x and sigma_y
    %% Calculating sigma_x
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
    %% Calculating sigma_y
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
    %% Calculating sigma_xy (I guess this doesn't work for circular orbits)
%     h=histogram(diff(x(2,:)).^2+diff(y(2,:)).^2)
%     perp=h.Values;
%     hori=h.BinEdges(1:end-1)+0.5*h.BinWidth;
%     % plot(x,y);
%     [fitresult,gof]=fit_Gaussian(hori,perp,plot_hist_fit_T_eff);
%     % a1*exp(-((x-b1)/c1)^2)
%     mu_xy=fitresult.b1;
%     sigma_xy=fitresult.c1/sqrt(2);
%     D_xy=sigma_xy^2/(2*dt);
%     
    
    %% Plotting D_eff v.s D/R^2
%     D=(D_x+D_y)/2;
    D=0.065; % Measured from Free_Diffusion.m from Xiangzun's old diffusive particle
    v_0 =2.1; % guessed
%     delta_t=1.5
%     theta_0=v_0/R_mean(2)*delta_t;
    D_eff;
    D/R_mean(2)^2;
    'D_eff/(D/R^2)'
    D_eff/(D/R_mean(2)^2)
%     D_eff/(4*D/(theta_0^2*R_mean(2))^2)
%     axis([-inf inf 0 inf])
    %% hist R
    figure(3)
    
    if exist('R','var')==0
        R(1:2,1:length(x))=0;
        R(2,:)=sqrt((x(2,:)-x(1,:)).^2+(y(2,:)-y(1,:)).^2);
        R(1,:)=0;
    end
    hist(R(2,:),1000)
    axis([0 inf 0 inf])