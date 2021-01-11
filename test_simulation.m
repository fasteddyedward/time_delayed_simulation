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
Date='2021.1.8' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
nth_take=300 % for a=5, particle 1 fixed
% nth_take=400 % for a=5, particle 1 not fixed
% nth_take=500 % for a=0, particle 1 not fixed
delta_t_matrix=2
T_matrix=[1]
v_0_matrix=[0.5:0.5:20]
dt=10^-1
intrinsic_delay=0.0 % Intrinsic delay
Obs_time_steps=10^5
% 
% Date='2021.1.8' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=400 % for a=5, particle 1 fixed
% delta_t_matrix=[1.8:0.1:4.0]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5

%%
D_eff_ratio_matrix_approx=[];
D_eff_ratio_matrix_full=[];
theta_0_matrix=[];
for delta_t_index=1:length(delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(v_0_matrix)
%% Input File Name
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
[movie_name,'.mat'];
load([movie_name,'.mat'])

%% Calculating angle
theta=[];
R=[];

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

%% Calculating D_eff
close all
recalculate_T_eff='yes'
    plot_hist_fit_T_eff='yes'
notation='theta'
switch recalculate_T_eff
    case 'yes'
        %% Finding the sigma of the omega (effective temperature)
        f=figure(80);clf;
        switch plot_hist_fit_T_eff
            case 'no'
                f.Visible='off';
        end
        switch notation
%             case 'omega' % This is for measuring sigma with histogram(diff(omega))
%                 h=histogram(diff(theta(2,:))/delta_t);  
            case 'theta' % This is for measuring sigma with histogram(diff(theta))
                
                h_theta=histogram(diff(theta));
                Values_theta=h_theta.Values;
                Bins_theta=h_theta.BinEdges(1:end-1)+0.5*h_theta.BinWidth;
                theta_0=v_0*delta_t/mean(R);
                if theta_0<=1
                    h_phi=histogram(abs(theta-theta_0*sin(theta)));
                else
                h_phi=histogram(abs(theta)-theta_0*sin(abs(theta)));
                end

        end
        Values_phi=h_phi.Values;
        Bins_phi=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
        Values_phi(Bins_phi<0.08)=[];
        Bins_phi(Bins_phi<0.08)=[];
        Bins_phi=[-flip(Bins_phi) Bins_phi];
        Values_phi=[flip(Values_phi) Values_phi];


        
        % plot(x,y);
        [fitresult,gof]=fit_Gaussian(Bins_phi,Values_phi,plot_hist_fit_T_eff);
        % a1*exp(-((x-b1)/c1)^2)
        mu_phi=fitresult.b1;
        sigma_phi=fitresult.c1/sqrt(2);
        
        [fitresult,gof]=fit_Gaussian(Bins_theta,Values_theta,plot_hist_fit_T_eff);
        mu_theta=fitresult.b1;
        sigma_theta=fitresult.c1/sqrt(2);

        %% Calculating the D_eff and T_eff
        D_phi=sigma_phi^2/(2);
        T_phi=D_phi;
        
        D_theta=sigma_theta^2/(2*dt);
        T_theta=D_theta;

    case 'no'
%         load([movie_name,'.mat'],'D_eff','T_eff','sigma','mu')
end


% theta_0=v_0*delta_t/mean(R);
theta_0_matrix=[theta_0_matrix theta_0];
D_eff_ratio_matrix_full=[D_eff_ratio_matrix_full D_phi/(D*delta_t^2/(mean(R)^2))];
D_eff_ratio_matrix_approx=[D_eff_ratio_matrix_approx D_theta/(4*D/(theta_0^2*mean(R)^2))];
%             D_phi
            nth_take=nth_take+1;
        end
        nth_take=nth_take+1;
    end
    nth_take=nth_take+1;
end
%%

close all
save(['nth_take=',num2str(nth_take),'.mat']);
nth_take
%% Calculating diffusion coefficients
if 1==0
    %% Calculating sigma_x for particle 1
%     figure(1)
    h_phi=histogram(diff(x(1,:)),100);
    Values_phi=h_phi.Values;
    Bins_phi=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(Bins_phi,Values_phi,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_x=fitresult.b1;
    sigma_x=fitresult.c1/sqrt(2);
    D_x=sigma_x^2/(2*dt);
    %% Calculating sigma_y for particle 1
%     figure(2)
    h_phi=histogram(diff(y(1,:)),100);
    Values_phi=h_phi.Values;
    Bins_phi=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(Bins_phi,Values_phi,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_y=fitresult.b1;
    sigma_y=fitresult.c1/sqrt(2);
    D_y=sigma_y^2/(2*dt);
    %% Calculating sigma_xy (I guess this doesn't work for circular orbits) for particle 1
    h_phi=histogram(diff(x(1,:)).^2+diff(y(1,:)).^2)
    perp=h_phi.Values;
    hori=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(hori,perp,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_xy=fitresult.b1;
    sigma_xy=fitresult.c1/sqrt(2);
    D_xy=sigma_xy^2/(2*dt);

    
    

    %% Calculating sigma_x for particle 2
%     figure(1)
    h_phi=histogram(diff(x(2,:)),100);
    Values_phi=h_phi.Values;
    Bins_phi=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(Bins_phi,Values_phi,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_x=fitresult.b1;
    sigma_x=fitresult.c1/sqrt(2);
    D_x=sigma_x^2/(2*dt);
    %% Calculating sigma_y for particle 2
%     figure(2)
    h_phi=histogram(diff(y(2,:)),100);
    Values_phi=h_phi.Values;
    Bins_phi=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(Bins_phi,Values_phi,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_y=fitresult.b1;
    sigma_y=fitresult.c1/sqrt(2);
    D_y=sigma_y^2/(2*dt);
    %% Calculating sigma_xy (I guess this doesn't work for circular orbits) for particle 2
    h_phi=histogram(diff(x(2,:)).^2+diff(y(2,:)).^2)
    perp=h_phi.Values;
    hori=h_phi.BinEdges(1:end-1)+0.5*h_phi.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(hori,perp,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_xy=fitresult.b1;
    sigma_xy=fitresult.c1/sqrt(2);
    D_xy=sigma_xy^2/(2*dt);
end
    
    
    
    
    
    
    