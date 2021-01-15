clear;
%% Parameters for importing mat files.
% Date='2020.12.10'
% nth_take=1
% delta_t_matrix=[0.5:0.5:16]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-1; % ms 
% % intrinsic_delay=0.0 % Intrinsic delay

% Date='2020.12.11'
% nth_take=1
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[3.5:0.1:10]
% dt=10^0
% intrinsic_delay=0.0 % Intrinsic delay


% Date='2020.11.24'
% nth_take=7
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[3.5:0.1:10]
% dt=10^-2




% Date='2020.12.11' % T_eff has been recalculated
% nth_take=1
% delta_t_matrix=[0.5:0.5:16]
% T_matrix=[1]
% v_0_matrix=5
% dt=5*10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5



% Date='2020.12.11'% T_eff has been recalculated % Matches Pretty well
% nth_take=7
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[3.5:0.1:10]
% dt=10^-2
% intrinsic_delay=0;
% Obs_time_steps=10^6;

% 
% Date='2020.12.11' % T_eff recalucluated, fits pretty well
% nth_take=100
% delta_t_matrix=[1.8:0.1:4]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5

% Date='2020.12.11' % T_eff recalucluated, fits not well
% nth_take=200
% Delta_t_matrix=[18:1:40]
% % Delta_t_matrix=2.3
% T_matrix=[1]
% v_0_matrix=0.5
% dt=10^0
% intrinsic_delay=0.0 % Intrinsic delay
% % Obs_time_steps=10^6
% Obs_time_steps=10^5

% Date='2020.12.15' %T_eff recalucluated, fits not well
% nth_take=1
% delta_t_matrix=[18:1:40]
% % Delta_t_matrix=2.3
% % Delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=0.5
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6

%%
% Date='2020.12.15' % T_eff has been recalculated
% nth_take=100
% delta_t_matrix=[1.8:0.1:4.0]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6
%%

% Date='2020.12.15' % T_eff has been recalculated 
% % The data look like a mess don't use
% nth_take=200
% delta_t_matrix=[1.8:0.1:4.0]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^0
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6

% Date='2020.12.15'% T_eff has been recalculated
% nth_take=300
% delta_t_matrix=[1.8:0.1:4.0]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^7

% Date='2020.12.16'
% nth_take=1
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[3.5:0.1:7]
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5

% Date='2020.12.16'
% nth_take=100
% delta_t_matrix=2
% T_matrix=[1]
% V_0_matrix=[3.5:0.1:10]
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^6


%%  Not in Laptop
% Date='2020.12.17'
% nth_take=1
% delta_t_matrix=2
% T_matrix=[1]
% V_0_matrix=[3.5:0.1:7.1]
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^7
% 

% Date='2021.1.6' % Note that the transition rates will be much higher than theoretical values because this is before bifurcation point
% nth_take=1
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[0.5:0.5:7.0]
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5


% Date='2021.1.6'
% nth_take=100
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[0.5:0.5:14]
% dt=10^-1
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5


% if length(Delta_t_matrix)>1
% load([Date,', delta_t_matrix, Obs_time=',num2str(Obs_time_steps),'.mat'])
% elseif length(V_0_matrix)>1
% load([Date,', v_0_matrix, Obs_time=',num2str(Obs_time_steps),'.mat'])
% end

%% Execution Parameters
recalculate_theta='yes' % normally just set to no
recalculate_hist='no' %
recalculate_R='no' % normally just set to no; R is calculated in hist already.
draw_hist='no'
recalculate_T_eff='yes'
    plot_hist_fit_T_eff='no'
fit_Viktor_method='no'
recalculate_T_Boltz='no'
    omega_plus_type='full_sine' % 'approximate'
    plot_Boltz_fit='off'

notation='theta';
% notation='omega';
%%
if exist('V_0_matrix','var')==0
    V_0_matrix=v_0_matrix;
end
if exist('Delta_t_matrix','var')==0
    Delta_t_matrix=delta_t_matrix;
end
%% Start Analzing
num_transitions_matrix=[];
theta_plus_matrix=[];
theta_minus_matrix=[];
R_matrix=[];

sigma_matrix=[];
mu_matrix=[];
T_eff_matrix=[];
D_eff_matrix=[];

T_Boltz_matrix=[];

std_R_matrix=[];
std_dR_matrix=[];
for delta_t_index=1:length(Delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(V_0_matrix)
            %             if  nth_take~=7
            if 1
                %             if ismember(nth_take,nth_interest)
close all
%% Input File Name
% movie_name=['test3']
% movie_name=['2020.11.25,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(V_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];
% movie_name=[Date,',dt=10e-3 take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];

[movie_name,'.mat'];
load([movie_name,'.mat'])

%% Calculating theta
switch recalculate_theta
    case 'yes'
        Analyze_theta=tic;
        moving_avg=1;
        plot_rot='no';
        Theta_Analysis_Fixed_Center(movie_name,partition_movie,N,v_0,Obs_time_steps,partition_time_steps,delta_t,dt,moving_avg,plot_rot);
        time_analyze_theta=toc(Analyze_theta)
    case 'no'
        load([movie_name,'.mat'],'theta')
        nth_take
end
%% Delete later
load([movie_name,'.mat'],'theta')
histogram(theta(2,:))
%% Plotting histogram
switch recalculate_hist
    case 'yes'
        moving_avg=1 ;
        num_bins=100  ;
        bin_limit=2;
        f=figure(81),clf % Note that hist_analysis only works for one particle orbitting a fixed particle at the moment
        switch draw_hist
            case 'no'
                f.Visible='off';
        end
        [k_trans,theta_plus,theta_minus,num_transitions]=hist_analysis(movie_name,moving_avg,num_bins,bin_limit,Obs_time_steps,delta_t,dt);
        set(gca, 'YScale', 'linear')
        saveas(gcf,[movie_name,' (hist).png'])
        
    case 'no'
        load([movie_name,'.mat'],'num_transitions','theta_plus','theta_minus')
        nth_take
end
%% Calculating R_mean. If recalculate_hist is 'yes' then recalculate_R can be set to 'no'
switch recalculate_R
    case 'yes'
        %% Calculating the histograms and stuff
%         theta=0; % 這個會把下面要用的theta蓋掉
        Analyze_theta=tic;
        moving_avg=1 ;
        plot_rot='no';
        Theta_Analysis_Fixed_Center(movie_name,partition_movie,N,v_0,Obs_time_steps,partition_time_steps,delta_t,dt,moving_avg,plot_rot);
        time_analyze_theta=toc(Analyze_theta)
        load([movie_name,'.mat'],'R_mean','R')
    case 'no'
        load([movie_name,'.mat'],'R_mean','R')
end
std_R=std(R(2,:));
std_dR=std(diff(R(2,:)));

%% Calculating Effective Temperature by fitting the histograwm of omega^dot
if exist('T_eff','var')==0
    recalculate_T_eff='yes'
end
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
        y=h.Values;
        x=h.BinEdges(1:end-1)+0.5*h.BinWidth;
        % plot(x,y);
        [fitresult,gof]=fit_Gaussian(x,y,plot_hist_fit_T_eff);
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

%% Calculating Boltzmann Temperature by fitting the histogram of omega
% if exist('T_Boltz','var')==0
%     recalculate_T_Boltz='yes'
% end
switch recalculate_T_Boltz % Modified from Boltzmann_dist_fit.m
    case 'yes'
        omega_0=v_0/R_mean(2);
        load([movie_name,'.mat'],'theta')
        bin_limit=1;
        num_bins=100;
        bin_interval=2*bin_limit/num_bins;
        bin_loc=-bin_limit:bin_interval:bin_limit;
        omega=theta(2,:)/delta_t;
        figure(82)
        h=histogram(omega,bin_loc);
        h.Visible='off';
        omega_count=h.Values;
        omega_bin=h.BinEdges(1:end-1)+0.5*h.BinWidth;
        dx=h.BinWidth;
        %         xline(theta_plus/delta_t)
        %         xline(theta_minus/delta_t)
        %% Fitting U=k_B*T*ln(p)+ln(Z)
        p=omega_count/length(omega); % Probability of the statistics
        ln_p=log(p/dx);
        omega_plus=sqrt(6/(omega_0*delta_t^3)*(omega_0*delta_t-1));
%         omega_plus_exact=sqrt(10-sqrt(1/42/(omega_0*delta_t)^6+120/(omega_0*delta_t)-20))/delta_t;
        U=@(omega) omega_0*delta_t^3/24*(omega.^2-2*omega_plus.^2).*omega.^2;
        
        f=figure(83);f.Visible=plot_Boltz_fit;
        
        [fitresult, gof] =Find_Boltzmann_Temp(omega_bin,ln_p,omega_0,delta_t,k_B,omega_plus_type);
        T_Boltz=fitresult.T;
        Z_fit=exp(fitresult.lnZ);
        xlabel('\omega')
        ylabel('ln(p (\omega))')
        title(['\omega_0 \deltat=',num2str(omega_0*delta_t),', T_{Boltzmann}/T_{eff}=',num2str(T_Boltz/T_eff),', sum(p(\omega))/Z_{fit}=',num2str(sum(exp(-U(bin_loc)/(k_B*T_Boltz))*dx)/Z_fit)])
        % Z_fit is because this is fitted with the fitting function, and it is
        % possible that the Z doens't give a completely normalized p(omega)
        save([movie_name,'.mat'],'T_Boltz','-append')
    case 'no'
        load([movie_name,'.mat'],'T_Boltz')
end
        
        
%% Appending the matrices
num_transitions_matrix=[num_transitions_matrix num_transitions];
theta_plus_matrix=[theta_plus_matrix, theta_plus];
theta_minus_matrix=[theta_minus_matrix, theta_minus];
R_matrix=[R_matrix R_mean(2)];

D_eff_matrix=[D_eff_matrix D_eff];
T_eff_matrix=[T_eff_matrix T_eff];
sigma_matrix=[sigma_matrix sigma];
mu_matrix=[mu_matrix mu];

T_Boltz_matrix=[T_Boltz_matrix T_Boltz];

std_R_matrix=[std_R_matrix std_R];
std_dR_matrix=[std_dR_matrix std_dR];

nth_take

            end
            nth_take=nth_take+1;
        end
        nth_take=nth_take+1;
    end
    nth_take=nth_take+1;
end
%% Saving results
time_duration=time(end);
clear R theta time v_omega v_x v_y x y 
% save(['plotting transitions v_0_matrix=3.5 0.1 10.mat'])
% save(['plotting transitions Delta_t_matrix=0.5 0.5 16.mat'])



%% Plotting and Analyzing the Bifurcation Diagram and Transision Rates with parameters a & b

switch fit_Viktor_method
    case 'yes'
% v_0_matrix(1)=[];
    %% Delta_t
    if length(Delta_t_matrix)>1
        close all
        %% Bifurcation Diagram
            %% Bifurcation Diagram (Original)
            figure(1);clf
            hold on
            plot(Delta_t_matrix,theta_plus_matrix,'.')
            plot(Delta_t_matrix,theta_minus_matrix,'.')
            title(['Bifurcation Diagram, v_0= ',num2str(v_0),', T=',num2str(T)])
            xlabel('\delta t')
            ylabel('\theta (rad)')
            x=@(delta_t,v_0,a)(v_0*delta_t)./(2*a);
            theta_theory=@(x)sqrt(10-sqrt(1/42./x.^6+120./x-20));
            theta_pl=@(x)sqrt(6./x.*(x-1));
            x_line=0:Delta_t_matrix(end)/1000:Delta_t_matrix(end);
            plot(x_line,+theta_theory(x(x_line,v_0,a)),'k')
            plot(x_line,-theta_theory(x(x_line,v_0,a)),'k')
            plot(x_line,theta_pl(x(x_line,v_0,a)),'g')
            plot(x_line,-theta_pl(x(x_line,v_0,a)),'g')
            legend('$\theta_+$','$\theta_-$','$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','','$\theta=\sqrt{\frac{6}{\omega_0 \delta t}(\omega_0 \delta t-1)}$','interpreter','latex','Location','northwest')
            saveas(gcf,['Bifurcation Diagram, v_0= ',num2str(v_0),', T=',num2str(T),'.png'])
            %% Bifurcation Diagram (Normalized)
            figure(2) ;clf

            hold on
            plot(v_0*Delta_t_matrix./R_matrix,theta_plus_matrix,'.')
            plot(v_0*Delta_t_matrix./R_matrix,theta_minus_matrix,'.')
            title(['Bifurcation Diagram, v_0= ',num2str(v_0),', T=',num2str(T)])
            xlabel('\omega_0 \delta t=v_0*\delta t/R')
            ylabel('\theta (rad)')
            x=@(delta_t,v_0,a)(v_0*delta_t)./(2*a);
            theta_theory=@(x)sqrt(10-sqrt(1/42./x.^6+120./x-20));
            theta_pl=@(x)sqrt(6./x.*(x-1));
            plot(x(x_line,v_0,a),+theta_theory(x(x_line,v_0,a)),'k')
            plot(x(x_line,v_0,a),-theta_theory(x(x_line,v_0,a)),'k')
            plot(x(x_line,v_0,a),theta_pl(x(x_line,v_0,a)),'g')
            plot(x(x_line,v_0,a),-theta_pl(x(x_line,v_0,a)),'g')
            legend('$\theta_+$','$\theta_-$','$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','','$\theta=\sqrt{\frac{6}{\omega_0 \delta t}(\omega_0 \delta t-1)}$','interpreter','latex','Location','northwest')
            saveas(gcf,['Bifurcation Diagram, v_0= ',num2str(v_0),', T=',num2str(T),'.png'])
        %% Transition Rates (modified from v_0)
            % Plotting original Data
            interest_flag=Delta_t_matrix>3
            delta_t_interest=Delta_t_matrix(interest_flag);
            num_transitions_interest=num_transitions_matrix(interest_flag);
            
            figure(3);clf
            hold on
            plot(delta_t_interest,num_transitions_interest,'.')
            plot(Delta_t_matrix(logical(1-interest_flag)),num_transitions_matrix(logical(1-interest_flag)),'x')
            xlabel('\delta t')
            ylabel('Number of Transitions')

            % Fitting with MATLAB tool
            figure(4)
            [fitresult, gof] = fit_trans_delta_t(delta_t_interest, num_transitions_interest)
            xlabel('\delta t')
            ylabel('Number of Transitions')
            fitresult.b

            % Transition Rates with normalized variables
            figure(5);clf
            hold on
            R_interest=R_matrix(interest_flag);
            omega_0_interest=delta_t_interest./(R_interest);
            plot(v_0*omega_0_interest,num_transitions_interest,'.')
            plot(v_0*Delta_t_matrix(logical(1-interest_flag))./R_matrix(logical(1-interest_flag)),num_transitions_matrix(logical(1-interest_flag)),'x')
            xlabel('v_0 \delta t/R')
            ylabel('Number of Transitions')

            % Fitting with MATLAB tool with normalized variables
            figure(6);clf
            omega_delta_t_interest=omega_0_interest*v_0; % x variable
            num_transitions_interest; % y variable
            [fitresult_norm, gof_norm] = fit_trans_delta_t_norm(omega_delta_t_interest, num_transitions_interest)
            xlabel('v_0 \delta t/R')
            ylabel('Number of Transitions')

    end


    %% V_0
    if length(V_0_matrix)>1
        close all
        %% Bifurcation Diagram
            %% Original
            figure(1) ;clf

            hold on
            plot(V_0_matrix,theta_plus_matrix,'.')
            plot(V_0_matrix,theta_minus_matrix,'.')
            title(['Bifurcation Diagram, \delta t= ',num2str(delta_t),', T=',num2str(T)])
            xlabel('v_0')
            ylabel('\theta (rad)')
            x=@(delta_t,v_0,a)  (v_0*delta_t)./(2*a);
            theta_theory=@(x)sqrt(10-sqrt(1/42./x.^6+120./x-20));
            theta_pl=@(x)sqrt(6./x.*(x-1));
            x_line=0:V_0_matrix(end)/1000:V_0_matrix(end);
            plot(x_line,+theta_theory(x(delta_t,x_line,a)),'k')
            plot(x_line,-theta_theory(x(delta_t,x_line,a)),'k')
            plot(x_line,theta_pl(x(delta_t,x_line,a)),'g')
            plot(x_line,-theta_pl(x(delta_t,x_line,a)),'g')
            legend('$\theta_+$','$\theta_-$','$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','','$\theta=\sqrt{\frac{6}{\omega_0 \delta t}(\omega_0 \delta t-1)}$','interpreter','latex','Location','northwest')
         
            
            %% Normalized
            figure(2) ;clf
            
            hold on
            plot(V_0_matrix*delta_t./(R_matrix),theta_plus_matrix,'.')
            plot(V_0_matrix*delta_t./(R_matrix),theta_minus_matrix,'.')
            title(['Bifurcation Diagram, \delta t= ',num2str(delta_t),', T=',num2str(T)])
            xlabel('\omega_0 \delta t=v_0*\delta t/R')
            ylabel('\theta (rad)')
            x=@(delta_t,v_0,a)(v_0*delta_t)./(2*a);
            theta_theory=@(x)sqrt(10-sqrt(1/42./x.^6+120./x-20));
            theta_pl=@(x)sqrt(6./x.*(x-1))
            plot(x(delta_t,x_line,a),theta_theory(x(delta_t,x_line,a)),'k')
            plot(x(delta_t,x_line,a),-theta_theory(x(delta_t,x_line,a)),'k')
            plot(x(delta_t,x_line,a),theta_pl(x(delta_t,x_line,a)),'g')
            plot(x(delta_t,x_line,a),-theta_pl(x(delta_t,x_line,a)),'g')
            legend('$\theta_+$','$\theta_-$','$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','','$\theta=\sqrt{\frac{6}{\omega_0 \delta t}(\omega_0 \delta t-1)}$','interpreter','latex','Location','northwest')

        %% Transition Rates

        % Plotting original Data
        interest_flag=V_0_matrix>5.5
        v_0_interest=V_0_matrix(interest_flag);
        num_transitions_interest=num_transitions_matrix(interest_flag);
        
        figure(3);clf
        hold on
        plot(v_0_interest,num_transitions_interest,'.')
        plot(V_0_matrix(logical(1-interest_flag)),num_transitions_matrix(logical(1-interest_flag)),'x')
        xlabel('v_0')
        ylabel('Number of Transitions')
        
        % Fitting with MATLAB tool
        figure(4)
        [fitresult, gof] = fit_trans_v_0(v_0_interest, num_transitions_interest)
        xlabel('v_0')
        ylabel('Number of Transitions')
%         fitresult.b
        
        % Transition Rates with normalized variables
        figure(5);clf
        hold on
        R_interest=R_matrix(interest_flag);
        omega_0_interest=v_0_interest./(R_interest);
        plot(delta_t*omega_0_interest,num_transitions_interest,'.')
        plot(delta_t*V_0_matrix(logical(1-interest_flag))./R_matrix(logical(1-interest_flag)),num_transitions_matrix(logical(1-interest_flag)),'x')
        xlabel('v_0 \delta t/R')
        ylabel('Number of Transitions')
        
        % Fitting with MATLAB tool with normalized variables
        figure(6);clf
        omega_delta_t_interest=omega_0_interest*delta_t; % x variable
        num_transitions_interest; % y variable
        [fitresult_norm, gof_norm] = fit_trans_v_0_norm(omega_delta_t_interest, num_transitions_interest)
        xlabel('v_0 \delta t/R')
        ylabel('Number of Transitions')
%         fitresult_norm.b/50
        
        
        fitresult.b
        fitresult_norm.b/50
    end
end

%% Save Variables for Viktor
% save('For_Viktor.mat','a','D','delta_t','dt','k_B','num_transitions_matrix','T','R_matrix','v_0_matrix')
%% D_eff and T_eff

D_eff_matrix=sigma_matrix.^2/(2*dt);
T_eff_matrix=D_eff_matrix;
%% Transition Rates
omega_0_matrix=V_0_matrix./R_matrix;
switch notation
    case 'omega'
        
        transition_rate_theory=2*sqrt(2)./(pi*omega_0_matrix.*Delta_t_matrix.^2).*(omega_0_matrix.*Delta_t_matrix-1).*exp(-3*(omega_0_matrix.*Delta_t_matrix-1).^2./(D_eff_matrix.*omega_0_matrix.^2.*Delta_t_matrix.^5));
% transition_rate_theory=2*sqrt(2)./(pi.*omega_0_matrix.*Delta_t_matrix.^2).*(omega_0_matrix.*Delta_t_matrix-1).*exp(-3/2*(omega_0_matrix.*Delta_t_matrix-1).^2./(omega_0_matrix.*k_B.*T_eff_matrix.*Delta_t_matrix.^3));
    case 'theta'
        theta_0_matrix=omega_0_matrix.*Delta_t_matrix;
        transition_rate_theory=2*sqrt(2)./(pi*theta_0_matrix.*Delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(D_eff_matrix.*Delta_t_matrix.*theta_0_matrix.^2));
end
        
        % transition_rate_theory=sqrt(2)./(pi.*omega_0_matrix.*Delta_t_matrix.^2).*(omega_0_matrix.*Delta_t_matrix-1).*exp(-3/2*(omega_0_matrix.*Delta_t_matrix-1).^2./(D_eff_matrix.*omega_0_matrix.*Delta_t_matrix.^2/2.*omega_0_matrix.*Delta_t_matrix.^3));


%% 2020.12.10 Plotting Transition Rates 
% time_duration=Obs_time_steps*dt+delta_t;
x0=omega_0_matrix.*Delta_t_matrix;
figure(7);clf;hold on;
plot(x0,transition_rate_theory,'o')
plot(x0,num_transitions_matrix/time_duration,'x')
xlabel('\omega_0 \delta t')
ylabel('Transition Rates (1/s)')


figure(8);clf;hold on;
flag=(x0>1.);
plot(x0(flag),transition_rate_theory(flag),'o')
plot(x0(flag),num_transitions_matrix(flag)/time_duration,'x')
xlabel('\omega_0 \delta t')
ylabel('Transition Rates (1/s)')

%% Plotting inverse transition rate
figure(11);clf;hold on;
plot(x0,1./transition_rate_theory,'o')
plot(x0,1./(num_transitions_matrix/time_duration),'x')
xlabel('\omega_0 \delta t')
ylabel('1/Transition Rates (s)')


figure(12);clf;hold on;
flag=(x0>1 );
plot(x0(flag),1./transition_rate_theory(flag),'o')
plot(x0(flag),1./(num_transitions_matrix(flag)/time_duration),'x')
xlabel('\omega_0 \delta t')
ylabel('1/Transition Rates (s)')

%% Comparing between 1-D and 2-D % It's the same, so just for check actually
% figure(13);clf;hold on;
% % title('R vs 2a')
% % omega_0_matrix=V_0_matrix./(R_matrix);
% gamma_eff_matrix=omega_0_matrix.*Delta_t_matrix.^2/2;
% T_eff_1dim=gamma_eff_matrix.*D_eff_matrix/k_B;
% ratio=T_eff_matrix./T_eff_1dim;
% plot(x0,ratio)
% xlabel('\omega_0 \delta t')


%% Comparing E_b and k_B*T
figure(14);clf;
plot(x0,-3/2*(omega_0_matrix.*Delta_t_matrix-1).^2./(omega_0_matrix.*k_B.*T_eff_matrix.*Delta_t_matrix.^3))
xlabel('\omega_0 \delta_t')
ylabel('E_b/k_B T')
%% Plotting R vs omega_0 delta_t
figure(10);clf;hold on;
errorbar(Delta_t_matrix*V_0_matrix./R_matrix , R_matrix/(2*a), std_dR_matrix/(2*a));
yline(1)
xlabel('\omega_0 \delta t ')
ylabel('R/(2a)')
axis([-inf inf 0 inf])
legend('(std of R)/2a','R=2a','Location','northwest')
%% Plotting std_R and std_dR
figure(22); clf; hold on;
title('Fluctuations on the Radial Direction')
plot(x0,std_R_matrix)
plot(x0,std_dR_matrix)
xlabel('omega_0 \delta t')
legend('r', 'dr')

%% Plotting T_eff and D_eff
figure(15);clf;
plot(x0,T_eff_matrix)
xlabel('\omega_0 \delta_t')
ylabel('T_{eff}')
figure(16);clf;
plot(x0,D_eff_matrix)
xlabel('\omega_0 \delta_t')
ylabel('D_{eff}')

%% Plotting T_eff and D_eff & R vs omega_0 * delta_t
figure(21);clf;hold on
% title('fixed \delta t')
% title('fixed \omega_0')
plot(x0,T_eff_matrix)
plot(x0,D_eff_matrix)
plot(x0,R_matrix/(2*a)/300)
xlabel('\omega_0 \delta t')
legend('T_{eff}','D_{eff}','R/(2a)/300','Location','northwest')

%% Plotting D_eff v.s D/R^2
figure(23);clf;hold on
title('D_{eff} v.s. D_0')
plot(x0,D_eff_matrix./(D./R_matrix.^2))
legend('D_{eff}/(D/R^2)','Location','northwest')
xlabel('\omega_0 \delta t')
axis([-inf inf 0 inf])
%% Storing Mat File
if length(Delta_t_matrix)>1
% save(['Transitions, ',Date,', delta_t_matrix, Obs_time=',num2str(Obs_time_steps),'.mat'])
elseif length(V_0_matrix)>1
% save(['Transitions, ',Date,', v_0_matrix, Obs_time=',num2str(Obs_time_steps),'.mat'])
end

