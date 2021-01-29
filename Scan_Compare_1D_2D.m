% 2020.1.27 This file is same as Compare_1D_2D_pc_Vik.m, but modified to 
%% scan over omega_0 and delta_t.
clear 
% close all
% cd('dt=0.1')

%% Running the compare functions
dt=0.01
D_0=10
delta_t_matrix_vik=[0.1 0.2 0.5 1 2 5 10]'  %% Make 1-D simulation
delta_t_matrix=[0.1 0.2 0.5 1 2 5 10]' %% for 2-D simulation
theta_0_1D_matrix=(logspace(log10(1),log10(3),30));

%%
file_name_pc=['2021.1.28_scan_pc_D_0=',num2str(D_0)];
file_name_vik=['2021.1.28_scan_vik_D_0=',num2str(D_0)];
%% Running the subfunctions
% compare_Viktor(D_0,dt,tV_max,theta_0_matrix,tau_matrix,runs)
% compare_Viktor(D_0,dt,500,linspace(1.1,1.6,15),tau_matrix,1)
% compare_Viktor(D_0,dt,50000,linspace(1.1,1.6,15),[1.8:0.1:4],1)

Scan_compare_Viktor(D_0,dt,50000,linspace(1.1,1.6,15),delta_t_matrix_vik,1,file_name_vik)


% compare_PC(D_0,dt,delta_t_matrix,v_0_matrix,Obs_time_steps)
% compare_PC(D_0,dt,tau_matrix,[11:0.5:20],10^4)
% compare_PC(D_0,dt,[2:0.2:4],5,10^6)

% Scan_compare_PC(D_0,dt,delta_t_matrix,theta_0_1D_matrix,10^6,file_name_pc)
Scan_compare_PC(D_0,dt,delta_t_matrix,theta_0_1D_matrix,10^5,file_name_pc)


%% Part 1, from modified_average_rate... from Viktor
    load([file_name_vik,'.mat'])
    figure(4);clf;hold on;

    % theta_0_1D_Vik=linspace(1.1,1.6,15)';
    for i=1:length(delta_t_matrix_vik)
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:),'k-x')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(1*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'r--')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(2*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'r--')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(1*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(2*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(1*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(2*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
    end
    v = [5 2 5];
    [caz,cel] = view(v);
    grid on
    xlabel('\delta t')
    ylabel('\theta_0')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'Zscale','log')
    axis([-inf inf 1.0 inf 1/500000 0.1])
    
    %% plot for only one plot
    figure(5);clf;hold on;
    tau=5/3
    
   
    title(['D_0=',num2str(D_0),', \tau=',num2str(tau)])
    for i=1:length(delta_t_matrix_vik)
        if delta_t_matrix_vik(i)==tau
        plot(AvN_plot,1./AratesV_plot(i,:),'k-x')
        plot(AvN_plot,(1*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'r--')
        plot(AvN_plot,(2*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'r--')
        plot(AvN_plot,(1*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
        plot(AvN_plot,(2*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
        plot(AvN_plot,(1*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
        plot(AvN_plot,(2*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
        end
    end

    set(gca,'YScale','log')
1
%%

%% Part 2: Pin-Chuan's 2-d simulation
    clearvars -except D_0 file_name_pc file_name_vik
    % load(['2021.1.25_compare_PC,D_0=',num2str(D_0),'.mat'])
    load([file_name_pc,'.mat'])
    %% Calculating the numeric solution for FPE
    iF = @(x) (sign(x)+1)./(sign(x)+1);
    figure(100);clf;hold on;
    [rateN_2D,rateA_2D,rateAeff_2D]=numeric_FPE(delta_t_matrix,theta_0_matrix,D_theta_matrix)
    for i=1:length(delta_t_matrix)
        trans_rate_sim_2D(i,:)=num_transitions_matrix(i,:)/time_duration(i);
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:),'xk-')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*rateN_2D(i,:),'-.g')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*(2*rateN_2D(i,:)),'-.g')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*rateA_2D(i,:),'-.r')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*(2*rateA_2D(i,:)),'-.r')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*rateAeff_2D(i,:),':b')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*(2*rateAeff_2D(i,:)),':b')
        
    end
        v = [5 2 5];
    [caz,cel] = view(v);
    grid on
    xlabel('\delta t')
    ylabel('\theta_0')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'Zscale','log')
    axis([-inf inf 1.0 inf -inf inf])
    %%  plotting theta_0_matrix
    figure(1);clf;hold on
    for i=1:length(delta_t_matrix)
        plot(theta_0_1D_matrix,theta_0_matrix(i,:),'x')
    end

    %% Plotting for the 2-D simulations in 3D
    % Calculating the Kramer's Rate
    % close all
    trans_rate_sim_2D(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
    trans_rate_measured(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
    transition_rate_Kramer(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;
    transition_rate_Kramer_eff(1:length(delta_t_matrix),1:length(theta_0_1D_matrix))=0;

    for i=1:length(delta_t_matrix)
        trans_rate_sim_2D(i,:)=num_transitions_matrix(i,:)/time_duration(i);
        trans_rate_measured(i,:)=1*sqrt(2)./(pi*theta_0_matrix(i,:).*delta_t_matrix(i)).*(theta_0_matrix(i,:)-1).*exp(-3*(theta_0_matrix(i,:)-1).^2./(D_theta_matrix(i,:).*theta_0_matrix(i,:).^2.*delta_t_matrix(i)));
        transition_rate_Kramer(i,:)=1*sqrt(2)./(pi*theta_0_matrix(i,:).*delta_t_matrix(i)).*(theta_0_matrix(i,:)-1).*exp(-3*(theta_0_matrix(i,:)-1).^2./(4*D_0./R_matrix(i,:).^2.*delta_t_matrix(i)));
        transition_rate_Kramer_eff(i,:)=1*sqrt(2)./(pi*theta_0_matrix(i,:).*delta_t_matrix(i)).*(theta_0_matrix(i,:)-1).*exp(-3*(theta_0_matrix(i,:)-1).^2./(2*D_0./R_matrix(i,:).^2.*delta_t_matrix(i).*theta_0_matrix(i,:).^2));
    end
    figure(2);clf;hold on;
    for i=1:length(delta_t_matrix)
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:),'xk-')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_measured(i,:),'-.g')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),2*trans_rate_measured(i,:),'-.g')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),transition_rate_Kramer(i,:),'-.r')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),2*transition_rate_Kramer(i,:),'-.r')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),transition_rate_Kramer_eff(i,:),':b')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),2*transition_rate_Kramer_eff(i,:),':b')
    end
    v = [5 2 5];
    [caz,cel] = view(v);
    grid on
    xlabel('\delta t')
    ylabel('\theta_0')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'Zscale','log')
    axis([-inf inf 1.0 inf -inf inf]) %% The minimum non-zero transition rate is 1/Obs_time_steps.
    %% plotting for bifurcation
    figure(3);clf;hold on;
    for i=1:length(delta_t_matrix)
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),theta_plus_matrix(i,:),'xb-')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),theta_minus_matrix(i,:),'xr-')
    end
    v = [5 2 5];
    [caz,cel] = view(v);
    grid on
    xlabel('\delta t')
    ylabel('\theta_0')
    set(gca,'XScale','log')
    % set(gca,'YScale','log')
    % set(gca,'Zscale','log')
    %% Plotting for 2D simulations for one specific data set
    tau=10
    i=find(delta_t_matrix==tau)
    figure(3);clf;hold on;
    plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:),'xk')
    plot(theta_0_matrix(i,:),trans_rate_measured(i,:),'-.g')
    plot(theta_0_matrix(i,:),2*trans_rate_measured(i,:),'-.g')
    plot(theta_0_matrix(i,:),transition_rate_Kramer(i,:),'-.r')
    plot(theta_0_matrix(i,:),2*transition_rate_Kramer(i,:),'-.r')
    plot(theta_0_matrix(i,:),transition_rate_Kramer_eff(i,:),':b')
    plot(theta_0_matrix(i,:),2*transition_rate_Kramer_eff(i,:),':b')
    set(gca,'YScale','log')
    set(gca,'Xscale','log')

    
    
    
    
    
    
    
    
    
    
    
    
    
%     
%     
%     
%     
%     
% %%
% 
% % load([file_name_vik,'.mat'])
% load(['2021.1.25_compare_vik,D_0=',num2str(D_0),'.mat'],'AvN_plot','AratesV_plot','rateA_plot','rateAeff_plot','rateN_plot','iF')
% % Plotting figures
% figure(2);clf;title(['1-D, D_0=',num2str(D_0)])
% 
% 
% hold on
% plot(AvN_plot, AratesV_plot,'x','DisplayName','simulation')
% plot(AvN_plot,1./(1*rateN_plot).*iF((-1 + AvN_plot)),'r','LineStyle','--','DisplayName','numerics')
% plot(AvN_plot,1./(1*rateA_plot).*iF((-1 + AvN_plot)),'g','LineStyle','-.','DisplayName','Kramers')
% plot(AvN_plot,1./(1*rateAeff_plot).*iF((-1 + AvN_plot)),'b','LineStyle',':','DisplayName','KramersEff')
% plot(AvN_plot,1./(2*rateN_plot).*iF((-1 + AvN_plot)),'r','LineStyle','--','DisplayName','numerics2')
% plot(AvN_plot,1./(2*rateA_plot).*iF((-1 + AvN_plot)),'g','LineStyle','-.','DisplayName','Kramers2')
% plot(AvN_plot,1./(2*rateAeff_plot).*iF((-1 + AvN_plot)),'b','LineStyle',':','DisplayName','KramersEff2')
% xlabel('$\omega\delta t$','Interpreter','latex');
% ylabel('$1/k$','Interpreter','latex');
% 
% 
% 
% figure(1);clf;title(['1-D, D_0=',num2str(D_0)])
% hold on
% plot(AvN_plot, 1./AratesV_plot,'x','DisplayName','simulation')
% plot(AvN_plot,(1*rateN_plot).*iF((-1 + AvN_plot)),'r','LineStyle','--','DisplayName','numerics')
% plot(AvN_plot,(1*rateA_plot).*iF((-1 + AvN_plot)),'g','LineStyle','-.','DisplayName','Kramers')
% plot(AvN_plot,(1*rateAeff_plot).*iF((-1 + AvN_plot)),'b','LineStyle',':','DisplayName','KramersEff')
% plot(AvN_plot,(2*rateN_plot).*iF((-1 + AvN_plot)),'r','LineStyle','--','DisplayName','numerics2')
% plot(AvN_plot,(2*rateA_plot).*iF((-1 + AvN_plot)),'g','LineStyle','-.','DisplayName','Kramers2')
% plot(AvN_plot,(2*rateAeff_plot).*iF((-1 + AvN_plot)),'b','LineStyle',':','DisplayName','KramersEff2')
% xlabel('$\omega\delta t$','Interpreter','latex');
% ylabel('$k$','Interpreter','latex');
% 
% 
% %% Part 1.2, plotting Kramer_eff and Kramer for 2-D, modified for Overleaf
% figure(1);
% set(gca,'YScale','log')
% set(gca,'Xscale','log')
% % legend('1-D simulation','2-D k_{D_{eff}}','2-D k{_D_{eff}} * 2','2-D k_{4D_0/(R^2\theta_0^2)}','2-D k_{4D_0/(R^2\theta_0^2)} * 2','2-D k_{2D_0/R^2}','2-D k_{2D_0/R^2} * 2','Location','southwest')
% % legend('simulation','numerics','Kramers','KramersEff','numerics*2','Kramers*2','KramersEff*2','2-D simulation','2-D measured','2-D measured*2','Location','northwest')
% legend('1-D simulation','1-D numerics', '1-D k_{4D_0/(R^2\theta_0^2)}','1-D k_{2D_0/(R^2)}','1-D numerics * 2','1-D k_{4D_0/(R^2\theta_0^2)} * 2','1-D k_{2D_0/(R^2)} * 2','Location','southwest')
% 
% figure(2)
% set(gca,'YScale','log')
% set(gca,'Xscale','log')
% % legend('2-D simulation','2-D k_{D_{eff}}','2-D k{_D_{eff}} * 2','2-D k_{4D_0/(R^2\theta_0^2)}','2-D k_{4D_0/(R^2\theta_0^2)} * 2','2-D k_{2D_0/R^2}','2-D k_{2D_0/R^2} * 2','Location','northwest')
% legend('1-D simulation','1-D numerics', '1-D k_{4D_0/(R^2\theta_0^2)}','1-D k_{2D_0/(R^2)}','1-D numerics * 2','1-D k_{4D_0/(R^2\theta_0^2)} * 2','1-D k_{2D_0/(R^2)} * 2','Location','northwest')
% 
% 
% 
% 
% %% Part 3, from plotting with D_eff measured with Gaussian
% trans_rate_measured=1*sqrt(2)./(pi*theta_0_matrix.*delta_t_matrix).*(theta_0_matrix-1).*exp(-3*(theta_0_matrix-1).^2./(D_theta_matrix.*theta_0_matrix.^2.*delta_t_matrix));
% figure(2)
% plot(theta_0_matrix,1./trans_rate_measured,'m--')
% plot(theta_0_matrix,1./(2*trans_rate_measured),'m--')
% figure(1)
% plot(theta_0_matrix,trans_rate_measured,'m--')
% plot(theta_0_matrix,2*trans_rate_measured,'m--')
% 
% 
