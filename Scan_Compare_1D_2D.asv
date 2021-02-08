% 2020.1.27 This file is same as Compare_1D_2D_pc_Vik.m, but modified to 
%% scan over omega_0 and delta_t.
clear 
close all
% cd('dt=0.1')

%% Running the compare functions
dt=0.01
D_0=1

% delta_t_matrix_vik=[0.1 0.2 0.5 1 2 5 10]'  %% Make 1-D simulation
% delta_t_matrix=[0.1 0.2 0.5 1 2 5 10]' %% for 2-D simulation
delta_t_matrix=2
% theta_0_1D_matrix=(logspace(log10(1),log10(3),30));

intrinsic_delay=1
theta_0_1D_matrix=(logspace(log10(1),log10(3),15));
file_name_pc=['2021.2.5_without_int_delay_D_0=',num2str(D_0)];


intrinsic_delay=0.01
% theta_0_1D_matrix=(logspace(log10(1),log10(3),15));
theta_0_1D_matrix=(linspace(0.1,3,15));
file_name_pc=['2021.2.8_with_int_delay_D_0=',num2str(D_0)];

% file_name_vik=['2021.1.28_scan_vik_D_0=',num2str(D_0)];


% D_0=20
% file_name_pc=['2021.2.4_scan_pc_D_0=',num2str(D_0)];
% file_name_vik=['2021.2.4_scan_vik_D_0=',num2str(D_0)];
%% Running the subfunctions
% compare_Viktor(D_0,dt,tV_max,theta_0_matrix,tau_matrix,runs)
% compare_Viktor(D_0,dt,500,linspace(1.1,1.6,15),tau_matrix,1)
% compare_Viktor(D_0,dt,50000,linspace(1.1,1.6,15),[1.8:0.1:4],1)

% Scan_compare_Viktor(D_0,dt,500,linspace(1.1,1.6,15),delta_t_matrix_vik,1,file_name_vik)


% compare_PC(D_0,dt,delta_t_matrix,v_0_matrix,Obs_time_steps)
% compare_PC(D_0,dt,tau_matrix,[11:0.5:20],10^4)
% compare_PC(D_0,dt,[2:0.2:4],5,10^6)

Scan_compare_PC(D_0,dt,delta_t_matrix,theta_0_1D_matrix,10^4,file_name_pc,intrinsic_delay)


%% Part 1: 1-dim simulations from Viktor
%% 
%     clearvars -except D_0 file_name_pc file_name_vik
    load([file_name_vik,'.mat'])
    %% 3D plot for all parameter space
    figure(1);clf;hold on;title(['1-D Simulation vs Kramer Escape Rate, D_0=',num2str(D_0)])

    % theta_0_1D_Vik=linspace(1.1,1.6,15)';
    for i=1:length(delta_t_matrix_vik)
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:),'k-x')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(1*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(2*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(1*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'r-.')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(2*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'r-.')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(1*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,(2*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
    end
    sim=plot(1,1,'k-x')
    numeric=plot(1,1,'g-.')
    Kramer=plot(1,1,'r-.')
    Kramer_eff=plot(1,1,'b:')
    legend([sim numeric Kramer Kramer_eff],{'Simulation', 'Numeric FPE with D_\theta=4D_0/(R^2\theta_0^2))', 'D_\theta=4D_0/(R^2\theta_0^2))','D_\theta=2D_0/R^2'})
    v = [5 2 5];
    [caz,cel] = view(v);
    grid on
    xlabel('\delta t')
    ylabel('\theta_0')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'Zscale','log')
    axis([-inf inf 1.0 inf 1/500000 inf])
    saveas(gcf,['D0=',num2str(D_0),', 1-D simulation 3D.png'])
    saveas(gcf,['D0=',num2str(D_0),', 1-D simulation 3D.fig'])
        %% 2D plot for all parameter space
    figure(2);clf;hold on;
    for i=1:length(delta_t_matrix_vik)
        subplot(2,4,i);hold on
        plot(AvN_plot,1./AratesV_plot(i,:),'kx')
        plot(AvN_plot,(1*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
        plot(AvN_plot,(2*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
%         plot(AvN_plot,(1*rateN_plot(i,:)),'r--')
%         plot(AvN_plot,(2*rateN_plot(i,:)),'r--')
        plot(AvN_plot,(1*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'r-.')
        plot(AvN_plot,(2*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'r-.')
        plot(AvN_plot,(1*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
        plot(AvN_plot,(2*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
        set(gca,'YScale','log')
        title(['\delta_t=',num2str(delta_t_matrix_vik(i))])
    end
    subplot(2,4,8);hold on;
    sim=plot(1,1,'kx')
    numeric=plot(1,1,'g-.')
    Kramer=plot(1,1,'r-.')
    Kramer_eff=plot(1,1,'b:')
    legend([sim numeric Kramer Kramer_eff],{'Simulation', 'Numeric FPE', 'D_\theta=4D_0/(R^2\theta_0^2))','D_\theta=2D_0/R^2'})
    suptitle(['1-D Simulation vs Kramer Escape Rate, D_0=',num2str(D_0)])
    saveas(gcf,['D0=',num2str(D_0),', 1-D simulation 2D.png'])

    %% plot for only one plot
%     figure(3);clf;hold on;
%     tau=0.1
%     
%    
%     title(['D_0=',num2str(D_0),', \tau=',num2str(tau)])
%     for i=1:length(delta_t_matrix_vik)
%         if delta_t_matrix_vik(i)==tau
%         plot(AvN_plot,1./AratesV_plot(i,:),'k-x')
%         plot(AvN_plot,(1*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'r--')
%         plot(AvN_plot,(2*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'r--')
%         plot(AvN_plot,(1*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
%         plot(AvN_plot,(2*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
%         plot(AvN_plot,(1*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
%         plot(AvN_plot,(2*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
%         end
%     end
%     set(gca,'YScale','log')
    
    %% 2D plot for all single plots
    close all;
    for i=1:length(delta_t_matrix_vik)
        %         subplot(2,4,i);hold on
        figure(10+i); clf;hold on;
        plot(AvN_plot,1./AratesV_plot(i,:),'kx')
        plot(AvN_plot,(1*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
        plot(AvN_plot,(2*rateN_plot(i,:)).*iF((-1 + AvN_plot)),'g-.')
        %         plot(AvN_plot,(1*rateN_plot(i,:)),'r--')
        %         plot(AvN_plot,(2*rateN_plot(i,:)),'r--')
        plot(AvN_plot,(1*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'r-.')
        plot(AvN_plot,(2*rateA_plot(i,:)).*iF((-1 + AvN_plot)),'r-.')
        plot(AvN_plot,(1*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
        plot(AvN_plot,(2*rateAeff_plot(i,:)).*iF((-1 + AvN_plot)),'b:')
        set(gca,'YScale','log')
        title(['1-D D_0=',num2str(D_0),', \delta t=',num2str(delta_t_matrix_vik(i))])
        saveas(gcf,['D0=',num2str(D_0),', 1-D, delta_t=',num2str(delta_t_matrix_vik(i)),'.png'])
    end
    %% Comparing the Effective Temperature and E_b
    figure(4);clf;hold on;title(['1-D Energy Barrier vs Thermal Fluctuation, D_0=',num2str(D_0)])
%     T_eff_matrix=2*AvN_plot/10^2;
    for i=1:length(delta_t_matrix_vik)
%                 plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),(1/12)/delta_t_matrix(i).*theta_plus_matrix(i,:).^4./(k_B*T_eff_matrix(i,:)),'-xg')
                plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,3*(AvN_plot-1).^2./(4*D_0./10^2.*delta_t_matrix_vik(i)),'-r')
                plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,3*(AvN_plot-1).^2./(2*D_0./10^2.*delta_t_matrix_vik(i).*AvN_plot.^2),'-b')
    end
    [x y] = meshgrid(-1:0.1:10); % Generate x and y data
    z = ones(size(x, 1))
    ; % Generate z data
    surf(x, y, z) % Plot the surface
    
%     Kramer_measured=plot(1,1,'-xg');
    Kramer=plot(1,1,'-r');
    Kramer_eff=plot(1,1,'-b');
    legend([  Kramer Kramer_eff],{'4D_0/(R^2\theta_0^2)','2D_0/R^2'})
    
    v = [5 2 5];
    [caz,cel] = view(v);
    grid on
    zlabel('E_b/(k_B T_{eff})')
    xlabel('\delta t')
    ylabel('\theta_0')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'Zscale','log')
    axis([-inf inf 1.05 pi/2 1e-2 inf])
    saveas(gcf,['D0=',num2str(D_0),', 1-D Energy Barrier.png'])
  
    %% Comparing the coefficients for the prefactor
%     figure(5);clf; hold on
%     %     for i=1:length(delta_t_matrix_vik)
%     for i=7
%         %         plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,AratesV_plot(i,:)./AratesV_plot(i,:),'kx')
%         plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:),'-.g')
%         plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateA_plot(i,:),'--r')
%         plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:),':b')
%         %         rateAeff_plot(i,:)).*iF((-1 + AvN_plot))
%         
%     end
% %     v = [5 2 5];
%         v=[1 0 0]
%     [caz,cel] = view(v);
%     grid on
%     xlabel('\delta t')
%     ylabel('\theta_0')
%     axis([0 inf 1 inf 0 3])
%             set(gca,'YScale','log')
%     %         set(gca,'ZScale','log')
%     %         title(['\delta_t=',num2str(delta_t_matrix_vik(i))])

%%







%% Part 2: Pin-Chuan's 2-d simulation
    clearvars -except D_0 file_name_pc file_name_vik
    % load(['2021.1.25_compare_PC,D_0=',num2str(D_0),'.mat'])
    load([file_name_pc,'.mat'])
    %% Additional Calculations and definitions
    [rateN_2D,rateN_eff_2D,rateN_measured_2D,rateA_2D,rateAeff_2D]=numeric_FPE(delta_t_matrix,theta_0_matrix,D_theta_matrix,D_0,R_matrix);
    iF = @(x) (sign(x)+1)./(sign(x)+1); % this function maps x to NaN for x<0, to 1 for x>0
            %% Simulations vs Kramer's Rate
            % Calculating the Kramer's Rate
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
    %%  plotting theta_0_matrix
    figure(101);clf;hold on
    for i=1:length(delta_t_matrix)
        plot(theta_0_1D_matrix,theta_0_matrix(i,:),'x')
    end
    yline(pi/2)
    %% Simulation vs Numerics
        %% Calculating the numeric solution for FPE
    
    figure(102);clf;hold on;
    subplot(1,2,1);hold on;title('Simulation vs Numerics')

        % Plotting Numerics vs Simulation
    for i=1:length(delta_t_matrix)
        trans_rate_sim_2D(i,:)=num_transitions_matrix(i,:)/time_duration(i);
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:),'xk-')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*rateN_measured_2D(i,:),'-.g')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*(2*rateN_measured_2D(i,:)),'-.g')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*rateN_2D(i,:),'-.r')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*(2*rateN_2D(i,:)),'-.r')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*rateN_eff_2D(i,:),':b')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*(2*rateN_eff_2D(i,:)),':b')
    end
    sim=plot(1,1,'x-k');
    Kramer_measured=plot(1,1,'-.g');
    Kramer=plot(1,1,'-.r');
    Kramer_eff=plot(1,1,':b');
    legend([sim Kramer_measured Kramer Kramer_eff],{'Simulations','D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
    
    v = [5 2 5];
    [caz,cel] = view(v);
    grid on
    
    xlabel('\delta t')
    ylabel('\theta_0')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'Zscale','log')
    axis([-inf inf 1.0 inf 1/min(time_duration) inf])
    
        %% Plotting all the 2-D 
        figure(103);clf;hold on;
        for i=1:length(delta_t_matrix)
            subplot(2,4,i);hold on
            plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'xk')
            plot(theta_0_matrix(i,:),rateN_measured_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.g')
            plot(theta_0_matrix(i,:),2*rateN_measured_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.g')
            plot(theta_0_matrix(i,:),rateN_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.r')
            plot(theta_0_matrix(i,:),2*rateN_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.r')
            plot(theta_0_matrix(i,:),rateN_eff_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),':b')
            plot(theta_0_matrix(i,:),2*rateN_eff_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),':b')
            set(gca,'YScale','log')
            axis([1.0 inf 1/min(time_duration) inf])
            title(['\delta_t=',num2str(delta_t_matrix(i))])
        end
        subplot(2,4,8);hold on;
        sim=plot(1,1,'xk');
        Kramer_measured=plot(1,1,'-.g');
        Kramer=plot(1,1,'-.r');
        Kramer_eff=plot(1,1,':b');
        legend([sim Kramer_measured Kramer Kramer_eff],{'Simulations','D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
        suptitle(['2-D Simulation vs FPE Numerics, D_0=',num2str(D_0)])
        saveas(gcf,['D0=',num2str(D_0),', 2-D simulation_numerics 2D.png'])
        %% plotting numeric solution for a single set 
%         tau=5
%         i=find(delta_t_matrix==tau);
%         figure(104);clf;hold on;
%         plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'xk')
%         plot(theta_0_matrix(i,:),rateN_measured_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.g')
%         plot(theta_0_matrix(i,:),2*rateN_measured_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.g')
%         plot(theta_0_matrix(i,:),rateN_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.r')
%         plot(theta_0_matrix(i,:),2*rateN_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.r')
%         plot(theta_0_matrix(i,:),rateN_eff_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),':b')
%         plot(theta_0_matrix(i,:),2*rateN_eff_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),':b')
%         set(gca,'YScale','log')
%         set(gca,'Xscale','log')
        
        for i=1:length(delta_t_matrix)
            figure(110+i);clf;hold on;
            sim=plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'xk')';
            Kramer_measured=plot(theta_0_matrix(i,:),rateN_measured_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.g');
            plot(theta_0_matrix(i,:),2*rateN_measured_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.g');
            Kramer=plot(theta_0_matrix(i,:),rateN_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.r');
            plot(theta_0_matrix(i,:),2*rateN_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),'-.r')
            Kramer_eff=plot(theta_0_matrix(i,:),rateN_eff_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),':b');
            plot(theta_0_matrix(i,:),2*rateN_eff_2D(i,:).*iF((-1 + theta_0_matrix(i,:))),':b')
            set(gca,'YScale','log')
            axis([1.0 inf 1/min(time_duration) inf])
            title(['2-D Numeric, D_0=',num2str(D_0),' \delta_t=',num2str(delta_t_matrix(i))])
            legend([sim Kramer_measured Kramer Kramer_eff],{'Simulations','D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'},'Location','Southwest')
            saveas(gcf,['D0=',num2str(D_0),', 2-D numeric, delta_t=',num2str(delta_t_matrix(i)),'.png'])
        end
       
    %% Simulation vs Kramer's rate

            %% Plotting
            figure(102);
            subplot(1,2,2);hold on;title('Simulation vs Kramer')
            for i=1:length(delta_t_matrix)
                plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:),'xk-')
                plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_measured(i,:),'-.g')
                plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),2*trans_rate_measured(i,:),'-.g')
                plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),transition_rate_Kramer(i,:),'-.r')
                plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),2*transition_rate_Kramer(i,:),'-.r')
                plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),transition_rate_Kramer_eff(i,:),':b')
                plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),2*transition_rate_Kramer_eff(i,:),':b')
            end
            sim=plot(1,1,'x-k');
            Kramer_measured=plot(1,1,'-.g');
            Kramer=plot(1,1,'-.r');
            Kramer_eff=plot(1,1,':b');
            legend([sim Kramer_measured Kramer Kramer_eff],{'Simulations','D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
            
            
            v = [5 2 5];
            [caz,cel] = view(v);
            grid on
            
            xlabel('\delta t')
            ylabel('\theta_0')
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            set(gca,'Zscale','log')
            axis([-inf inf 1.0 inf 1/min(time_duration) inf]) %% The minimum non-zero transition rate is 1/Obs_time_steps.
            suptitle(['2-D, D_0=',num2str(D_0)])
            saveas(gcf,['D0=',num2str(D_0),', 2-D simulation 3D.png'])
            %% Plotting all the 2-D
            figure(106);clf;hold on;
            for i=1:length(delta_t_matrix)
                subplot(2,4,i);hold on
                plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:),'xk')
                plot(theta_0_matrix(i,:),trans_rate_measured(i,:),'-.g')
                plot(theta_0_matrix(i,:),2*trans_rate_measured(i,:),'-.g')
                plot(theta_0_matrix(i,:),transition_rate_Kramer(i,:),'-.r')
                plot(theta_0_matrix(i,:),2*transition_rate_Kramer(i,:),'-.r')
                plot(theta_0_matrix(i,:),transition_rate_Kramer_eff(i,:),':b')
                plot(theta_0_matrix(i,:),2*transition_rate_Kramer_eff(i,:),':b')
                set(gca,'YScale','log')
                set(gca,'Xscale','log')
                set(gca,'YScale','log')
                axis([1.0 inf 1/min(time_duration) inf])
                title(['\delta_t=',num2str(delta_t_matrix(i))])
            end
            subplot(2,4,8);hold on;
            sim=plot(1,1,'xk');
            Kramer_measured=plot(1,1,'-.g');
            Kramer=plot(1,1,'-.r');
            Kramer_eff=plot(1,1,':b');
            legend([sim Kramer_measured Kramer Kramer_eff],{'Simulations','D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
            suptitle(['2-D Simulation vs Kramers Escape Rate, D_0=',num2str(D_0)])
            saveas(gcf,['D0=',num2str(D_0),', 2-D simulation_Kramer 2D.png'])
        %% Plotting for 2D simulations for one specific data set
%             tau=10
%             i=find(delta_t_matrix==tau)
%             figure(107);clf;hold on;
%             plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:),'xk')
%             plot(theta_0_matrix(i,:),trans_rate_measured(i,:),'-.g')
%             plot(theta_0_matrix(i,:),2*trans_rate_measured(i,:),'-.g')
%             plot(theta_0_matrix(i,:),transition_rate_Kramer(i,:),'-.r')
%             plot(theta_0_matrix(i,:),2*transition_rate_Kramer(i,:),'-.r')
%             plot(theta_0_matrix(i,:),transition_rate_Kramer_eff(i,:),':b')
%             plot(theta_0_matrix(i,:),2*transition_rate_Kramer_eff(i,:),':b')
%             set(gca,'YScale','log')
%             set(gca,'Xscale','log')

            for i=1:length(delta_t_matrix)
                figure(120+i);clf;hold on
                sim=plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:),'xk');
                Kramer_measured=plot(theta_0_matrix(i,:),trans_rate_measured(i,:),'-.g');
                plot(theta_0_matrix(i,:),2*trans_rate_measured(i,:),'-.g')
                Kramer=plot(theta_0_matrix(i,:),transition_rate_Kramer(i,:),'-.r');
                plot(theta_0_matrix(i,:),2*transition_rate_Kramer(i,:),'-.r');
                Kramer_eff=plot(theta_0_matrix(i,:),transition_rate_Kramer_eff(i,:),':b');
                plot(theta_0_matrix(i,:),2*transition_rate_Kramer_eff(i,:),':b');
                set(gca,'YScale','log')
                set(gca,'Xscale','log')
                set(gca,'YScale','log')
                axis([1.0 inf 1/min(time_duration) inf])
                title(['2-D Kramer, D_0=',num2str(D_0),' \delta_t=',num2str(delta_t_matrix(i))])
                legend([sim Kramer_measured Kramer Kramer_eff],{'Simulations','D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
                saveas(gcf,['D0=',num2str(D_0),', 2-D Kramer, delta_t=',num2str(delta_t_matrix(i)),'.png'])
            end
          
    %% Bifurcations
    figure(108);clf;hold on;title(['2-D Bifurcation diagram, D_0=',num2str(D_0)])
%     subplot(1,2,2);hold on
    for i=1:length(delta_t_matrix)
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),theta_plus_matrix(i,:),'xk')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),theta_minus_matrix(i,:),'xk')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),+sqrt(6./theta_0_matrix(i,:).*(theta_0_matrix(i,:)-1)),'b-')
        plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),-sqrt(6./theta_0_matrix(i,:).*(theta_0_matrix(i,:)-1)),'r-')
    end
    sim=plot3(1,1,1,'xk');
    blue=plot3(1,1,1,'b-');
    red=plot3(1,1,1,'r-');
    legend([sim blue red],{'Simulations','\theta_+ = sqrt{(6(\theta_0-1)/\theta_0)}','-\theta_+'},'Location','Northwest')
    v = [5 2 5];
    v=[ 1 0 0]
    [caz,cel] = view(v);
    grid on
    
    xlabel('\delta t')
    ylabel('\theta_0')
    set(gca,'XScale','log')
    saveas(gcf,['D0=',num2str(D_0),', 2-D Bifurcation.png'])

    %% Numerics vs Kramer
        %% Kramer's rate (original: D_theta=4*D_0/(R^2 theta_0^2)
        figure(109);clf;hold on;
        subplot(2,2,1); hold on;
        title('D_{\theta}=4D_0/(R^2 \theta_0^2)')
        for i=1:length(delta_t_matrix)
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_measured(i,:),'x-r')
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),2*trans_rate_measured(i,:),'x-r')
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*rateN_measured_2D(i,:),'-.b')
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*(2*rateN_measured_2D(i,:)),'-.b')
        end
%         FPE=plot(1,1,'x-r')
%         Kramer=plot(1,1,'-.b')
%         legend([FPE Kramer],{'FPE Numerics','Kramer'})
        
        v = [5 2 5];
        [caz,cel] = view(v);
        grid on
        
        
        xlabel('\delta t')
        ylabel('\theta_0')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca,'Zscale','log')
        axis([-inf inf 1.0 inf 1/min(time_duration) inf]) %% The minimum non-zero transition rate is 1/Obs_time_steps.
        %% Kramer_eff's rate (D_theta=2D_0/R^2
        %         figure(110);clf;hold on;
        subplot(2,2,2); hold on;title('D_{\theta}=2D_0/R^2')
        for i=1:length(delta_t_matrix)
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),transition_rate_Kramer(i,:),'x-r')
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),2*transition_rate_Kramer(i,:),'x-r')
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*rateN_2D(i,:),'-.b')
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*(2*rateN_2D(i,:)),'-.b')
        end
%         FPE=plot(1,1,'x-r')
%         Kramer=plot(1,1,'-.b')
%         legend([FPE Kramer],{'FPE Numerics','Kramer'})
        
        v = [5 2 5];
        [caz,cel] = view(v);
        grid on
        
        xlabel('\delta t')
        ylabel('\theta_0')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca,'Zscale','log')
        axis([-inf inf 1.0 inf 1/min(time_duration) inf]) %% The minimum non-zero transition rate is 1/Obs_time_steps.
        %% Kramer_measured (D_theta=D_eff)title('D_{\theta}=D_{eff}')
        %         figure(111);clf;hold on;
        subplot(2,2,3);hold on;title('D_{\theta}=D_{eff}')
        for i=1:length(delta_t_matrix)
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_measured(i,:),'x-r')
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),2*trans_rate_measured(i,:),'x-r')
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*rateN_measured_2D(i,:),'-.b')
            plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),iF((-1 + theta_0_matrix(i,:))).*(2*rateN_measured_2D(i,:)),'-.b')
        end
%         FPE=plot(1,1,'x-r')
%         Kramer=plot(1,1,'-.b')
%         legend([FPE Kramer],{'FPE Numerics','Kramer'})
        
        v = [5 2 5];
        [caz,cel] = view(v);
        grid on
        
        xlabel('\delta t')
        ylabel('\theta_0')
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        set(gca,'Zscale','log')
        axis([-inf inf 1.0 inf 1/min(time_duration) inf]) %% The minimum non-zero transition rate is 1/Obs_time_steps.
        %%
        subplot(2,2,4);hold on;
        FPE=plot(1,1,'x-r')
        Kramer=plot(1,1,'-.b')
        legend([FPE Kramer],{'FPE Numerics','Kramer'})
        suptitle(['2-D: Numerics vs Kramer Escape Rate, D_0=',num2str(D_0)])
        saveas(gcf,['D0=',num2str(D_0),', 2-D Numeric_Kramer 3D.png'])
    %% Comparing the Effective Temperature and E_b
    T_eff_matrix=D_theta_matrix;
    figure(110);clf;hold on;title(['2-D Energy Barrier vs Thermal Fluctuation, D_0=',num2str(D_0)]);
    for i=1:length(delta_t_matrix)
%                 plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),(1/12)./delta_t_matrix(i).*theta_plus_matrix(i,:).^4,'xk-')
%                 plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),k_B*T_eff_matrix(i,:),'-.g')

                plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),(1/12)/delta_t_matrix(i).*theta_plus_matrix(i,:).^4./(k_B*T_eff_matrix(i,:)),'-xg')
                plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),3*(theta_0_matrix(i,:)-1).^2./(4*D_0./R_matrix(i,:).^2.*delta_t_matrix(i)),'-r')
                plot3(ones(length(theta_0_1D_matrix))*delta_t_matrix(i),theta_0_matrix(i,:),3*(theta_0_matrix(i,:)-1).^2./(2*D_0./R_matrix(i,:).^2.*delta_t_matrix(i).*theta_0_matrix(i,:).^2),'-b')
    end
    [x y] = meshgrid(-1:0.1:10); % Generate x and y data
    z = ones(size(x, 1)); % Generate z data
    surf(x, y, z) % Plot the surface
    
    Kramer_measured=plot(1,1,'-xg');
    Kramer=plot(1,1,'-r');
    Kramer_eff=plot(1,1,'-b');
    legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
    
    v = [5 2 5];
    [caz,cel] = view(v);
    grid on
%     title('E_b/(k_B T_{eff})')
    zlabel('E_b/(k_B T_{eff})')
    xlabel('\delta t')
    ylabel('\theta_0')
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'Zscale','log')
    axis([-inf inf 1.05 pi/2 1e-3 inf])
    saveas(gcf,['D0=',num2str(D_0),', 2-D Energy Barrier.png'])
    

        
    
    
