%% Comparing the coefficients for the prefactor with all different D_0's

close all
clear

single_plots='off'

% D_0_matrix=[1:15,20] % for 1-D
% D_0_matrix=[1 2 5 10 15]
D_0_matrix=[1 3 4 5 6 7 8 10 11 12 13 15] % for 1-D
part1='on'
part2='off'


% D_0_matrix=[1 3 4 5 6 7 8 10 11 12 13 15] % for 2-D
% % D_0_matrix=[15]
% part1='off'
% part2='on'


%% For fitting the prefactors
Kramers_ratio1D=[];
theta_ratio1D=[];
D_delta_t_fit1D=[];
theta_fit1D=[];
prefactor_fit1D=[];

Kramers_ratio2D=[];
theta_ratio2D=[];
D_delta_t_fit2D=[];
theta_fit2D=[];
prefactor_fit2D=[];
%%
for D_0_index=1:length(D_0_matrix)
%     dt=0.01;
    D_0=D_0_matrix(D_0_index);
%     delta_t_matrix_vik=[0.1 0.2 0.5 1 2 5 10]  %% Make 1-D simulation
%     delta_t_matrix=[0.1 0.2 0.5 1 2 5 10] %% for 2-D simulation

%     theta_0_1D_matrix=(logspace(log10(1),log10(3),30));
%     theta_0_1D_matrix=(logspace(log10(1),log10(3),15))
    
    
    
%     file_name_pc=['2021.1.28_scan_pc_D_0=',num2str(D_0)];
%     file_name_vik=['2021.1.28_scan_vik_D_0=',num2str(D_0)];
    file_name_pc=['2021.2.18_scan_pc_D_0=',num2str(D_0)];
    file_name_vik=['2021.2.17_scan_vik_D_0=',num2str(D_0)];
    
    %% Determining figure numbers for D_0*delta_t 
%     D_times_delta_t_matrix=reshape(D_0_matrix.*delta_t_matrix',[length(D_0_matrix)*length(delta_t_matrix),1]);
%     D_times_delta_t_matrix=sort(unique(D_times_delta_t_matrix));
%     D_times_delta_t_matrix_vik=reshape(D_0_matrix.*delta_t_matrix_vik',[length(D_0_matrix)*length(delta_t_matrix_vik),1]);
%     D_times_delta_t_matrix_vik=sort(unique(D_times_delta_t_matrix));
    
    
%% Part 1    
    switch part1
        case 'on'
    %% Loading for part 1
%     tV_max=50000; % set to be 50000 in the previous simulations, could be loaded from mat file on the next file
    load([file_name_vik,'.mat'])
    delta_t_matrix_vik=reshape(tau_matrix,[1 length(tau_matrix)]);
    theta_0_1D_matrix=AvN_plot;
    %% Determining figure numbers for D_0*delta_t 
    D_times_delta_t_matrix_vik=reshape(D_0_matrix.*delta_t_matrix_vik',[length(D_0_matrix)*length(delta_t_matrix_vik),1]);
    D_times_delta_t_matrix_vik=sort(unique(D_times_delta_t_matrix_vik));
    %% 1-D (3D)
        
        figure(1);hold on;
        for i=1:length(delta_t_matrix_vik)
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:),'-.g')
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateA_plot(i,:),'--r')
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:),':b')
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateNeff_plot(i,:),':m')

            
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:),'--b')
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateNeff_plot(i,:),'.m')
        end
        if D_0_index==length(D_0_matrix)
        v = [5 2 5];
        [caz,cel] = view(v);
        grid on;
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        axis([0 inf 1 inf 0 3])
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        title('1-D ratio')
%         sim=plot(1,1,'kx')
        numeric=plot(1,1,'g-.');
        Kramer=plot(1,1,'r-.');
        Kramer_eff=plot(1,1,'b:');
        numeric_eff=plot(1,1,'m:');
        legend([ numeric Kramer Kramer_eff,numeric_eff],{'Numeric FPE, D_\theta=4D_0/(R^2\theta_0^2)', 'Kramers, D_\theta=4D_0/(R^2\theta_0^2))','Kramers, D_\theta=2D_0/R^2','Numeric FPE, D_\theta=2D_0/R^2'},'Location','southwest')
        saveas(gcf,['1-D prefactor 3D.png'])
    
        end
        %% 1-D (3D) with mask with both D_theta's
        figure(97);hold on;
        for i=1:length(delta_t_matrix_vik)
            error=5*10^-2;
%             cutoff_mask(i,:)=cutoff(1./AratesV_plot(i,:),1/tV_max);
            cutoff_mask(i,:)=cutoff(rateNeff_plot(i,:),1/tV_max);
            pre_mask_3D(i,:)=pre_filter(rateA_plot(i,:),rateN_plot(i,:),error);
            pre_mask_3Deff(i,:)=pre_filter(rateAeff_plot(i,:),rateNeff_plot(i,:),error);
            mask_3D(i,:)=pre_mask_3D(i,:).*cutoff_mask(i,:);
            mask_3Deff(i,:)=pre_mask_3Deff(i,:).*cutoff_mask(i,:);
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:).*mask_3D(i,:),'-.g')
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateA_plot(i,:).*mask_3D(i,:),'--r')
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:).*mask_3Deff(i,:),':b')
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateNeff_plot(i,:).*mask_3Deff(i,:),':m')
        end
        if D_0_index==length(D_0_matrix)
        v = [5 2 5];
%         v = [1 0 0]
        [caz,cel] = view(v);
        grid on;
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        axis([0 inf 1 inf 0 3])
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        title(['1-D ratio (masked with error=',num2str(error),')'])
%         sim=plot(1,1,'kx')
        numeric=plot(1,1,'g-.');
        Kramer=plot(1,1,'r-.');
        Kramer_eff=plot(1,1,'b:');
        numeric_eff=plot(1,1,'m:');
        legend([ numeric Kramer Kramer_eff,numeric_eff],{'Numeric FPE, D_\theta=4D_0/(R^2\theta_0^2)', 'Kramers, D_\theta=4D_0/(R^2\theta_0^2))','Kramers, D_\theta=2D_0/R^2','Numeric FPE, D_\theta=2D_0/R^2'},'Location','southwest')
        saveas(gcf,['1-D prefactor 3D Masked.png'])
        end
        %% 1-D (3D) with mask with both D_theta's and hide theta_0(i)>1.3+0.4*log(D_0*delta_t)/log(10)
        figure(96);hold on;
        for i=1:length(delta_t_matrix_vik)
            error=5*10^-2;
%             cutoff_mask(i,:)=cutoff(1./AratesV_plot(i,:),1/tV_max);
%             cutoff_mask(i,:)=cutoff(rateNeff_plot(i,:),1/tV_max);
%             pre_mask_3D(i,:)=pre_filter(rateA_plot(i,:),rateN_plot(i,:),error);
            mask_diag1D(i,:)=diagonal1D(AvN_plot,delta_t_matrix_vik(i),D_0);
            mask_3D_full(i,:)=mask_3D(i,:).*mask_diag1D(i,:);
            mask_3D_fulleff(i,:)=mask_3Deff(i,:).*mask_diag1D(i,:);
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:).*mask_3D_full(i,:),'-.g')
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateA_plot(i,:).*mask_3D_full(i,:),'--r')
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:).*mask_3D_full(i,:),':b')
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateNeff_plot(i,:).*mask_3D_full(i,:),':m')
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:).*mask_3D_full(i,:),'-og')
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateA_plot(i,:).*mask_3D_full(i,:),'-or')
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:).*mask_3D_fulleff(i,:),':ob')
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateNeff_plot(i,:).*mask_3D_fulleff(i,:),':om')
        
            %% Appending to Fit the ratios, 2D fit
            Kramers_ratio1D=[Kramers_ratio1D 1./AratesV_plot(i,:)./rateAeff_plot(i,:).*mask_3D_fulleff(i,:)];
            theta_ratio1D=[theta_ratio1D AvN_plot]
            %% 3D fit
            D_delta_t_fit1D=[D_delta_t_fit1D ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix_vik(i)];
            theta_fit1D=[theta_fit1D AvN_plot];
            prefactor_fit1D=[prefactor_fit1D 1./AratesV_plot(i,:)./rateAeff_plot(i,:).*mask_3D_fulleff(i,:)]
            
            
        end
        if D_0_index==length(D_0_matrix)
        v = [5 2 5];
%         v = [1 0 0]
        [caz,cel] = view(v);
        grid on;
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        axis([0 inf 1 inf 0 3])
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        title(['1-D ratio (masked with error=',num2str(error),'), \theta_0>1.3+0.4*log(D_0*\delta t)/log(10) is hidden'])
        numeric=plot(1,1,'g-.');
        Kramer=plot(1,1,'r-.');
        Kramer_eff=plot(1,1,'b:');
        numeric_eff=plot(1,1,'m:');
        legend([ numeric Kramer Kramer_eff,numeric_eff],{'Numeric FPE, D_\theta=4D_0/(R^2\theta_0^2)', 'Kramers, D_\theta=4D_0/(R^2\theta_0^2))','Kramers, D_\theta=2D_0/R^2','Numeric FPE, D_\theta=2D_0/R^2'},'Location','southwest')
        saveas(gcf,['1-D prefactor 3D Masked, trimmed.png'])
        end
        %% 1-D (3D) with mask with D_theta = 4D/theta_0^2
        
        figure(98);hold on;
        for i=1:length(delta_t_matrix_vik)
            error=5*10^-2;
            pre_mask(i,:)=pre_filter(rateA_plot(i,:),rateN_plot(i,:),error);
            mask(i,:)=pre_mask(i,:).*cutoff_mask(i,:).*mask_diag1D(i,:);
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:).*mask(i,:),'-.g')
            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateA_plot(i,:).*mask(i,:),'--r')
        end
        if D_0_index==length(D_0_matrix)
%         v = [5 2 5];
        v = [1 0 0];
        [caz,cel] = view(v);
        grid on;
        
        ylabel('\theta_0')
        axis([0 inf 1 inf 0 3])
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        title(['1-D ratio (masked with error=',num2str(error),'), D=4D_0/(R^2\theta_0^2)'])
%         sim=plot(1,1,'kx')
        numeric=plot(1,1,'g-.');
        Kramer=plot(1,1,'r-.');
        Kramer_eff=plot(1,1,'b:');
        numeric_eff=plot(1,1,'m:');
        legend([ numeric Kramer Kramer_eff,numeric_eff],{'Numeric FPE, D_\theta=4D_0/(R^2\theta_0^2)', 'Kramers, D_\theta=4D_0/(R^2\theta_0^2))','Kramers, D_\theta=2D_0/R^2','Numeric FPE, D_\theta=2D_0/R^2'},'Location','southwest')
        saveas(gcf,['1-D prefactor 3D Masked, D=4D_theta_0.png'])
        end
                %% 1-D (3D) with mask with D_theta = 2D
       
        figure(99);hold on;
        for i=1:length(delta_t_matrix_vik)
            error=5*10^-2;
            pre_mask_eff(i,:)=pre_filter(rateAeff_plot(i,:),rateNeff_plot(i,:),error);
            mask_eff(i,:)=pre_mask_eff(i,:).*cutoff_mask(i,:).*mask_diag1D(i,:);
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:).*mask_eff(i,:),':b')
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateNeff_plot(i,:).*mask_eff(i,:),':m')

            plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:).*mask_eff(i,:),'-ob')
%             plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateNeff_plot(i,:).*mask_eff(i,:),'-om')
        end
        if D_0_index==length(D_0_matrix)
%         v = [5 2 5];
        v = [1 0 0]
        [caz,cel] = view(v);
        grid on;
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        axis([0 inf 1 inf 0 3])
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        title(['1-D ratio (masked with error=',num2str(error),'), D=2D_0/R^2'])
%         sim=plot(1,1,'kx')
        numeric=plot(1,1,'g-.');
        Kramer=plot(1,1,'r-.');
        Kramer_eff=plot(1,1,'b:');
        numeric_eff=plot(1,1,'m:');
        legend([ numeric Kramer Kramer_eff,numeric_eff],{'Numeric FPE, D_\theta=4D_0/(R^2\theta_0^2)', 'Kramers, D_\theta=4D_0/(R^2\theta_0^2))','Kramers, D_\theta=2D_0/R^2','Numeric FPE, D_\theta=2D_0/R^2'},'Location','southwest')
        saveas(gcf,['1-D prefactor 3D Masked, D=2D.png'])
        end
    %% Plotting 1-D with one specific
    switch single_plots
        case 'on'
            %% 1-D Kramer's Rate 2D plots
        for i=1:length(delta_t_matrix_vik)
            j=find(D_0*delta_t_matrix_vik(i)==D_times_delta_t_matrix_vik);
            figure(10+j);hold on;
    %         plot(AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:),'-.g')
            plot(AvN_plot,1./AratesV_plot(i,:)./rateA_plot(i,:),'--r')
            plot(AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:),':b')
        end
            %% 1-D Numerics Rate 2D plots
        for i=1:length(delta_t_matrix_vik)
            j=find(D_0*delta_t_matrix_vik(i)==D_times_delta_t_matrix_vik);
            figure(50+j);hold on;
    %         plot(AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:),'-.g')
            plot(AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:),'--r')
            plot(AvN_plot,1./AratesV_plot(i,:)./rateNeff_plot(i,:),':b')
        end
    end
    end
            
            
            
            
%% Part 2
    switch part2
        case 'on'
    load([file_name_pc,'.mat'])
        %% Determining figure numbers for D_0*delta_t 
        D_times_delta_t_matrix=reshape(D_0_matrix.*delta_t_matrix,[length(D_0_matrix)*length(delta_t_matrix),1]);
        D_times_delta_t_matrix=sort(unique(D_times_delta_t_matrix));
        %% Additional Calculations and definitions
        [rateN_2D,rateN_eff_2D,rateN_measured_2D,rateA_2D,rateAeff_2D]=numeric_FPE(delta_t_matrix,theta_0_matrix,D_theta_matrix,D_0,R_matrix);
        % rateA_2D ~= transition_rate_Kramer; rateAeff_2D ~= transition_rate_Kramer_eff
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

        %% 2-D Kramers rate 3D plot
        figure(2);hold on
        for i=1:length(delta_t_matrix) 
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(trans_rate_measured(i,:).*iF(-1+theta_0_matrix(i,:))),'-.g')
%             plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer(i,:).*iF(-1+theta_0_matrix(i,:))),'-.r')
%             plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer_eff(i,:).*iF(-1+theta_0_matrix(i,:))),':b')
        end
        if D_0_index==length(D_0_matrix)
        v = [5 2 5];
        view(v);
        grid on
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        set(gca,'YScale','log')
        set(gca,'XScale','log')
        set(gca,'Zscale','linear')
        title('2-D Kramer')
        axis([-inf inf 1.0 1.6 1/min(time_duration) 5]) %% The minimum non-zero transition rate is 1/Obs_time_steps.

        Kramer_measured=plot(1,1,'-.g');
        Kramer=plot(1,1,'-.r');
        Kramer_eff=plot(1,1,':b');
        legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
        saveas(gcf,['2-D prefactor Kramers 3D.png'])
        end
        %% 2-D Kramers rate measured 3D plot
        figure(7);hold on
        for i=1:length(delta_t_matrix) 
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(trans_rate_measured(i,:).*iF(-1+theta_0_matrix(i,:))),'-.k')
%             plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer(i,:).*iF(-1+theta_0_matrix(i,:))),'-.r')
%             plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer_eff(i,:).*iF(-1+theta_0_matrix(i,:))),':b')
        end
        if D_0_index==length(D_0_matrix)
        v = [5 2 5];
        view(v);
        grid on
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        set(gca,'YScale','log')
        set(gca,'XScale','log')
        set(gca,'Zscale','linear')
        title('2-D Kramer, D_\theta=D_{eff}')
        axis([-inf inf 1.0 1.6 1/min(time_duration) 5]) %% The minimum non-zero transition rate is 1/Obs_time_steps.

%         Kramer_measured=plot(1,1,'-.g');
%         Kramer=plot(1,1,'-.r');
%         Kramer_eff=plot(1,1,':b');
%         legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
        saveas(gcf,['2-D prefactor Kramers 3D measured.png'])
        end
        %% 2-D Kramers rate 3D plot with mask
        figure(4);hold on
        for i=1:length(delta_t_matrix) 
%             error=10^-1;
            error=5*10^-2
            cutoff_mask(i,:)=cutoff(transition_rate_Kramer_eff(i,:),1/time_duration(1));
            cutoff_mask(i,:)=cutoff(trans_rate_sim_2D(i,:),1/time_duration(1));
%             cutoff_mask(i,:)=(transition_rate_Kramer_eff(i,:)>1/time_duration(1)); % masking for rate > 1/time_duration
            pre_mask_2Kramers(i,:)=pre_filter(rateN_eff_2D(i,:),transition_rate_Kramer_eff(i,:),error);
            mask_diag2D(i,:)=diagonal2D(theta_0_matrix(i,:),delta_t_matrix(i),D_0);
            mask_2Kramers(i,:)=pre_mask_2Kramers(i,:).*cutoff_mask(i,:).*mask_diag2D(i,:);
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(trans_rate_measured(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2Kramers(i,:),'-og')
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2Kramers(i,:),'-or')
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer_eff(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2Kramers(i,:),'-ob')
        end
        if D_0_index==length(D_0_matrix)
%         v = [5 2 5];
        v = [1 0 0];
        view(v);
        grid on
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        set(gca,'YScale','log')
        set(gca,'XScale','log')
        set(gca,'Zscale','linear')
        title('2-D Kramer masked')
        axis([-inf inf 1.0 1.6 1/min(time_duration) 5]) %% The minimum non-zero transition rate is 1/Obs_time_steps.

        Kramer_measured=plot(1,1,'-.g');
        Kramer=plot(1,1,'-.r');
        Kramer_eff=plot(1,1,':b');
        legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
        saveas(gcf,['2-D prefactor Kramers 3D masked.png'])
        end
        %% 2-D Kramers rate 3D plot (trans_rate_measured) with mask
        figure(6);hold on
        for i=1:length(delta_t_matrix) 
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(trans_rate_measured(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2Kramers(i,:),'-ok')

            %% Appending to Fit the ratios, 2D fit
           
%             Kramers_ratio2D=[Kramers_ratio2D trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer_eff(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2Kramers(i,:)]; %% D_theta=2D
            Kramers_ratio2D=[Kramers_ratio2D trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(trans_rate_measured(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2Kramers(i,:)]; % D_theta=D_eff
            theta_ratio2D=[theta_ratio2D theta_0_matrix(i,:)]
            size(Kramers_ratio2D)
            size(theta_ratio2D)

            
            %% 3D fit
            D_delta_t_fit2D=[D_delta_t_fit2D ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i)];
            theta_fit2D=[theta_fit2D theta_0_matrix(i,:)];

%             prefactor_fit2D=[prefactor_fit2D trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer_eff(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2Kramers(i,:)] %% D_theta=2D
            prefactor_fit2D=[prefactor_fit2D trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(trans_rate_measured(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2Kramers(i,:)] % D_theta=D_eff 
        end
        if D_0_index==length(D_0_matrix)
        v = [5 2 5];
%         v = [1 0 0];
        view(v);
        grid on
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        set(gca,'YScale','log')
        set(gca,'XScale','log')
        set(gca,'Zscale','linear')
        title('2-D Kramer masked, D_\theta=D_{eff}')
        axis([-inf inf 1.0 1.6 1/min(time_duration) 5]) %% The minimum non-zero transition rate is 1/Obs_time_steps.

%         Kramer_measured=plot(1,1,'-.g');
%         Kramer=plot(1,1,'-.r');
%         Kramer_eff=plot(1,1,':b');
%         legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
        saveas(gcf,['2-D prefactor Kramers 3D masked measured.png'])
        end
        
            %% 2-D Kramer's rate 2D plots
            switch single_plots
                case 'on'
            for i=1:length(delta_t_matrix)
                j=find(D_0*delta_t_matrix(i)==D_times_delta_t_matrix);
                figure(100+j);hold on;
                plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(trans_rate_measured(i,:).*iF(-1+theta_0_matrix(i,:))),'-.g')
                plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer(i,:).*iF(-1+theta_0_matrix(i,:))),'-.r')
                plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer_eff(i,:).*iF(-1+theta_0_matrix(i,:))),':b')
%                 if D_0_index==length(D_0_matrix)
%                     xlabel('\theta_0')
%                     ylabel('ratio')
%                     title(['2-D prefactor sim/Kramer, D_0 \delta t=',num2str(D_times_delta_t_matrix(j))])
%                     axis([1 inf 0 5])
%                     Kramer_measured=plot(1,1,'-.g');
%                     Kramer=plot(1,1,'-.r');
%                     Kramer_eff=plot(1,1,':b');
%                     legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'},'Location','southwest')
%                     saveas(gcf,['2-D prefactor Kramers, D_0 delta t=',num2str(D_times_delta_t_matrix(j)),'.png'])
%                 end
            end
            end
        %% 2-D Numerics rate 3D plot
                figure(3);hold on
        for i=1:length(delta_t_matrix)
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_measured_2D(i,:).*iF(-1+theta_0_matrix(i,:))),'-.g')
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_2D(i,:).*iF(-1+theta_0_matrix(i,:))),'-.r')
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_eff_2D(i,:).*iF(-1+theta_0_matrix(i,:))),':b')
            
        end
        if D_0_index==length(D_0_matrix)
        v = [5 2 5];
        view(v);
        grid on
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        set(gca,'YScale','log')
        set(gca,'XScale','log')
        set(gca,'Zscale','linear')
        title('2-D Numeric')
        axis([-inf inf 1.0 1.6 1/min(time_duration) 5]) %% The minimum non-zero transition rate is 1/Obs_time_steps.
        Kramer_measured=plot(1,1,'-.g');
        Kramer=plot(1,1,'-.r');
        Kramer_eff=plot(1,1,':b');
        legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
        saveas(gcf,['2-D prefactor Numeric 3D.png'])
        end
        %% 2-D Numerics rate 3D plot with mask
        figure(5);hold on
        for i=1:length(delta_t_matrix)
            error=5*10^-2;
            pre_mask_2numeric(i,:)=pre_filter(rateN_eff_2D(i,:),transition_rate_Kramer_eff(i,:),error);
%             mask_diag2D(i,:)=diagonal2D(theta_0_matrix(i,:),delta_t_matrix(i),D_0);
            mask_2numeric(i,:)=pre_mask_2numeric(i,:).*cutoff_mask(i,:).*mask_diag2D(i,:);
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_measured_2D(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2numeric(i,:),'-.g')
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_2D(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2numeric(i,:),'-.r')
            plot3(ones(size(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_eff_2D(i,:).*iF(-1+theta_0_matrix(i,:))).*mask_2numeric(i,:),':b')
            
        end
        
        if D_0_index==length(D_0_matrix)
            %         v = [5 2 5];
        v = [1 0 0]
        view(v);
        grid on
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        set(gca,'YScale','log')
        set(gca,'XScale','log')
        set(gca,'Zscale','linear')
        title('2-D Numeric masked')
        axis([-inf inf 1.0 1.6 1/min(time_duration) 5]) %% The minimum non-zero transition rate is 1/Obs_time_steps.
        Kramer_measured=plot(1,1,'-.g');
        Kramer=plot(1,1,'-.r');
        Kramer_eff=plot(1,1,':b');
        legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
        saveas(gcf,['2-D prefactor Numeric 3D masked.png'])
        end
            %% 2-D Numeric's rate 2D plots
            switch single_plots
                case 'on'
            for i=1:length(delta_t_matrix)
                j=find(D_0*delta_t_matrix(i)==D_times_delta_t_matrix);
                figure(200+j);hold on;
                plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_measured_2D(i,:).*iF(-1+theta_0_matrix(i,:))),'-.g')
                plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_2D(i,:).*iF(-1+theta_0_matrix(i,:))),'-.r')
                plot(theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_eff_2D(i,:).*iF(-1+theta_0_matrix(i,:))),':b')
%                 if D_0_index==length(D_0_matrix)
%                     xlabel('\theta_0')
%                     ylabel('ratio')
%                     title(['2-D ratio sim/Numeric, D_0 \delta t=',num2str(D_times_delta_t_matrix(j))])
%                     axis([1 inf 0 5])
%                     Kramer_measured=plot(1,1,'-.g');
%                     Kramer=plot(1,1,'-.r');
%                     Kramer_eff=plot(1,1,':b');
%                     legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'},'Location','northwest')
%                     saveas(gcf,['2-D prefactor Numeric, D_0 delta t=',num2str(D_times_delta_t_matrix(j)),'.png'])
%                 end
            end
            end
    end
end








%% Legending and Saving (The commented codes above don't cover all j's, so we need this section of codes)
switch single_plots
    case 'on'
% for j=1:length(D_times_delta_t_matrix)
    switch part1
        case 'on'
           for j=1:length(D_times_delta_t_matrix_vik)
    %% 1-D Kramer's Rate 2D plots
    figure(10+j)
    xlabel('\theta_0')
    ylabel('ratio')
    title(['1-D prefactor sim/Kramers, D_0 \delta t=',num2str(D_times_delta_t_matrix_vik(j))])
    axis([1 inf 0 5])
    %     numeric=plot(1,1,'g-.');
    Kramer=plot(1,1,'r-.');
    Kramer_eff=plot(1,1,'b:');
    
    legend([Kramer Kramer_eff],{'Kramers, D_\theta=4D_0/(R^2\theta_0^2))','Kramers, D_\theta=2D_0/R^2'},'Location','northwest')
    saveas(gcf,['1-D prefactor Kramers, D_0 delta t=',num2str(D_times_delta_t_matrix_vik(j)),'.png'])

    %% 1-D Numerics Rate 2D plots
    figure(50+j)
    xlabel('\theta_0')
    ylabel('ratio')
    title(['1-D prefactor sim/Numeric, D_0 \delta t=',num2str(D_times_delta_t_matrix_vik(j))])
    axis([1 inf 0 5])
    %     numeric=plot(1,1,'g-.');
    numeric=plot(1,1,'r-.');
    numeric_eff=plot(1,1,'b:');
    %     legend([ numeric Kramer Kramer_eff],{'Numeric FPE, D_\theta=4D_0/(R^2\theta_0^2)', 'Kramers, D_\theta=4D_0/(R^2\theta_0^2))','Kramers, D_\theta=2D_0/R^2'},'Location','northwest')
    
    legend([numeric numeric_eff],{'Numeric FPE, D_\theta=4D_0/(R^2\theta_0^2))','Numeric FPE, D_\theta=2D_0/R^2'},'Location','northwest')
    saveas(gcf,['1-D prefactor Numeric, D_0 delta t=',num2str(D_times_delta_t_matrix_vik(j)),'.png'])
           end
    end
    
    
    
    %%
    switch part2
        case 'on'
            for j=1:length(D_times_delta_t_matrix)
    %% 2-D Kramer's rate 2D plots
    figure(100+j);hold on;
    xlabel('\theta_0')
    ylabel('ratio')
    title(['2-D prefactor sim/Kramer, D_0 \delta t=',num2str(D_times_delta_t_matrix(j))])
    axis([1 inf 0 5])
    Kramer_measured=plot(1,1,'-.g');
    Kramer=plot(1,1,'-.r');
    Kramer_eff=plot(1,1,':b');
    legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'},'Location','southwest')
    saveas(gcf,['2-D prefactor Kramers, D_0 delta t=',num2str(D_times_delta_t_matrix(j)),'.png'])

    %% 2-D Numeric's rate 2D plots
    figure(200+j);hold on;
    xlabel('\theta_0')
    ylabel('ratio')
    title(['2-D ratio sim/Numeric, D_0 \delta t=',num2str(D_times_delta_t_matrix(j))])
    axis([1 inf 0 5])
    Kramer_measured=plot(1,1,'-.g');
    Kramer=plot(1,1,'-.r');
    Kramer_eff=plot(1,1,':b');
    legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'},'Location','northwest')
    saveas(gcf,['2-D prefactor Numeric, D_0 delta t=',num2str(D_times_delta_t_matrix(j)),'.png'])
            end

    end
% end
end


%%


figure(97)
axis([-inf inf 1 3.5 0 10])
view([0 0 1 ])
% view([5 2 5])
% view([1 0 0 ])
t=linspace(1,100,100)
% plot3(t,1.2+0.3*log(t)/log(10),ones(size(t)))

%%
figure(97)
view([0 0 1])
% view([5 2 5])
t=linspace(1,100,100)
% plot3(t,1.3+0.4*log(t)/log(10),ones(size(t)))
plot3(t,3.2+1.7/(2-log10(4))*(log10(t)-2),ones(size(t)))
%% 
figure(4)
view([1 0 0])
% view([5 2 5])
%%
% figure(5)
% view([0 0 1])
% t=linspace(1,100,100)
% plot3(t,1.2+0.3*log(t)/log(10),ones(size(t)))

%% Fitting for the masked data

    switch part1 
        case 'on'
    figure(11)

    % fit3D(delta_t_fit2D, theta_fit2D, prefactor_fit2D)
    [fitresult, gof] = fit3D(D_delta_t_fit1D, theta_fit1D, prefactor_fit1D)
    xlabel('D_0 \delta t')
    ylabel('\theta_0')
    zlabel('prefactor')
    title(['1-D prefactor=',num2str(fitresult.p00+fitresult.p01),'+(',num2str(fitresult.p10),'D_0 \delta t)+(',num2str(fitresult.p01),'(\theta_0-1))'])
    view([5 2 5])
    saveas(gcf,['1-D prefactor fitted 3D.png'])
    saveas(gcf,['1-D prefactor fitted 3D.fig'])
    %% 
    figure(12)
    % fit2D(theta_ratio2D,Kramers_ratio2D)
    [fitresult, gof] = fit2D(theta_ratio1D,Kramers_ratio1D)
    xlabel('\theta_0')
    ylabel('prefactor')
    title(['1-D prefactor=',num2str(fitresult.p2+fitresult.p1),'+(',num2str(fitresult.p1),'(\theta_0-1))'])
    saveas(gcf,['1-D prefactor fitted 2D.png'])
    end
    
    switch part2
        case 'on'
    figure(8)

    [fitresult, gof] = fit3D(D_delta_t_fit2D, theta_fit2D, prefactor_fit2D);
    % fit3D(D_delta_t_fit1D, theta_fit1D, prefactor_fit1D)
    xlabel('D_0 \delta t')
    ylabel('\theta_0')
    zlabel('prefactor')
    title(['2-D prefactor=',num2str(fitresult.p00+fitresult.p01),'+(',num2str(fitresult.p10),'D_0 \delta t)+(',num2str(fitresult.p01),'(\theta_0-1))'])
    view([5 2 5])
    saveas(gcf,['2-D prefactor Kramer fitted 3D.png'])
    saveas(gcf,['2-D prefactor Kramer fitted 3D.fig'])
    %% 
    figure(9)
    [fitresult, gof] = fit2D(theta_ratio2D,Kramers_ratio2D)
    % fit2D(theta_ratio1D,Kramers_ratio1D)
    xlabel('\theta_0')
    ylabel('prefactor')
    title(['2-D prefactor=',num2str(fitresult.p2+fitresult.p1),'+(',num2str(fitresult.p1),'(\theta_0-1))'])
    saveas(gcf,['2-D prefactor Kramer fitted 2D.png'])
    end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sub Functions
%% Masking functions
function mask=pre_filter(a,b,err)
mask(1:length(a))=0;
for i=1:length(a)
    if abs(a(i)-b(i))/abs(a(i))<err
        mask(i)=1;
    else
        mask(i)=NaN;
    end
end
end

function cutoff_mask=cutoff(data,threshold)
cutoff_mask(1:length(data))=0;
for i=1:length(data)
    if data(i)>threshold
        cutoff_mask(i)=1;
    else
        cutoff_mask(i)=NaN;
    end
end
end

function mask_diag=diagonal2D(theta_0,delta_t,D_0)
mask_diag(1:length(theta_0))=0;
for i=1:length(theta_0)
    if theta_0(i)>1.2+0.3*log(D_0*delta_t)/log(10)
        mask_diag(i)=NaN;
    else
        mask_diag(i)=1;
    end
end
end

function mask_diag=diagonal1D(theta_0,delta_t,D_0)
mask_diag(1:length(theta_0))=0;
for i=1:length(theta_0)
%     plot3(t,3.2+1.7/(2-log10(4))*(log10(t)-2),ones(size(t)))
    %     if theta_0(i)>1.3+0.4*log(D_0*delta_t)/log(10)
    if theta_0(i)>3.2+1.7/(2-log10(4))*log10(D_0*delta_t)
        mask_diag(i)=NaN;
    else
        mask_diag(i)=1;
    end
end
end

%% Fitting 3D plot for 2-dim case
function [fitresult, gof] = fit3D(x, y, z)
[xData, yData, zData] = prepareSurfaceData( x, y, z );

% Set up fittype and options.
ft = fittype( 'poly11' );

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, [xData, yData], zData );
% legend( h, 'untitled fit 1', 'z_fit vs. x_fit, y_fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes

grid on
end
%% Fitting 2D plot for 2-dim case
function [fitresult, gof] = fit2D(x, y)


[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
% figure( 'Name', 'untitled fit 2' );
h = plot( fitresult, xData, yData );
% legend( h, 'Kramers_ratio vs. theta_ratio', 'untitled fit 2', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
% title(['prefactor=',num2str(fitresult.p00),'+',num2str(fitresult.p10),'D_0 \delta t+',num2str(fitresult.p01),'\theta_0'])


grid on
end


