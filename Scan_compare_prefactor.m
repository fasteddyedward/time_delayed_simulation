    %% Comparing the coefficients for the prefactor with all different D_0's
    
    close all
    clear
D_0_matrix=[1 5 10 20]
%%
% D_0_matrix=[ 20]
for D_0_index=1:length(D_0_matrix)
    dt=0.01;
    D_0=D_0_matrix(D_0_index);
    delta_t_matrix_vik=[0.1 0.2 0.5 1 2 5 10]'  %% Make 1-D simulation
    delta_t_matrix=[0.1 0.2 0.5 1 2 5 10]' %% for 2-D simulation
    theta_0_1D_matrix=(logspace(log10(1),log10(3),30));
    file_name_pc=['2021.1.28_scan_pc_D_0=',num2str(D_0)];
    file_name_vik=['2021.1.28_scan_vik_D_0=',num2str(D_0)];
    
    %% Determining figure numbers for D_0*delta_t
    D_times_delta_t_matrix=reshape(D_0_matrix.*delta_t_matrix,[length(D_0_matrix)*length(delta_t_matrix),1]);
    D_times_delta_t_matrix=sort(unique(D_times_delta_t_matrix));
    

    %% 1-D 
        load([file_name_vik,'.mat'])
        figure(1);hold on;
        for i=1:length(delta_t_matrix_vik)
            plot3(ones(length(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:),'-.g')
            plot3(ones(length(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateA_plot(i,:),'--r')
            plot3(ones(length(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:),':b')
            
        end
        if D_0_index==length(D_0_matrix)
        v = [5 2 5];
        [caz,cel] = view(v);
        grid on;
        xlabel('\delta t')
        ylabel('\theta_0')
        axis([0 inf 1 inf 0 3])
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        title('1-D ratio')
%         sim=plot(1,1,'kx')
        numeric=plot(1,1,'g-.');
        Kramer=plot(1,1,'r-.');
        Kramer_eff=plot(1,1,'b:');
        legend([ numeric Kramer Kramer_eff],{'Numeric FPE, D_\theta=4D_0/(R^2\theta_0^2)', 'Kramers, D_\theta=4D_0/(R^2\theta_0^2))','Kramers, D_\theta=2D_0/R^2'},'Location','southwest')
        saveas(gcf,['1-D prefactor 3D.png'])
        end
    %% Plotting 1-D with one specific
%     load([file_name_vik,'.mat'])
%     figure(11);hold on;
    for i=1:length(delta_t_matrix_vik)
        j=find(D_0*delta_t_matrix_vik(i)==D_times_delta_t_matrix);
        figure(10+j);hold on;
        plot(AvN_plot,1./AratesV_plot(i,:)./rateN_plot(i,:),'-.g')
        plot(AvN_plot,1./AratesV_plot(i,:)./rateA_plot(i,:),'--r')
        plot(AvN_plot,1./AratesV_plot(i,:)./rateAeff_plot(i,:),':b')
%         if D_0_index==length(D_0_matrix)
%             xlabel('\theta_0')
%             ylabel('ratio')
%             title(['1-D prefactor, D_0 \delta t=',num2str(D_times_delta_t_matrix(j))])
%             axis([1 inf 0 5])
%             numeric=plot(1,1,'g-.');
%             Kramer=plot(1,1,'r-.');
%             Kramer_eff=plot(1,1,'b:');
%             legend([ numeric Kramer Kramer_eff],{'Numeric FPE, D_\theta=4D_0/(R^2\theta_0^2)', 'Kramers, D_\theta=4D_0/(R^2\theta_0^2))','Kramers, D_\theta=2D_0/R^2'},'Location','northwest')
%             saveas(gcf,['1-D prefactor, D_0 delta t=',num2str(D_times_delta_t_matrix(j)),'.png'])
%         end
    end

            
            
            
            
    %% Part 2
    
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

        %% 2-D Kramers rate 3D plot
        figure(2);hold on
        for i=1:length(delta_t_matrix) 
            plot3(ones(length(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(trans_rate_measured(i,:).*iF(-1+theta_0_matrix(i,:))),'-.g')
            plot3(ones(length(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer(i,:).*iF(-1+theta_0_matrix(i,:))),'-.r')
            plot3(ones(length(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(transition_rate_Kramer_eff(i,:).*iF(-1+theta_0_matrix(i,:))),':b')
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
        title('Kramer')
        axis([-inf inf 1.0 1.6 1/min(time_duration) 5]) %% The minimum non-zero transition rate is 1/Obs_time_steps.

        Kramer_measured=plot(1,1,'-.g');
        Kramer=plot(1,1,'-.r');
        Kramer_eff=plot(1,1,':b');
        legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
        saveas(gcf,['2-D prefactor Kramers 3D.png'])
        end
            %% 2-D Kramer's rate 2D plots
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
        %% 2-D Numerics rate
                figure(3);hold on
        for i=1:length(delta_t_matrix)
            plot3(ones(length(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_measured_2D(i,:).*iF(-1+theta_0_matrix(i,:))),'-.g')
            plot3(ones(length(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_2D(i,:).*iF(-1+theta_0_matrix(i,:))),'-.r')
            plot3(ones(length(theta_0_1D_matrix))*D_0*delta_t_matrix(i),theta_0_matrix(i,:),trans_rate_sim_2D(i,:).*iF(-1+theta_0_matrix(i,:))./(rateN_eff_2D(i,:).*iF(-1+theta_0_matrix(i,:))),':b')
            
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
        title('Numeric')
        axis([-inf inf 1.0 1.6 1/min(time_duration) 5]) %% The minimum non-zero transition rate is 1/Obs_time_steps.
        Kramer_measured=plot(1,1,'-.g');
        Kramer=plot(1,1,'-.r');
        Kramer_eff=plot(1,1,':b');
        legend([ Kramer_measured Kramer Kramer_eff],{'D_{eff}','4D_0/(R^2\theta_0^2)','2D_0/R^2'})
        saveas(gcf,['2-D prefactor Numeric 3D.png'])
        end
            %% 2-D Numeric's rate 2D plots
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

%% Legending and Saving (The commented codes above don't cover all j's, so we need this section of codes)
for j=1:length(D_times_delta_t_matrix)

    %% Plotting 1-D with one specific
    figure(10+j)
    xlabel('\theta_0')
    ylabel('ratio')
    title(['1-D prefactor, D_0 \delta t=',num2str(D_times_delta_t_matrix(j))])
    axis([1 inf 0 5])
    numeric=plot(1,1,'g-.');
    Kramer=plot(1,1,'r-.');
    Kramer_eff=plot(1,1,'b:');
    legend([ numeric Kramer Kramer_eff],{'Numeric FPE, D_\theta=4D_0/(R^2\theta_0^2)', 'Kramers, D_\theta=4D_0/(R^2\theta_0^2))','Kramers, D_\theta=2D_0/R^2'},'Location','northwest')
    saveas(gcf,['1-D prefactor, D_0 delta t=',num2str(D_times_delta_t_matrix(j)),'.png'])
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





