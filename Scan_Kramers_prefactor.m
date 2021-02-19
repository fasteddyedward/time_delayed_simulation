%% This is basically similar to Scan_compare_prefactor, with the difference that now we use the potential U(theta) = -k_B T ln(p(theta)). 
%% scan over omega_0 and delta_t.
clear
close all

D_0_matrix=[ 5 10 20 50 100 200]
% D_0_matrix=[10 100]
for D_0_index=1:length(D_0_matrix)

    dt=0.01;
    D_0=D_0_matrix(D_0_index);
    % D_0=10 is the best
    % D_0=1 the rate is too low
    
    intrinsic_delay=0.
    
    % delta_t_matrix=1
    % theta_0_1D_matrix=1.2
    file_name_vik=['2021.2.18_scan_vik_D_0=',num2str(D_0)];
    
    runs=1
    delta_t_matrix_vik=[0.1 0.2 0.5 1 2 5 10]
    
    theta_0_1D_matrix=linspace(1.1,5,30);
    % theta_0=linspace(1.1,3,5)
    %% Simulating data
%     Scan_compare_Viktor(D_0,dt,50000,theta_0_1D_matrix,delta_t_matrix_vik,runs,file_name_vik)
    %%
    
    load([file_name_vik,'.mat'])
    
    %% Plotting the figures
    % close all;
    for i=1:length(delta_t_matrix_vik)
        figure(10);hold on
        plot(AvN_plot,1./AratesV_plot(i,:),'-o')
        plot(AvN_plot,rateA_potential_plot(i,:),'-o')
        plot(AvN_plot,rateN_potential_plot(i,:),'-o')
        % plot(Bins,ln_p)
        legend('simulation','Kramers','Numerical FPE');
        set(gca, 'YScale', 'log')
        xlabel('\theta_0')
        
        
        figure(11);hold on
        plot(AvN_plot,1./AratesV_plot(i,:)./rateA_potential_plot(i,:),'-o')
        set(gca, 'YScale', 'log')
        title('sim data/Kramers')
        axis([-inf inf 0 inf ])
        xlabel('\theta_0')
        
        figure(12);hold on
        plot(AvN_plot,rateN_potential_plot(i,:)./rateA_potential_plot(i,:),'-o')
        set(gca, 'YScale', 'linear')
        title('Numerical/Kramers')
        axis([-inf inf 0 inf ])
        xlabel('\theta_0')
        
        figure(13);hold on
        plot(AvN_plot,1./AratesV_plot(i,:)./rateN_potential_plot(i,:),'-o')
        set(gca, 'YScale', 'linear')
        title('sim data/Numerical')
        axis([-inf inf 0 inf ])
        xlabel('\theta_0')
        %% 3D plots w.r.t to delta_t
        figure(14);hold on
        plot3(ones(length(AvN_plot))*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateA_potential_plot(i,:),'k-o')
        set(gca, 'ZScale', 'log')
        xlabel('\delta t')
        ylabel('\theta_0')
        zlabel('sim data/Kramers')
        view([5 2 5])
        grid on
        
        
        %% 3D plots w.r.t to D_0 * delta_t
        figure(15);hold on;
        plot3(ones(size(AvN_plot))*D_0*delta_t_matrix_vik(i),AvN_plot,1./AratesV_plot(i,:)./rateA_potential_plot(i,:),'k-o')
        set(gca, 'ZScale', 'log')
        set(gca, 'XScale', 'log')
        xlabel('D_0 \delta t')
        ylabel('\theta_0')
        zlabel('sim data/Kramers')
        grid on
        axis([1 10^3  -inf inf -inf inf ])
        view([5 5 5])
        saveas(gcf, 'Prefactor_corrected.png')
        view([1 0 0 ])
        saveas(gcf, 'Prefactor_corrected_sideview.png')
      
    end
end






