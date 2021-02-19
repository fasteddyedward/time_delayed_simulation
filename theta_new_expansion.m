%% 2021.1.22 This files simulate the dot(theta)=2*omega_0*cos( )*sin( ) with language used in the previous simulations
clear;
% close all


%% This setup doesn't give a good match (the scale plot neither), but still D_theta=2D_0/R^2
%% Could it be that at 
% nth_take=1
% delta_t_matrix=0.2
% T_matrix=[1]
% v_0_matrix=[5:0.5:14]
% dt=10^-3  % Don't set to be 10^-1 or else the simulation breaks
% R_mean=1; % actually R doesn't matter, we just put this to fit the 2-D simulations.y
% Obs_time_steps=10^7


%% This setup gives very good match for full simulation, and also D_theta=2D_0/R^2
% nth_take=2
% delta_t_matrix=2
% T_matrix=[1]
% v_0_matrix=[5.2:0.2:8]
% dt=10^-1
% Obs_time_steps=10^7
% R_mean=10;

%% This setup fits incredibly well for full simulation
% nth_take=3 % for a=5, particle 1 fixed
% delta_t_matrix=[1.8:0.1:3.4]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-2
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^5
% R_mean=10;


%%
nth_take=4
delta_t_matrix=0.02
T_matrix=[1]
v_0_matrix=[5:0.5:14]
dt=10^-2  % Don't set to be 10^-1 or else the simulation breaks
R_mean=0.1; % actually R doesn't matter, we just put this to fit the 2-D simulations
Obs_time_steps=10^6
D_0=0.1
%% This setup starts from theta_0>1
% nth_take=5 % for a=5, particle 1 fixed
% delta_t_matrix=[2.5:0.1:3.4]
% T_matrix=[1]
% v_0_matrix=5
% dt=10^-2
% intrinsic_delay=0.0 % Intrinsic delay
% Obs_time_steps=10^7
% R_mean=10;


%%
plot_hist_fit_T_eff='yes';
%% Matrices for analysis

%%
for delta_t_index=1:length(delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(v_0_matrix)
            if 1
                delta_t=delta_t_matrix(delta_t_index);
                v_0= v_0_matrix(v_0_index);
%                 D_0=1;
                
                D=D_0/R_mean^2;
                omega_0=v_0/R_mean;
                omega_0*delta_t
                % omega_0=1;
                % delta_t=100;
                
                %% Simulation setup
                phi(1:1+Obs_time_steps+round(delta_t/dt))=0;
                theta_new=normrnd(0,sqrt(2*D*dt),[Obs_time_steps,1]);
%                 theta_new(1:1+Obs_time_steps+round(delta_t/dt))=0;
                
                theta_0=omega_0*delta_t;
                theta_plus_approx_theory=sqrt(6*(theta_0-1)/theta_0);
                r_theta_sim=normrnd(0,sqrt(8*D*dt/theta_0^2),[Obs_time_steps,1]);
                r_theta_new=normrnd(0,sqrt(2*D*dt),[Obs_time_steps+round(delta_t/dt),1]);
                r_phi=normrnd(0,sqrt(2*D*dt),[Obs_time_steps,1]);

                %% For full simulation (This is what I'm more interested in, deal with this first)
                for k=1:Obs_time_steps
                    phi(1+k+round(delta_t/dt))=phi(k+round(delta_t/dt))+omega_0*sin(phi(k+round(delta_t/dt))-phi(k))*dt+r_phi(k);
                    
                end
                theta=phi(1+round(delta_t/dt):1+Obs_time_steps+round(delta_t/dt))-phi(1:1+Obs_time_steps);
                %% For New simulio (This is what I'm more interested in, deal with this first)
                for k=1:Obs_time_steps
                    theta_new(1+k+round(delta_t/dt))=theta_new(k+round(delta_t/dt))+2*omega_0*cos((theta_new(k+round(delta_t/dt))+theta_new(k))/2)*sin((theta_new(k+round(delta_t/dt))-theta_new(k))/2)*dt+r_theta_new(k+round(delta_t/dt))-r_theta_new(k);
%                 theta_new(1+k)=theta_new(k)+1/(theta_new(k)*theta_0)*(-(theta_0-1/2*theta_0*theta(k)^2-1)-sqrt((theta(k)-1/2*theta_0*theta(k)^2-1)^2-2*theta_0*theta(k)*r_theta_new(k)))*dt;
                end
                
                
                %% Testing for transition rates
                [V_parameters,~,~,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0)
                [~,~,theta_new_plus,theta_new_minus,num_transitions_theta_new]=find_theta_plus(theta_new,theta_0)
                %% testing for theta_new
%                 close all
                theta_new_avg=movmean(theta_new,round(delta_t/dt));
                figure(1);clf;
                plot(theta_new)
                figure(2);clf;
                histogram(theta_new)
%                 figure(3);clf;
%                 plot(theta_new_avg);
%                 figure(4);clf;
%                 histogram(theta_new_avg);
                figure(5);clf;
                plot(theta);
                figure(6);clf;
                histogram(theta)
                %                 histogram(diff(theta_new))
                 theta_0
                 
                1
            end
            nth_take=nth_take+1
        end
        nth_take=nth_take+1
    end
    nth_take=nth_take+1
end
%%


