clear;
Date='2020.12.16'
nth_take=1
delta_t_matrix=2
T_matrix=[1]
v_0_matrix=[3.5:0.1:7]
dt=10^-1
intrinsic_delay=0.0 % Intrinsic delay
Obs_time_steps=10^5
omega_plus_type='full_sine' % 'approximate'
%%
if exist('V_0_matrix','var')==0
    V_0_matrix=v_0_matrix;
end
if exist('Delta_t_matrix','var')==0
    Delta_t_matrix=delta_t_matrix;
end
%% Start Analyzing
for delta_t_index=1:length(Delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(V_0_matrix)
            %% Input File Name
% movie_name=['test3']
% movie_name=['2020.11.25,dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(V_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];
% movie_name=[Date,',dt=10e-3 take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];

            %% just_evaluate_this

            %             if  nth_take~=7
            if 1

close all


[movie_name,'.mat'];
load([movie_name,'.mat'])
load([movie_name,'.mat'],'theta')
num_bins=100;
bin_limit=1;
bin_interval=2*bin_limit/num_bins;
bin_loc=-bin_limit:bin_interval:bin_limit;
omega=theta(2,:)/delta_t;
h=histogram(omega,bin_loc);
% h.Visible='off'
omega_count=h.Values;
omega_bin=h.BinEdges(1:end-1)+0.5*h.BinWidth;
dx=h.BinWidth;
xline(theta_plus/delta_t)
xline(theta_minus/delta_t)

%% Plotting the Prob. distribution of sim and theory. Note that p(omega) is prob. density
f=figure(1);clf;hold on;f.Visible='off';
title('p(\omega) of simulation')
p=omega_count/length(omega); % Probability of the statistics
plot(omega_bin,p)
xline(theta_plus/delta_t)
xline(theta_minus/delta_t)


f=figure(2);clf;hold on;f.Visible='off';
title('p(\omega) of theory')
omega_0=v_0/R_mean(2);
    %% Omega_plus for approximate theory
omega_plus=sqrt(6/(omega_0*delta_t^3)*(omega_0*delta_t-1));
    %% Omega_plus_exact for expansion to 7th power in the original sin
% omega_plus_exact=sqrt(10-sqrt(1/42/(omega_0*delta_t)^6+120/(omega_0*delta_t)-20))/delta_t;
U=@(omega) omega_0*delta_t^3/24*(omega.^2-2*omega_plus.^2).*omega.^2;
p_omega=@(omega)exp(-U(omega)/(k_B*T_eff));
x_line=-bin_limit:0.01:bin_limit;
Z=integral(p_omega,-inf,inf) ;
plot(x_line,p_omega(x_line)/Z*dx) % dx has to be multiplied because p is prob. density
if exist('omega_plus','var')
xline(+real(omega_plus))
xline(-real(omega_plus))
end

%% Comparing in terms of Probability
% x_line=-bin_limit:0.01:bin_limit;
figure(3); clf; hold on
title('Comparing p(\omega) of Simulation and Theory')
plot(x_line,p_omega(x_line)/Z*dx)
plot(omega_bin,p)
legend('sim','theory')
xlabel('\omega')
ylabel('p(\omega)')

%% Comparing in terms of Potential
% T_eff=1
figure(4);clf; hold on
title('Comparing U(\omega) of Simulation and Theory')
plot(omega_bin,-k_B*T_eff*(log(p)-log(p(end/2))),'o') % Setting the zero point to match, p(end/2) is p(x=0)
x_line=-1:0.01:1;
plot(x_line,U(x_line))
legend('sim','theory')
xlabel('\omega')
ylabel('U(\omega)')
%% Fitting U=k_B*T*ln(p)+ln(Z)
ln_p=log(p/dx);
figure(5);clf;hold on;
[fitresult, gof] =Find_Boltzmann_Temp(omega_bin,ln_p,omega_0,delta_t,k_B,omega_plus_type);
T_Boltz=fitresult.T
Z_fit=exp(fitresult.lnZ)
xlabel('\omega')
ylabel('ln(p (\omega))')
title(['\omega_0 \deltat=',num2str(omega_0*delta_t),', T_{Boltzmann}/T_{eff}=',num2str(T_Boltz/T_eff),', sum(p(\omega))/Z_{fit}=',num2str(sum(exp(-U(bin_loc)/(k_B*T_Boltz))*dx)/Z_fit)])
% Z_fit is because this is fitted with the fitting function, and it is
% possible that the Z doens't give a completely normalized p(omega)
%% Fitting U=k_B*T*ln(p)+ln(Z) with only a small range
ln_p=log(p/dx);
figure(6);clf;hold on;
% interest_flag=(ln_p>ln_p(end/2));
threshold=0.5
interest_flag=(abs(omega_bin)<threshold)

% allowed_half_width=0.2
% interest_flag=(omega_plus-allowed_half_width<abs(omega_bin) & abs(omega_bin)<omega_plus+allowed_half_width)

[fitresult, gof] =Find_Boltzmann_Temp(omega_bin(interest_flag),ln_p(interest_flag),omega_0,delta_t,k_B,omega_plus_type);
plot(omega_bin(end/2),ln_p(end/2),'og')
T_Boltz=fitresult.T
Z_fit=exp(fitresult.lnZ)
xlabel('\omega')
ylabel('ln(p (\omega))')
title(['\omega_0 \deltat=',num2str(omega_0*delta_t),', T_{Boltzmann}/T_{eff}=',num2str(T_Boltz/T_eff),', sum(p(\omega))/Z_{fit}=',num2str(sum(exp(-U(bin_loc)/(k_B*T_Boltz))*dx)/Z_fit)])
% Z_fit is because this is fitted with the fitting function, and it is
% possible that the Z doens't give a completely normalized p(omega)
    %% extending fitted line
    figure(7);clf;hold on;
    plot(omega_bin,ln_p,'.')
    plot(x_line,-U(x_line)/(k_B*T_Boltz)-fitresult.lnZ);
    xline(real(omega_plus))
    xlabel('\omega')
    ylabel('ln(p (\omega))')
    title(['\omega_0 \deltat=',num2str(omega_0*delta_t),', T_{Boltzmann}/T_{eff}=',num2str(T_Boltz/T_eff),', sum(p(\omega))/Z_{fit}=',num2str(sum(exp(-U(bin_loc)/(k_B*T_Boltz))*dx)/Z_fit)])
    
%% 2020.12.18 Is Z okay for dist. p?
%     sum(exp(-U(bin_loc)/(k_B*T_fit))*dx)/Z_fit
            end
            nth_take=nth_take+1;
        end
        nth_take=nth_take+1;
    end
    nth_take=nth_take+1;
end
            
