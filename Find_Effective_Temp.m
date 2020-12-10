clear;
close all;


%%
Date='2020.11.24'
nth_take=7
delta_t_matrix=2
T_matrix=[1]
v_0_matrix=[3.5:0.1:10]
dt=10^-2




%%
plot_hist_fit='no';
%%
V_0_matrix=v_0_matrix;
Delta_t_matrix=delta_t_matrix;
%%
sigma_matrix=[];
mu_matrix=[];
T_eff_matrix=[];
D_eff_matrix=[];
for delta_t_index=1:length(Delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(V_0_matrix)
            %             if  nth_take~=7
            if 1
                %             if ismember(nth_take,nth_interest)
close all
%% Input File Name
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(V_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];
load([movie_name,'.mat'])

%% Finding the sigma of the omega (effective temperature)
f=figure(1),clf;
switch plot_hist_fit
    case 'no'
        f.Visible='off';
end
h=histogram(diff(theta(2,:)));
y=h.Values;
x=h.BinEdges(1:end-1)+0.5*h.BinWidth;
% plot(x,y);
[fitresult,gof]=fit_Gaussian(x,y,plot_hist_fit);
% a1*exp(-((x-b1)/c1)^2)
mu=fitresult.b1;
sigma=fitresult.c1/sqrt(2);

%% Calculating the D_eff and T_eff
D_eff=sigma^2/(4*dt);
T_eff=(v_0/(2*a)*delta_t^2*sigma^2)/(8*k_B*dt);

save([movie_name,'.mat'],'D_eff','T_eff','-append')
%% Appending the matrices
D_eff_matrix=[D_eff_matrix D_eff];
T_eff_matrix=[T_eff_matrix T_eff];
sigma_matrix=[sigma_matrix sigma];
mu_matrix=[mu_matrix mu];
nth_take
            end
            nth_take=nth_take+1;
        end
        nth_take=nth_take+1;
    end
    nth_take=nth_take+1;
end
%% Transition Rates
omega_0_matrix=v_0_matrix/(2*a);
transition_rate_theory=sqrt(2)./(pi.*omega_0_matrix.*delta_t_matrix.^2).*(omega_0_matrix.*delta_t_matrix-1).*exp(-3/2*(omega_0_matrix.*delta_t_matrix-1).^2./(k_B.*T_eff_matrix.*delta_t_matrix.^3));

%% For delta_t_matrix
if length(delta_t_matrix)>1
figure(2)
plot(delta_t_matrix,sigma_matrix)
title(['Gaussian \sigma vs \delta t, v_0= ',num2str(v_0),', T=',num2str(T)])
% title('Gaussian \sigma vs \delta t')
xlabel('\delta t')
ylabel('\sigma')

% mu might be very small and negligible
figure(3)
plot(delta_t_matrix,mu_matrix)
title('Gaussian \mu vs \delta t')
xlabel('\delta t')
ylabel('\mu')

figure(4)
plot(delta_t_matrix.^2,v_0_matrix*delta_t_matrix.^2.*sigma_matrix.^2,'o')
title(['v_0*\delta t^2*\sigma^2 vs \delta t, v_0= ',num2str(v_0),', T=',num2str(T)])
% title('Gaussian \sigma vs v_0')
xlabel('\delta_t^2')
ylabel('v_0*\delta t^2*\sigma^2 ')

figure(5)
plot(delta_t_matrix,T_eff_matrix,'o')
title(['T_{eff} vs \delta t, v_0= ',num2str(v_0),', T=',num2str(T)])
xlabel('\delta t')
ylabel('T_{eff}')

figure(6)
plot(delta_t_matrix,D_eff_matrix,'o')
title(['D_{eff} vs \delta t, v_0= ',num2str(v_0),', T=',num2str(T)])
xlabel('\delta t')
ylabel('D_{eff}')

figure(7)
plot(delta_matrix,transition_rate_theory,'o')
end
%% For v_0_matrix
if length(v_0_matrix)>1
figure(2)
plot(v_0_matrix,sigma_matrix)
title(['Gaussian \sigma vs v_0 t, \delta t= ',num2str(delta_t),', T=',num2str(T)])
% title('Gaussian \sigma vs v_0')
xlabel('v_0')
ylabel('\sigma')

% mu might be very small and negligible
figure(3)
plot(v_0_matrix,mu_matrix)
title('Gaussian \mu vs v_0')
xlabel('v_0')
ylabel('\mu')

figure(4)
plot(v_0_matrix,v_0_matrix*delta_t_matrix.^2.*sigma_matrix.^2,'o')
title(['v_0*\delta t^2*\sigma^2 vs v_0 t, \delta t= ',num2str(delta_t),', T=',num2str(T)])
% title('Gaussian \sigma vs v_0')
xlabel('v_0')
ylabel('v_0*\delta t^2*\sigma^2 ')

figure(5)
plot(v_0_matrix,T_eff_matrix,'o')
title(['T_{eff} vs v_0, \delta t= ',num2str(delta_t),', T=',num2str(T)])
xlabel('v_0')
ylabel('T_{eff}')

figure(6)
plot(v_0_matrix,D_eff_matrix,'o')
title(['D_{eff} vs v_0, \delta t= ',num2str(delta_t),', T=',num2str(T)])
xlabel('v_0')
ylabel('D_{eff}')

figure(7)
plot(v_0_matrix,transition_rate_theory,'o')
plot(omega_0_matrix.*delta_t_matrix,transition_rate_theory,'o')
end



