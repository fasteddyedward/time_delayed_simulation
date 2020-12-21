Date='2020.12.16'
nth_take=1
delta_t_matrix=2
T_matrix=[1]
v_0_matrix=[3.5:0.1:7]
dt=10^-1
intrinsic_delay=0.0 % Intrinsic delay
Obs_time_steps=10^5


%%
% Nth_take=nth_take;
V_0_matrix=v_0_matrix;
Delta_t_matrix=delta_t_matrix;
%%

%%
auto_exp_matrix=[];
%%
for delta_t_index=1:length(Delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(V_0_matrix)
            %             if  nth_take~=7
            if 1
%% Input File Name
movie_name=[Date,',dt=',num2str(dt),' take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(V_0_matrix(v_0_index)),', delta_t=',num2str(Delta_t_matrix(delta_t_index))];
[movie_name,'.mat'];
load([movie_name,'.mat'])
%% Start Calculating Autocorrelation Functions
close all
num_lag=10000;
[acf,lags,bounds] =autocorr((theta(2,:)),'NumLags',num_lag);
time_lags=0:dt:(num_lag)*dt;
plot(time_lags,acf)
xline(delta_t,'g')
xlabel('\tau')
ylabel('g(\tau)')
title('Autocorrelation function of \theta')
legend('g(\tau)',['\delta t=',num2str(delta_t)])
%% Fitting the Autocorrelation Functions
[fitresult, gof] = fit_exp(time_lags, acf);
auto_exp_matrix=[auto_exp_matrix -fitresult.b];
            end
            nth_take=nth_take+1;
        end
        nth_take=nth_take+1;
    end
    nth_take=nth_take+1;
end

%% Plotting the exponenets  % These can be changed into omega_0*delta_t if we import R_matrix as well
figure
plot(v_0_matrix,1./auto_exp_matrix)
xlabel('v_0')
ylabel('Correlation exponent time scale ')

function [fitresult, gof] = fit_exp(time_lags, acf)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( time_lags, acf );

% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.019407441306517 -0.0373033083348603];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'acf vs. time_lags', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'time_lags', 'Interpreter', 'none' );
ylabel( 'acf', 'Interpreter', 'none' );
grid on
end


