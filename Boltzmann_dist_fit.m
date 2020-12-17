close all
load([movie_name,'.mat'],'theta')
bin_limit=1;
bin_interval=2*bin_limit/num_bins;
bin_loc=-bin_limit:bin_interval:bin_limit;
omega=theta(2,:)/delta_t;
h=histogram(omega,bin_loc);
y=h.Values;
x=h.BinEdges(1:end-1)+0.5*h.BinWidth;
dx=h.BinWidth;
xline(theta_plus/delta_t)
xline(theta_minus/delta_t)

%% Plotting the Prob. distribution of sim and theory. Note that p(omega) is prob. density
f=figure(1);clf;hold on;f.Visible='off';
title('p(\omega) of simulation')
p=y/length(omega); % Probability of the statistics
plot(x,p)
xline(theta_plus/delta_t)
xline(theta_minus/delta_t)


f=figure(2);clf;hold on;f.Visible='off';
title('p(\omega) of theory')
omega_0=v_0/R_mean(2);
omega_plus=sqrt(6/(omega_0*delta_t^3)*(omega_0*delta_t-1));
U=@(omega) omega_0*delta_t^3/24*(omega.^2-2*omega_plus.^2).*omega.^2;
p_omega=@(omega)exp(-U(omega)/(k_B*T_eff));
x_line=-bin_limit:0.01:bin_limit;
Z=integral(p_omega,-inf,inf) ;
plot(x_line,p_omega(x_line)/Z*dx) % dx has to be multiplied because p is prob. density
xline(+omega_plus)
xline(-omega_plus)

%% Comparing in terms of Probability
x_line=-bin_limit:0.01:bin_limit;
figure(3); clf; hold on
title('Comparing p(\omega) of Simulation and Theory')
plot(x_line,p_omega(x_line)/Z*dx)
plot(x,p)
legend('sim','theory')

%% Comparing in terms of Potential
% T_eff=1
figure(4);clf; hold on
title('Comparing U(\omega) of Simulation and Theory')
plot(x,-k_B*T_eff*(log(p)-log(p(end/2)))) % Setting the zero point to match, p(end/2) is p(x=0)
x_line=-1:0.01:1;
plot(x_line,U(x_line))
legend('sim','theory')

%%
ln_p=log(p/dx);
[fitresult, gof] =createFit(x,ln_p);

T_fit=1/(k_B*fitresult.a)*omega_0*delta_t^3/24
omega_plus_fit=sqrt(fitresult.b/2)
Z_fit=exp(-fitresult.c)
%% Error
(T_fit-T_eff)/T_eff
(omega_plus_fit-omega_plus)/omega_plus
(Z_fit-Z)/Z
%%
function [fitresult, gof] = createFit(x, ln_p)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, ln_p );

% Set up fittype and options.
ft = fittype( '-a*(x.^2-b).*x.^2+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.754686681982361 0.276025076998578 0.679702676853675];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'ln_p vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'ln_p', 'Interpreter', 'none' );
grid on
end