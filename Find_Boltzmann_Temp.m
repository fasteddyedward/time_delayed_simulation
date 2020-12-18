function [fitresult, gof] = Find_Boltzmann_Temp(x, ln_p , omega_0, delta_t,k_B)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, ln_p );

% Set up fittype and options.
% -1./(k_B*T).*(omega_0*delta_t^3/24.*omega.^4-1/2*(omega_0*delta_t-1)*omega.^2)-lnZ




% ft = fittype( [num2str(delta_t),'/a*x.^4-c'], 'independent', 'x', 'dependent', 'y' );
ft = fittype( [num2str(-1./k_B),'./T.*(',num2str(omega_0*delta_t^3/24),'.*omega.^4-1/2*',num2str(omega_0*delta_t-1),'.*omega.^2)-lnZ'], 'independent', 'omega', 'dependent', 'ln_p' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.001 log(8)];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% figure(1)
h = plot( fitresult, xData, yData );
% h.Visible=plot_Boltz_fit;
legend( h, 'ln_p vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
% xlabel( 'x', 'Interpreter', 'none' );
% ylabel( 'ln_p', 'Interpreter', 'none' );
grid on
end