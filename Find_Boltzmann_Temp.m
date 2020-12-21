function [fitresult, gof] = Find_Boltzmann_Temp(x, ln_p , omega_0, delta_t,k_B,omega_plus_type)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, ln_p );

% Set up fittype and options.
switch omega_plus_type
    %% Type omega_plus from approximate theory
    case 'approximate'
               omega_plus=sqrt(6/(omega_0*delta_t^3)*(omega_0*delta_t-1));
        %% Type omega_plus from full sine expansion
    case 'full_sine'
        omega_plus=sqrt(10-sqrt(1/42/(omega_0*delta_t)^6+120/(omega_0*delta_t)-20))/delta_t;
end

% -1./(k_B*T).*(omega_0*delta_t^3/24.*omega.^4-1/2*(omega_0*delta_t-1)*omega.^2)-lnZ
%         ft = fittype( [num2str(-1./k_B),'./T.*(',num2str(omega_0*delta_t^3/24),'.*omega.^4-1/2*',num2str(omega_0*delta_t-1),'.*omega.^2)-lnZ'], 'independent', 'omega', 'dependent', 'ln_p' );

% -1./(k_B*T).*omega_0*delta_t^3/24*(omega.^2-2*omega_plus^2).*omega.^2-lnZ
ft = fittype( [num2str(-1./k_B),'./T.*(',num2str(omega_0*delta_t^3/24),'.*(omega.^2-2*',num2str(omega_plus^2),').*omega.^2)-lnZ'], 'independent', 'omega', 'dependent', 'ln_p' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.001 log(8)];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


h = plot( fitresult, xData, yData );
legend( h, 'ln_p vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
grid on
end