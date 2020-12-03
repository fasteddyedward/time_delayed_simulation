clear;
close all
theta=@(x)sqrt(10-sqrt(1./(42.*x.^6)+120./x-20));
theta_pl=@(x)sqrt(6./x.*(x-1));
x_line=0:0.01:10;
figure
hold on
plot(x_line,theta(x_line),'b')
plot(x_line,theta_pl(x_line),'r')
title('Two theoretical lines')
ylabel('\theta')
xlabel('\omega_0 \delta t')
legend('$\theta=\sqrt{10-\sqrt{\frac{1}{42\omega_0 \delta t}+\frac{120}{\omega_0 \delta t }-20}}$','$\theta=\sqrt{\frac{6}{\omega_0 \delta t}(\omega_0 \delta t-1)}$','interpreter','latex')
plot(x_line,-theta(x_line),'b')
plot(x_line,-theta_pl(x_line),'r')
% legend('$\frac{1}{2}$','interpreter','latex')

%%
lags = [1 0.2];
tspan = [0 5];
sol = dde23(@ddefun, lags, @history, tspan);
plot(sol.x,sol.y,'-o')
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2','y_3','Location','NorthWest');
function dydt = ddefun(t,y,Z) % equation being solved
  ylag1 = Z(:,1);
  ylag2 = Z(:,2);

  dydt = [ylag1(1); 
          ylag1(1)+ylag2(2); 
          y(2)];
end
%-------------------------------------------
function s = history(t) % history function for t <= 0
  s = ones(3,1);
end
%-------------------------------------------




% 
% %%
% m=[1 2 3 4;5 6 7 8;9 10 11 12]
% fcol=@(x)deal(x(:,1:2),x(:,3:4))
% [a b]=fcol(m)
% 
% deal(1)
% %%
% 
% x = 1:50;
% y = rand(1,50);
% s = 5; %Marker width in units of X
% h = scatter(x,y); % Create a scatter plot and return a handle to the 'hggroup' object
% %Obtain the axes size (in axpos) in Points
% currentunits = get(gca,'Units');
% set(gca, 'Units', 'Points');
% axpos = get(gca,'Position');
% set(gca, 'Units', currentunits);
% markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
% set(h, 'SizeData', markerWidth^2)
% %%
% ax = gca;
% AR = get(gca, 'dataaspectratio');
% if ~isequal(AR(1:2), [1 1])
%   error('Units are not equal on X and Y, cannot create marker size that is one unit on both');
% end
% oldunits = get(ax, 'Units');
% set(ax, 'Units', 'points');
% pos = get(ax, 'Position');    %[X Y width height]
% set(ax, 'Units', oldunits');
% XL = xlim(ax);
% points_per_unit = pos(3) / (XL(2) - XL(1));
% marker_size = points_per_unit .^2 * pi / 4;    %pi*r^2 but remember points_per_unit is d not r
% %%
% currentunits=ax.Units;
% ax.Position;
% ax.Units='Points';
% axpos=ax.Position;
% pba=ax.PlotBoxAspectRatio;
% ax.Units=currentunits;
% scale=min([axpos(3),pba(1)/pba(2)*axpos(4)])/diff(xlim);
% a=markerwidth.*scale;
% set(h,'SizeData',a.^2);
% %%
% x = linspace(0,10);
% y = sin(4*x);
% scatter(x,y,'o')
% %%
% ax = gca; % current axes
% ax.FontSize = 12;
% ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% ax.YLim = [-2 2];
% fig = gcf;
% ax = fig.CurrentAxes;
% ax.Units=1
% %%
% figure;    
% close all
%     h=scatter(0.5,0.5);
%     axis square
%     axis([0 1 0 1])
%     xlabel('x (\mum)')
%     ylabel('y (\mum)')
%     title 'cleaned data'
% %     set(gcf, 'Color', 'w');
% %     set(gcf, 'PaperPositionMode', 'auto');
%     set( findobj(gca,'type','line'), 'LineWidth', 2)
%     s = 0.4; %Marker width in units of X
%     currentunits = get(gca,'Units');
%     set(gca, 'Units', 'Points');
%     axpos = get(gca,'Position');
%     set(gca, 'Units', currentunits);
%     markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points
%     set(h, 'Marker', 'square')
%     set(h, 'SizeData', markerWidth^2)
%     hold on
%     h1=scatter(0.5,0.5);
%     s = 0.4;%Marker width in units of X
%     markerWidth = s/diff(xlim)*axpos(4); % Calculate Marker width in points
%     set(h1, 'Marker', 'o')
%     set(h1, 'SizeData', (markerWidth)^2)
%     %axpos =   54.6000   34.6500  325.5000  256.7250
%     ax=get(gcf);
% %     ax.Position(1)=1000
%     ax.Position
    
