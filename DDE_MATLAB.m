clear;
close all
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
%   s = ones(3,1);
  s=normrnd(2,10,3,1);
end
%-------------------------------------------