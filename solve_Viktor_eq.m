clear
close all
delta_t=10^-7;
omega_0=1;

omega_plus=sqrt(6/(omega_0*delta_t^3)*(omega_0*delta_t-1));
[t,y] = ode45(@(t,y) vdp1(t,y,omega_0,delta_t),[0 2],[1;-1]);
figure(1)
% plot(t,y(:,1),'-')
plot(t(1:end),y(1:end,1),'-')
%%
clf
plot(y(1:end-150,1))
title('Solution of van der Pol Equation (\mu = 1) with ODE45');
xlabel('Time t');
ylabel('Solution y');
% set(gca, 'YScale', 'log')

function dydt = vdp1(t,y,omega_0,delta_t)


dydt(1:2,1)=0;
dydt(1)=y(2);
% dydt(2)=(y(1)*delta_t-y(1)/omega_0-delta_t^3*y(1)^3/6-delta_t^2*y(2)/2)/(-delta_t^3/6);
dydt(2)=(1/2*y(2)*delta_t^2-y(1)*delta_t+asin(y(1)/omega_0)/(delta_t^3/6));



% x_plus=10;
% dydt(1)=y(2);
% dydt(2)=y(2)-(x_plus^2-y(1)^2)*y(1)
% dydt(2)=y(2)+0.001*(x_plus^2-y(1)^2)*y(1);
% omega_plus_sqr=(6/(omega_0*delta_t^3)*(omega_0*delta_t-1));
% dydt(2)=-1/(omega_0*delta_t^3/6)*(-omega_0*delta_t^2/2*y(2)-omega_0*delta_t^3/6*(omega_plus_sqr-y(1)^2)*y(1));

end



