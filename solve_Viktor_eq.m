clear
close all
delta_t=1;
omega_0=1.1;
%% Choose which omega_plus: 1. For const theta model; 2. For the second order ODE approx.
omega_plus_1=sqrt(10-sqrt(1/42./(omega_0*delta_t).^6+120./(omega_0*delta_t)-20));
omega_plus_2=sqrt(6/(omega_0*delta_t^3)*(omega_0*delta_t-1));

%% Tolerance setup
opts = odeset('RelTol',1e-5);
%% First Order
[t,y] = ode45(@(t,y) vdp1(t,y,omega_0,delta_t),[0 100],[10],opts);
%% Second Order
% [t,y] = ode45(@(t,y) vdp1(t,y,omega_0,delta_t),[0 100],[1;-1]);
%% Third Order
% [t,y] = ode45(@(t,y) vdp1(t,y,omega_0,delta_t),[0 10],[abs(omega_0);0.01;0]);
%% Fourth Order
% [t,y] = ode45(@(t,y) vdp1(t,y,omega_0,delta_t),[0 10],[abs(omega_0);0.01;0;0]);
%%

figure(1)
plot(t(1:end),y(1:end,1),'-')
axis([-inf inf 0 inf])
%%
clf
plot(y(1:end,1))
title('Solution of Viktor''s Approximation with ODE45');
xlabel('Time t');
ylabel('Solution \omega');
% set(gca, 'YScale', 'log')
if imag(omega_plus_1)==0
    yline(omega_plus_1,'r')
else
    yline(real(omega_plus_1),'r')
end
if imag(omega_plus_2)==0
    yline(omega_plus_2,'g');
else
    yline(real(omega_plus_1),'g')
end
legend('Numeric solution','\omega_+ (full eqn.)','\omega_+ (approx. eqn.)')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')







%%
function dydt = vdp1(t,y,omega_0,delta_t)
%% First order ODE
omega_plus=sqrt(6/(omega_0*delta_t^3)*(omega_0*delta_t-1));
dydt(1,1)=0;
dydt(1)=delta_t/3*(omega_plus^2-y(1)^2)*y(1);
% dydt(2)=(y(1)*delta_t-y(1)/omega_0-delta_t^3*y(1)^3/6-delta_t^2*y(2)/2)/(-delta_t^3/6);
% dydt(2)=3/delta_t*y(2)-(omega_plus^2-y(1)^2)*y(1); %same as line above
% dydt(2)=(1/2*y(2)*delta_t^2-y(1)*delta_t+asin(y(1)/omega_0)/(delta_t^3/6));
% dydt(2)=(y(1)*delta_t-y(1)/omega_0-delta_t^3*y(1)^3/6-delta_t^2*y(2)/2)/(+delta_t^3/6);%  Reverted sign for mass



%% Second order ODE
% % omega_plus=sqrt(6/(omega_0*delta_t^3)*(omega_0*delta_t-1));
% dydt(1:2,1)=0;
% dydt(1)=y(2);
% % dydt(2)=(y(1)*delta_t-y(1)/omega_0-delta_t^3*y(1)^3/6-delta_t^2*y(2)/2)/(-delta_t^3/6);
% % dydt(2)=3/delta_t*y(2)-(omega_plus^2-y(1)^2)*y(1); %same as line above
% % dydt(2)=(1/2*y(2)*delta_t^2-y(1)*delta_t+asin(y(1)/omega_0)/(delta_t^3/6));
% dydt(2)=(y(1)*delta_t-y(1)/omega_0-delta_t^3*y(1)^3/6-delta_t^2*y(2)/2)/(+delta_t^3/6);%  Reverted sign for mass
% 
% 
% % x_plus=10;
% % dydt(1)=y(2);
% % dydt(2)=y(2)-(x_plus^2-y(1)^2)*y(1)
% % dydt(2)=y(2)+0.001*(x_plus^2-y(1)^2)*y(1);
% % omega_plus_sqr=(6/(omega_0*delta_t^3)*(omega_0*delta_t-1));
% % dydt(2)=-1/(omega_0*delta_t^3/6)*(-omega_0*delta_t^2/2*y(2)-omega_0*delta_t^3/6*(omega_plus_sqr-y(1)^2)*y(1));

%% Third Order ODE: original 
% dydt(1:3,1)=0;
% dydt(1)=y(2);
% dydt(2)=y(3);
% dydt(3)=y(3)^2/y(2);
%% Third Order ODE: second method
% dydt(1:3,1)=0;
% dydt(1)=y(2);
% dydt(2)=y(3);
% dydt(3)=24/delta_t^4*(-y(1)/omega_0+y(1)*delta_t-1/2*y(2)*delta_t^2+1/6*y(3)+1/4*y(1)^2*y(2)*delta_t^4-1/6*y(1)*delta_t^3);
%% Fourth Order ODE
% dydt(1:4,1)=0;
% dydt(1)=y(2);
% dydt(2)=y(3);
% dydt(3)=y(4);
% dydt(4)=(-1/2*y(3)^2+1/6*y(3)*y(4)*delta_t+1/2*y(2)*y(4)+y(1)*y(2)^3*delta_t)/(1/6*y(2)*delta_t);
% % dydt(4)=0;
end



