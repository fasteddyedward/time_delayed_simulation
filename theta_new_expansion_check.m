%% 2021.1.22 This files simulate the dot(theta)=2*omega_0*cos( )*sin( ) with another language to check if the codes are okay.
%% Interestingly the codes work sometimes, sometimes not. When they work, they look very much alike the phi equation 
%% dot(phi)=omega_0*sin()
clear;

dt=10^-4;
steps=10^7;
ntau=0.02/dt;
omega_0=60;
D=1
phi(1:1+steps+ntau,1)=0;
theta_new(1:1+steps+ntau,1)=0;
theta_new(2+ntau:end)=NaN;
% theta_new=1;
% r=sqrt(2*D*dt)*normrnd(0,1,steps+ntau+1,1);
r=normrnd(0,sqrt(2*D*dt),steps+ntau+1,1);


%%
for k=1:steps
    phi(ntau+k+1)=phi(ntau+k)+omega_0*sin(phi(k+ntau)-phi(k))*dt+r(k);
%     theta_new(ntau+k+1)=theta_new(ntau+k)+2*omega_0*cos((theta_new(ntau+k)+theta_new(k))/2)*sin((theta_new(ntau+k)-theta_new(k))/2)*dt+r(ntau+k)-r(k);
    theta_new(ntau+k+1)=theta_new(ntau+k)+omega_0*(sin(theta_new(ntau+k))-sin(theta_new(k)))*dt+r(ntau+k)-r(k);
    
%     plot(theta_new)
%     if mod(k,1000)==0
%         clf;plot(theta_new)
%         k
%     end
    
end
theta_phi=phi(1+ntau:end)-phi(1:end-ntau);
%%
figure(1) ;clf
plot(theta_new)

figure(2);clf;
plot(theta_phi)

figure(3);clf
histogram(theta_new)

figure(4);clf;
histogram(theta_phi)

% figure(5);clf;
% histogram(diff(theta_new))

% figure(6);clf;
% histogram(diff(theta_phi))
figure(7);clf;
hold on;
histogram(theta_phi)
histogram(theta_new)
theta_0=omega_0*ntau*dt

