function [rateN,rateN_eff,rateN_measured,rateA,rateAeff]=numeric_FPE(tau_matrix,theta_0_matrix,D_theta_matrix,D_0,R_matrix)
%% rateN: numeric FPE with D_theta=4*D_0/(R^2 theta_0^2)
%% rateN_eff: numeric FPE with D_theta=2*D_0/R^2
%% rateN_measured: numeric FPE with D_theta= sigma_theta^2/(2dt) (measured with hist(diff(theta))
%% rateA     = transition_rate_Krame:       Kramer's escape rate with D_theta=4*D_0/(R^2 theta_0^2)
%% rateA_eff = transition_rate_Kramer_eff:  Kramer's escape rate with D_theta=2*D_0/R^2


Av=theta_0_matrix; %% just subsitiute this guy;
% p2=0.005
% p1=100

%% This is for the original numerical FPE with D_theta = 4*D_0/(R^2 theta_0^2)
% Dx = 4*p2/p1^2;
for itau=1:length(tau_matrix)
    for ia = 1:1:size(theta_0_matrix,2)
        
        % parameters
        tau=tau_matrix(itau);
        %         p1 = Av(ia);
        p1=Av(itau,ia); % theta_0
        p2=D_0/R_matrix(itau,ia)^2; % D=D_0/R^2

        % potential
        thp2 = 6*(p1-1)/p1;
        Dx = 4*p2/p1^2; % D_theta = 4*D/theta_0^2 = 4*D_0/(R^2 theta_0^2) as in the general Kramer's formula
        
        % rate from FPE
        
        a = -1*sqrt(thp2)/2;
        N = 500;
        b = 3;
        xV = linspace(a,b,N + 1);
        Delta = xV(2) - xV(1);
        iX = zeros(1, (length(xV) - 2)*3 + 4);
        iY = zeros(1, (length(xV) - 2)*3 + 4);
        iL = zeros(1, (length(xV) - 2)*3 + 4);
        
        U = @(x) (x^3 - 2*x*thp2)*x/(12*tau);
        
        % rate matrix
        nux = Dx/Delta^2;
        Bx = 1/Dx;
        
        % left boundary (left jump allowed)
        x = a;
        Ui = U(x);
        UiL = U(x - Delta);
        UiR = U(x + Delta);
        rL = nux*exp(- Bx*(UiL - Ui)/2);
        rR = nux*exp(- Bx*(UiR - Ui)/2);
        iX(1) = 1;
        iY(1) = 1;
        iL(1) = -(rL + rR);
        rR = nux*exp(Bx*(UiR - Ui)/2);
        iX(2) = 1;
        iY(2) = 2;
        iL(2) = rR;
        
        % right boundary (right jump forbidden)
        x = xV(N+1);
        Ui = U(x);
        UiL = U(x - Delta);
        rL = nux*exp(- Bx*(UiL - Ui)/2);
        iX(length(iX)) = N + 1;
        iY(length(iX)) = N + 1;
        iL(length(iX)) = - rL;
        rL = nux*exp(Bx*(UiL - Ui)/2);
        iX(length(iX) - 1) = N + 1;
        iY(length(iX) - 1) = N;
        iL(length(iX) - 1) = rL;
        
        % all the rest
        for ix = 2:1:N
            middleC = 1 + (ix - 1)*3;
            x = xV(ix);
            Ui = U(x);
            UiL = U(x - Delta);
            UiR = U(x + Delta);
            rL = nux*exp(- Bx*(UiL - Ui)/2);
            rR = nux*exp(- Bx*(UiR - Ui)/2);
            iX(middleC) = ix;
            iY(middleC) = ix;
            iL(middleC) = -(rL + rR);
            rL = nux*exp(Bx*(UiL - Ui)/2);
            rR = nux*exp(Bx*(UiR - Ui)/2);
            iX(middleC - 1) = ix;
            iY(middleC - 1) = ix - 1;
            iL(middleC - 1) = rL;
            iX(middleC + 1) = ix;
            iY(middleC + 1) = ix + 1;
            iL(middleC + 1) = rR;
        end
        
        LM = sparse(iX,iY,iL);
        [V,e]=eigs(LM,1,0);
        jstatFPSS = e;
        
        rateAnal = sqrt(2)/pi*abs(p1-1)/(tau*p1)*exp(-3/4*(p1-1)^2/(p2*tau));
        rateAnalEff = sqrt(2)/pi*abs(p1-1)/(tau*p1)*exp(-3*(p1-1)^2./(2*p2*p1^2*tau)); %% in p.c.'s 2-D simulation, the D_eff is about 2*D/R^2=2*p2
        
        rateN(itau,ia) = - jstatFPSS;
        rateA(itau,ia) = rateAnal;
        rateAeff(itau,ia) = rateAnalEff;
        
    end
end

%% This is for FPE with D_theta = 2*D_0/R^2
% Dx = 2*p2;

for itau=1:length(tau_matrix)
    for ia = 1:1:size(theta_0_matrix,2)
        
        % parameters
        tau=tau_matrix(itau);
        p1=Av(itau,ia);
        p2=D_0/R_matrix(itau,ia)^2; % D_theta = 2*D = 2*D_0/R^2 as in the general Kramer's formula

        % potential
        thp2 = 6*(p1-1)/p1;
%         Dx = 4*p2/p1^2;
        Dx = 2*p2;
        % rate from FPE
        
        a = -1*sqrt(thp2)/2;
        N = 500;
        b = 3;
        xV = linspace(a,b,N + 1);
        Delta = xV(2) - xV(1);
        iX = zeros(1, (length(xV) - 2)*3 + 4);
        iY = zeros(1, (length(xV) - 2)*3 + 4);
        iL = zeros(1, (length(xV) - 2)*3 + 4);
        
        U = @(x) (x^3 - 2*x*thp2)*x/(12*tau);
        
        % rate matrix
        nux = Dx/Delta^2;
        Bx = 1/Dx;
        
        % left boundary (left jump allowed)
        x = a;
        Ui = U(x);
        UiL = U(x - Delta);
        UiR = U(x + Delta);
        rL = nux*exp(- Bx*(UiL - Ui)/2);
        rR = nux*exp(- Bx*(UiR - Ui)/2);
        iX(1) = 1;
        iY(1) = 1;
        iL(1) = -(rL + rR);
        rR = nux*exp(Bx*(UiR - Ui)/2);
        iX(2) = 1;
        iY(2) = 2;
        iL(2) = rR;
        
        % right boundary (right jump forbidden)
        x = xV(N+1);
        Ui = U(x);
        UiL = U(x - Delta);
        rL = nux*exp(- Bx*(UiL - Ui)/2);
        iX(length(iX)) = N + 1;
        iY(length(iX)) = N + 1;
        iL(length(iX)) = - rL;
        rL = nux*exp(Bx*(UiL - Ui)/2);
        iX(length(iX) - 1) = N + 1;
        iY(length(iX) - 1) = N;
        iL(length(iX) - 1) = rL;
        
        % all the rest
        for ix = 2:1:N
            middleC = 1 + (ix - 1)*3;
            x = xV(ix);
            Ui = U(x);
            UiL = U(x - Delta);
            UiR = U(x + Delta);
            rL = nux*exp(- Bx*(UiL - Ui)/2);
            rR = nux*exp(- Bx*(UiR - Ui)/2);
            iX(middleC) = ix;
            iY(middleC) = ix;
            iL(middleC) = -(rL + rR);
            rL = nux*exp(Bx*(UiL - Ui)/2);
            rR = nux*exp(Bx*(UiR - Ui)/2);
            iX(middleC - 1) = ix;
            iY(middleC - 1) = ix - 1;
            iL(middleC - 1) = rL;
            iX(middleC + 1) = ix;
            iY(middleC + 1) = ix + 1;
            iL(middleC + 1) = rR;
        end
        
        LM = sparse(iX,iY,iL);
        [V,e]=eigs(LM,1,0);
        jstatFPSS = e;
        
%         rateAnal = sqrt(2)/pi*abs(p1-1)/(tau*p1)*exp(-3/4*(p1-1)^2/(p2*tau));
%         rateAnalEff = sqrt(2)/pi*abs(p1-1)/(tau*p1)*exp(-3*(p1-1)^2./(2*p2*p1^2*tau)); %% in p.c.'s 2-D simulation, the D_eff is about 2*D/R^2=2*p2
        
        rateN_eff(itau,ia) = - jstatFPSS;
%         rateA(itau,ia) = rateAnal;
%         rateAeff(itau,ia) = rateAnalEff;
        
    end
end

%% THis is for FPE with D_theta = D_eff
for itau=1:length(tau_matrix)
    for ia = 1:1:size(theta_0_matrix,2)
        
        % parameters
        tau=tau_matrix(itau);
        p1=Av(itau,ia);
%         p2=D_theta_matrix(itau,ia);

        % potential
        thp2 = 6*(p1-1)/p1;

%         Dx = p2; 
Dx = D_theta_matrix(itau,ia); % D_theta = D_eff, measured with the sigma_theta^2 
        % rate from FPE
        
        a = -1*sqrt(thp2)/2;
        N = 500;
        b = 3;
        xV = linspace(a,b,N + 1);
        Delta = xV(2) - xV(1);
        iX = zeros(1, (length(xV) - 2)*3 + 4);
        iY = zeros(1, (length(xV) - 2)*3 + 4);
        iL = zeros(1, (length(xV) - 2)*3 + 4);
        
        U = @(x) (x^3 - 2*x*thp2)*x/(12*tau);
        
        % rate matrix
        nux = Dx/Delta^2;
        Bx = 1/Dx;
        
        % left boundary (left jump allowed)
        x = a;
        Ui = U(x);
        UiL = U(x - Delta);
        UiR = U(x + Delta);
        rL = nux*exp(- Bx*(UiL - Ui)/2);
        rR = nux*exp(- Bx*(UiR - Ui)/2);
        iX(1) = 1;
        iY(1) = 1;
        iL(1) = -(rL + rR);
        rR = nux*exp(Bx*(UiR - Ui)/2);
        iX(2) = 1;
        iY(2) = 2;
        iL(2) = rR;
        
        % right boundary (right jump forbidden)
        x = xV(N+1);
        Ui = U(x);
        UiL = U(x - Delta);
        rL = nux*exp(- Bx*(UiL - Ui)/2);
        iX(length(iX)) = N + 1;
        iY(length(iX)) = N + 1;
        iL(length(iX)) = - rL;
        rL = nux*exp(Bx*(UiL - Ui)/2);
        iX(length(iX) - 1) = N + 1;
        iY(length(iX) - 1) = N;
        iL(length(iX) - 1) = rL;
        
        % all the rest
        for ix = 2:1:N
            middleC = 1 + (ix - 1)*3;
            x = xV(ix);
            Ui = U(x);
            UiL = U(x - Delta);
            UiR = U(x + Delta);
            rL = nux*exp(- Bx*(UiL - Ui)/2);
            rR = nux*exp(- Bx*(UiR - Ui)/2);
            iX(middleC) = ix;
            iY(middleC) = ix;
            iL(middleC) = -(rL + rR);
            rL = nux*exp(Bx*(UiL - Ui)/2);
            rR = nux*exp(Bx*(UiR - Ui)/2);
            iX(middleC - 1) = ix;
            iY(middleC - 1) = ix - 1;
            iL(middleC - 1) = rL;
            iX(middleC + 1) = ix;
            iY(middleC + 1) = ix + 1;
            iL(middleC + 1) = rR;
        end
        
        LM = sparse(iX,iY,iL);
        [V,e]=eigs(LM,1,0);
        jstatFPSS = e;
        
%         rateAnal = sqrt(2)/pi*abs(p1-1)/(tau*p1)*exp(-3/4*(p1-1)^2/(p2*tau));
%         rateAnalEff = sqrt(2)/pi*abs(p1-1)/(tau*p1)*exp(-3*(p1-1)^2./(2*p2*p1^2*tau)); %% in p.c.'s 2-D simulation, the D_eff is about 2*D/R^2=2*p2
        
        rateN_measured(itau,ia) = - jstatFPSS;
%         rateA(itau,ia) = rateAnal;
%         rateAeff(itau,ia) = rateAnalEff;
        
    end
end
end

