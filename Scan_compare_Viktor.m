%% Same as modified_average_rate_and_decay_time_delay_angle.m. Functionalized to compare the transition rate with Pin_Chuan's in 
%% Compare_1D_2D_pc_Vik.m

function Scan_compare_Viktor(D_0,dt,tV_max,Av,tau_matrix,runs,file_name_vik)
% close all
rng('shuffle')

% dt = 0.01;
% dt=0.1
tV=0:dt:tV_max;
% tV = 0:dt:50000;
% tV = 0:dt:500000;
% tV = 0:dt:1000;
% runs = 1;

dDtT = 0.1;
transitionDTimes = 0:dDtT:1000;
dtT = 1;
transitionTimes = 0:dtT:length(tV);

%Av is theta_0
% Av = linspace(1.001,1.2,15);
% Av = linspace(1.1,1.6,15);
AratesV = zeros(length(tau_matrix),length(Av));
ArelaxV = zeros(length(tau_matrix),length(Av));
sigma2V = zeros(length(tau_matrix),length(Av));
TeffV = zeros(length(tau_matrix),length(Av));
DeffV = zeros(length(tau_matrix),length(Av));
num_transitions_matrix=zeros(length(tau_matrix),length(Av));
for itau=1:length(tau_matrix)
    for iAv = 1:1:length(Av)

        fprintf('Omega tau %i/%i\n\n',iAv,length(Av));

        p1 = Av(iAv); %V0*tau/2R % This is just theta_0
        %     p2 = 0.01; %D*tau^2/4R^2 % But Viktor says it's actually D_0/(r^2), I think r is 2R
        p2=D_0/10^2;
        %     tau = 1;
        if length(tau_matrix)>1
            tau=tau_matrix(itau);
        else
            tau=tau_matrix;
        end
        for iruns = 1:runs

            PtDT = zeros(1,length(transitionDTimes));
            PtT = zeros(1,length(transitionTimes));

            if mod(iruns,10) == 0
                fprintf('done %i/%i\n',iruns,runs);
            end

            ntau = round(tau/dt);
            xV = zeros(1,length(tV));

            xini = p1;

            for it = 1:length(tV)-1

                if it - ntau < 1
                    xd = xini; %xV(1);
                else
                    xd = xV(it - ntau);
                end

                xV(it+1) = xV(it) + p1/tau*sin(xV(it) - xd)*dt + sqrt(2*p2*dt)*randn;

            end

            Vim = (xV(2:end)-xV(1:end-1))/dt; % this looks like omega=diff(phi)/dt

            Vavr = xV(ntau + 1:end) - xV(1:end-ntau); % this looks like phi(t+delta_t)-phi(t)=theta(t)

            figure(1);
            subplot(1,2,1)
            VH = histogram(Vim);
            dOmegaH = (VH.BinEdges(2) - VH.BinEdges(1))/2;
            OmegaH = VH.BinEdges(1:end-1) + dOmegaH;
            valH = VH.Values;        
            xlabel('$\omega$','Interpreter','latex');
            ylabel('$P(\omega)$','Interpreter','latex');
            subplot(1,2,2)
            VavrH = histogram(Vavr);


            theta_hist=histogram(diff(xV));
            dBins = (theta_hist.BinEdges(2) - theta_hist.BinEdges(1))/2;
            Bins = theta_hist.BinEdges(1:end-1) + dBins;
            Values=theta_hist.Values;
    %         xlabel('$\omega$','Interpreter','latex'); % Should be theta instead
    %         ylabel('$P(\omega)$','Interpreter','latex');
            xlabel('$\theta$','Interpreter','latex'); 
            ylabel('$P(\theta)$','Interpreter','latex');

            %% get effective temperature
    %         [fitresult,gof]=fit_Gaussian(OmegaH,valH,'no');
    %         mu_pc=fitresult.b1;
    %         sigma_pc=fitresult.c1/sqrt(2);
    %         D_eff=sigma_pc^2*dt/(2)
            %% Calculate for theta
    %         Bins=Bins/trapz(Bins,Values) 
    %         [fitresult,gof]=fit_Gaussian(Bins,Values,'no');
    %         mu_pc=fitresult.b1;
    %         sigma_pc=fitresult.c1/sqrt(2);
    %         D_eff=sigma_pc^2/(2*dt)


            valH = valH/trapz(OmegaH,valH);
            sigma2 = trapz(OmegaH,OmegaH.^2.*valH) - (trapz(OmegaH,OmegaH.*valH))^2; % same as sigma_pc^2 if we instead set Dcalc=sigma^2/(2*dt)
            Dcalc = sigma2*dt/2;
            Teff = Dcalc*p1/2;


            %% get stacionary values of omega
            VavrH = histogram(Vavr);
            [pmax,valmax] = max(VavrH.Values);
            OmegaStc = abs((VavrH.BinEdges(valmax+1) + VavrH.BinEdges(valmax+1))/2);
            Vdichotomic = zeros(1,length(Vavr));
            numberOfTransitions = 0;
            for it = 1:length(Vavr)
                if Vavr(it) <  0
                    if it > 1 % because we don't know sing at it==1
                        if sing == 1
                            numberOfTransitions = numberOfTransitions + 1;
                            indicesOfTransitions(numberOfTransitions) = it; % something like k_trans
                        end
                    end
                    sing =  - 1;
                else
                    if it > 1
                        if sing == -1
                            numberOfTransitions = numberOfTransitions + 1;
                            indicesOfTransitions(numberOfTransitions) = it;
                        end
                    end
                    sing = 1;
                end
                Vdichotomic(it) = OmegaStc*sing;
            end

            if numberOfTransitions > 0

                for iNt = 1:numberOfTransitions-1

                    it0 = indicesOfTransitions(iNt);
                    itMax = it0;
                    itMin = it0;
                    stop = 1;
                    while abs(Vavr(itMax)) < OmegaStc && stop
                        if itMax == length(Vavr)
                            stop = 0;
                        else
                            itMax = itMax + 1;
                        end
                    end
                    stop = 1;
                    while abs(Vavr(itMin)) < OmegaStc && stop
                        if itMin == 1
                            stop = 0;
                        else
                            itMin = itMin - 1;
                        end
                    end

                    if itMax >= indicesOfTransitions(iNt + 1)
                        numberOfTransitions = numberOfTransitions - 1;
                        for itk = it0:indicesOfTransitions(iNt + 1) - 1
                            Vdichotomic(itk) = Vdichotomic(itk - 1);
                        end
                    end

                end

                indicesOfTransitions = zeros(1, numberOfTransitions);
                numberOfTransitions = 0;
                for it = 1:length(Vdichotomic)
                    if Vdichotomic(it) <  0
                        if it > 1
                            if sing == 1
                                numberOfTransitions = numberOfTransitions + 1;
                                indicesOfTransitions(numberOfTransitions) = it;
                            end
                        end
                        sing =  - 1;
                    else
                        if it > 1
                            if sing == -1
                                numberOfTransitions = numberOfTransitions + 1;
                                indicesOfTransitions(numberOfTransitions) = it;
                            end
                        end
                        sing = 1;
                    end
                end

                durationsOfTransitions = zeros(1,numberOfTransitions);
                rates = zeros(1,numberOfTransitions);

                for iNt = 1:numberOfTransitions

                    it0 = indicesOfTransitions(iNt);
                    if iNt > 1
                        rates(iNt) = tV(it0) - tV(indicesOfTransitions(iNt - 1));
                    else
                        rates(iNt) = tV(it0) - tV(1);
                    end

                    itMax = it0;
                    itMin = it0;
                    stop = 1;
                    while abs(Vavr(itMax)) < OmegaStc && stop
                        if itMax == length(Vavr)
                            stop = 0;
                        else
                            itMax = itMax + 1;
                        end
                    end
                    stop = 1;
                    while abs(Vavr(itMin)) < OmegaStc && stop
                        if itMin == 1
                            stop = 0;
                        else
                            itMin = itMin - 1;
                        end
                    end

                    durationsOfTransitions(iNt) = tV(itMax) - tV(itMin);
                end

            end

            %% collect data
            for iNt = 1:numberOfTransitions
                DtT = durationsOfTransitions(iNt);
                tT = rates(iNt);

                iDtT = floor(DtT/dDtT) + 1;
                itT = floor(tT/dtT) + 1;

                PtDT(iDtT) = PtDT(iDtT) + 1;
                PtT(itT) = PtT(itT) + 1;
            end

        end

        PtDT = PtDT/sum(PtDT)/dDtT;
        PtT = PtT/sum(PtT)/dtT;

        ArelaxV(itau,iAv) = trapz(transitionDTimes, transitionDTimes.*PtDT);
        AratesV(itau,iAv) = trapz(transitionTimes, transitionTimes.*PtT);

        TeffV(itau,iAv) = Teff;
        DeffV(itau,iAv) = Dcalc;
    %     D_eff_theta(iAv)=D_eff;
        sigma2V(itau,iAv) = sigma2;
        num_transitions_matrix(itau,iAv)=numberOfTransitions;
        fprintf('\n');
    end
end

clear Vim Vavr Vdichotomic tV transitionTimes PtT

% nameSave = ['rate_p2_',num2str(p2),'_tau_',num2str(tau),'.mat'];

% save(nameSave)

%% process data
iF = @(x) (sign(x)+1)./(sign(x)+1);

gammaEff = TeffV./DeffV;

rateA = zeros(length(tau_matrix), length(Av));
rateAeff = zeros(length(tau_matrix), length(Av));
rateN = zeros(length(tau_matrix), length(Av));


%% I guess this is for solving the formula analytically
%% The tau is not changed for solving  the FPE!
for itau=1:length(tau_matrix)
    for ia = 1:1:length(Av)

        % parameters
        % tau = 1;
        tau=tau_matrix(itau);
        p1 = Av(ia);
        %p2 = 0.005;

        % potential
        thp2 = 6*(p1-1)/p1;
        Dx = 4*p2/p1^2;

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

%% This part I don't know what they mean, I guess just getting rid of NaN values.
%% Maybe they don't really matter for plotting.


    AvN_plot = Av;
    AratesV_plot = AratesV;
    rateN_plot = rateN;
    rateA_plot = rateA;
    rateAeff_plot = rateAeff;


% if ~isempty(find(isnan(AratesV), 1))
%     iAmax = find(isnan(AratesV), 1)-1;
%     AvN_plot = Av(1:1:iAmax);
%     AratesV_plot = AratesV(1:1:iAmax);
%     rateN_plot = rateN(1:1:iAmax);
%     rateA_plot = rateA(1:1:iAmax);
%     rateAeff_plot = rateAeff(1:1:iAmax);
% else
%     AvN_plot = Av;
%     AratesV_plot = AratesV;
%     rateN_plot = rateN;
%     rateA_plot = rateA;
%     rateAeff_plot = rateAeff;
% end


%% 
% clear iL iX iY LM OmegaH PtDT transitionDTimes V xV
% close all
% save([file_name_vik,'.mat'])
save([file_name_vik,'.mat'],'AvN_plot','AratesV_plot','rateA_plot','rateAeff_plot','rateN_plot','iF')
end