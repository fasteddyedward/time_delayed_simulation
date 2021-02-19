%% Calculating theta_+ and k_trans with hist(theta)
function [V_parameters,theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0)
% theta_stable is theta_+-
theta_stable(1:length(theta))=0;

    theta_stable(1:length(theta))=0;
    k_trans_theta=[];
    theta_plus=0;
    theta_minus=0;
    num_transitions_theta=0;

        %% 2021.1.21 Modification
%     h=histogram(theta);
% %     Values=h.Values;
%     Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
%     Values_reverted=flip(h.Values);
%     Values=(h.Values+Values_reverted)/2;
% %     plot(Bins,Values)
% %     1
        %% 2021.1.21 Mod 2nd (making the histogram symmetric
        theta_symmetric=[theta -theta];
        h=histogram(theta_symmetric);
        Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
        Values=h.Values;
    %% This will not work with 2021.1.21 Mod 2nd because the index for ind(Values==max(Values(Bins>0))) will be wrong for Bins>0
%     if sum(Bins>0)==0
%         theta_plus=0;
%     else
%         %         theta_plus_index=find(Values==max(Values(Bins>0))); % end/2+2 instead of end/2 to avoid theta=0;
%         theta_plus_index=find(Values(Bins>0)==max(Values(Bins>0))); % end/2+2 instead of end/2 to avoid theta=0;
%         if length(theta_plus_index)==1
%             theta_plus=Bins(theta_plus_index);
%         elseif isempty(theta_plus_index)
%             theta_plus=0;
%         else
%             theta_plus=mean(Bins(theta_plus_index));
%             warning('There is more than one maximum for theta_plus')
%         end
%     end
%     
%     if sum(Bins<0)==0
%         theta_minus=0;
%     else
%         %         theta_minus_index=find(Values==max(Values(Bins<0)));  % end/2-2 instead of end/2 to avoid theta=0;
%         theta_minus_index=find(Values(Bins<0)==max(Values(Bins<0)));  % end/2-2 instead of end/2 to avoid theta=0;
%         if  length(theta_minus_index)==1
%             theta_minus=Bins(theta_minus_index);
%         elseif isempty(theta_minus_index)
%             theta_minus=0;
%         else
%             theta_minus=mean(Bins(theta_minus_index));
%             warning('There is more than one maximum for theta_minus')
%         end
%     end
%% 2021.1.21 Modified

    Bins_plus=Bins(length(Bins)/2:end); % The positive half of Bins
    Values_plus=Values(length(Values)/2:end); 
    if sum(Bins_plus>0)==0
        theta_plus=0;
    else
        %         theta_plus_index=find(Values==max(Values(Bins>0))); % end/2+2 instead of end/2 to avoid theta=0;
        theta_plus_index=find(Values_plus==max(Values_plus)); % end/2+2 instead of end/2 to avoid theta=0;
        if length(theta_plus_index)==1
            theta_plus=Bins_plus(theta_plus_index);
        elseif isempty(theta_plus_index)
            theta_plus=0;
        else
            theta_plus=mean(Bins_plus(theta_plus_index));
            theta_plus_index=round(mean(theta_plus_index)); % This theta_plus_index is wrong, but I use it just for convenience;
            warning('There is more than one maximum for theta_plus')
        end
    end
    
%     theta_plus
    theta_minus=-theta_plus;

%     1
    %% Start Calculating Transition Rates
    sign_old=0; % Initial 'order parameter' for the orbit. +1 for stable orbit with theta_plus, -1 for stable orbit with theta_miunus
    sign_current=0;
    num_transitions_theta=0;
    k_trans_theta=[];
    %% 2020.2.8 Update: Calculating transitions only after bifurcation; threshold still needs to be specified
    if theta_plus>0.3
        for k=1:length(theta)
            if theta(k)>theta_plus
                sign_current=1;
                %             theta_stable(k)=theta_plus;
            elseif theta(k)<theta_minus
                sign_current=-1;
                %             theta_stable(k)=theta_minus;
            end
            % Updating sign: sign_old -> sign_current
            if sign_current*sign_old==-1
                num_transitions_theta=num_transitions_theta+1;
                k_trans_theta=[k_trans_theta k];
            end
            
            if sign_current==1
                theta_stable(k)=theta_plus;
            elseif sign_current==-1
                theta_stable(k)=theta_minus;
            end
            
            sign_old=sign_current;
        end
    end
    
    
    %% Start fitting the real effective potential
%     figure(98)
%     histogram(theta)
%     
%     figure(99)
%     plot(Bins,Values/sum(Values))
%    
%     
%     ln_p=log(Values/sum(Values));
%     figure(100);clf;
%     plot(Bins,-ln_p)
    
    %% Cutoff values: Since entries that are too far away from the bell curve have high stat error,
    %% We only allow entries > 0.05*(entries at theta_plus) be included
    err=5*10^-2; 
    Bins_cutoff=Bins(Values>err*Values_plus(theta_plus_index));
    Values_cutoff=Values(Values>err*Values_plus(theta_plus_index));
    
%     figure(101);clf
%     plot(Bins_cutoff,Values_cutoff/sum(Values_cutoff))
    
    ln_p_cutoff=log(Values_cutoff/sum(Values_cutoff));
%     figure(102);clf;
%     plot(Bins_cutoff,-ln_p_cutoff)

    %% Fitting Potential 
    
    [V_parameters, ~] = potential_fit(Bins_cutoff, -ln_p_cutoff);

    
    
%     V=@(x) fitresult.p0 + fitresult.p2.*x.^2 + fitresult.p4.*x.^4 + fitresult.p6.*x.^6 + fitresult.p8.*x.^8; % V(theta) = U(theta)/(k_B T)
    
%     V_prime_prime = @(x) 2*fitresult.p2 + 12*fitresult.p4.*x.^2 + 30*fitresult.p6.*x.^4 + 56*fitresult.p8.*x.^6; % V''(theta)=U''(theta)/(k_B T)
    





function [fitresult, gof] = potential_fit(Bins, ln_p)

%% Fit: 'Potential'.
[xData, yData] = prepareCurveData( Bins, ln_p );

% Set up fittype and options.
ft = fittype( 'p0+p2*x^2+p4*x^4+p6*x^6+p8*x^8', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.321136379000634 0.754272380536349 0.815593506580009 0.0863539417530238 0.030584328217947];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% 
% legend( h, 'ln_p vs. Bins', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'Bins', 'Interpreter', 'none' );
% ylabel( 'ln_p', 'Interpreter', 'none' );
% grid on


end

end