%% Calculating theta_+ and k_trans with hist(theta)
function [theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0)
% theta_stable is theta_+-
theta_stable(1:length(theta))=0;
if theta_0<0.9
    theta_stable(1:length(theta))=0;
    k_trans_theta=[];
    theta_plus=0;
    theta_minus=0;
    num_transitions_theta=0;
else
    h=histogram(theta);
    Values=h.Values;
    Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
%     theta_plus_index=find(Values==max(Values(end/2:end))); % end/2+2 instead of end/2 to avoid theta=0;
%     theta_minus_index=find(Values==max(Values(1:end/2)));  % end/2-2 instead of end/2 to avoid theta=0;
    if sum(Bins>0)==0
        theta_plus=0;
    else
        theta_plus_index=find(Values==max(Values(Bins>0))); % end/2+2 instead of end/2 to avoid theta=0;
        if length(theta_plus_index)==1
            theta_plus=Bins(theta_plus_index);
        elseif isempty(theta_plus_index)
            theta_plus=0;
        else
            theta_plus=mean(Bins(theta_plus_index));
            warning('There is more than one maximum for theta_plus')
        end
    end
    
    if sum(Bins<0)==0
        theta_minus=0;
    else
        theta_minus_index=find(Values==max(Values(Bins<0)));  % end/2-2 instead of end/2 to avoid theta=0;
        if  length(theta_minus_index)==1
            theta_minus=Bins(theta_minus_index);
        elseif isempty(theta_minus_index)
            theta_minus=0;
        else
            theta_minus=mean(Bins(theta_minus_index));
            warning('There is more than one maximum for theta_minus')
        end
    end
    %% Start Calculating Transition Rates
    sign_old=0; % Initial 'order parameter' for the orbit. +1 for stable orbit with theta_plus, -1 for stable orbit with theta_miunus
    sign_current=0;
    num_transitions_theta=0;
    k_trans_theta=[];
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
end