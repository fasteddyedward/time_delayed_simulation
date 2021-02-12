%% Calculating theta_+ and k_trans with hist(theta)
function [theta_stable,k_trans_theta,theta_plus,theta_minus,num_transitions_theta]=find_theta_plus(theta,theta_0)
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
        %% 2021.1.21 Mod 2nd
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

    Bins_plus=Bins(length(Bins)/2:end);
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
    

end