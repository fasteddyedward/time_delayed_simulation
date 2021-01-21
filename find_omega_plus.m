function [omega_stable,k_trans_omega,omega_plus,omega_minus,num_transitions]=find_omega_plus(omega,theta_0)
% omega_stable is omega_+-
omega_stable(1:length(omega))=0;
if theta_0<0.9
    omega_stable(1:length(omega))=0;
    k_trans_omega=[];
    omega_plus=0;
    omega_minus=0;
    num_transitions=0;
else
    %% 2021.1.21 Modified
    
    %     h=histogram(omega);
    %     Values=h.Values;
    %     Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    
    omega_symmetric=[omega -omega];
    h=histogram(omega_symmetric);
    Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    Values=h.Values;


        
%         %%
%     if sum(Bins>0)==0
%         omega_plus=0;
%     else
%         omega_plus_index=find(Values==max(Values(Bins>0))); % end/2+2 instead of end/2 to avoid theta=0;
%         if length(omega_plus_index)==1
%             omega_plus=Bins(omega_plus_index);
%         elseif isempty(omega_plus_index)
%             omega_plus=0;
%         else
%             omega_plus=mean(Bins(omega_plus_index));
%             warning('There is more than one maximum for theta_plus')
%         end
%     end
%     
%     if sum(Bins<0)==0
%         omega_minus=0;
%     else
%         omega_minus_index=find(Values==max(Values(Bins<0)));  % end/2-2 instead of end/2 to avoid theta=0;
%         if  length(omega_minus_index)==1
%             omega_minus=Bins(omega_minus_index);
%         elseif isempty(omega_minus_index)
%             omega_minus=0;
%         else
%             omega_minus=mean(Bins(omega_minus_index));
%             warning('There is more than one maximum for theta_minus')
%         end
%     end
%% 2021.1.21 Modified

    Bins_plus=Bins(length(Bins)/2:end);
    Values_plus=Values(length(Values)/2:end);
    if sum(Bins_plus>0)==0
        omega_plus=0;
    else
        %         theta_plus_index=find(Values==max(Values(Bins>0))); % end/2+2 instead of end/2 to avoid theta=0;
        theta_plus_index=find(Values_plus==max(Values_plus)); % end/2+2 instead of end/2 to avoid theta=0;
        if length(theta_plus_index)==1
            omega_plus=Bins_plus(theta_plus_index);
        elseif isempty(theta_plus_index)
            omega_plus=0;
        else
            omega_plus=mean(Bins_plus(theta_plus_index));
            warning('There is more than one maximum for theta_plus')
        end
    end
    
%     theta_plus
    omega_minus=-omega_plus;

%     1
    %% Start Calculating Transition Rates
    sign_old=0; % Initial 'order parameter' for the orbit. +1 for stable orbit with theta_plus, -1 for stable orbit with theta_miunus
    sign_current=0;
    num_transitions=0;
    k_trans_omega=[];
    for k=1:length(omega)
        if omega(k)>omega_plus
            sign_current=1;
            %             theta_stable(k)=theta_plus;
        elseif omega(k)<omega_minus
            sign_current=-1;
            %             theta_stable(k)=theta_minus;
        end
        % Updating sign: sign_old -> sign_current
        if sign_current*sign_old==-1
            num_transitions=num_transitions+1;
            k_trans_omega=[k_trans_omega k];
        end
        
        if sign_current==1
            omega_stable(k)=omega_plus;
        elseif sign_current==-1
            omega_stable(k)=omega_minus;
        end
        
        sign_old=sign_current;
    end
    
end
end
