nth_take=7
% delta_t_matrix=[0.01 0.1 1 5 10]
% T_matrix=[0.01 0.1 1 10]
% v_0_matrix=[0.01 0.1 1 10]
delta_t_matrix=[2]
T_matrix=[1]
v_0_matrix=[7:0.1:9]
intrinsic_delay=0 % Intrinsic delay
num_transitions_matrix=[];
for delta_t_index=1:length(delta_t_matrix)
    for T_index=1:length(T_matrix)
        for v_0_index=1:length(v_0_matrix)
            if 1
                %             if ismember(nth_take,nth_interest)
close all
%% Output File Name
% movie_name=['2020.11.19,dt=10e-3 take ',num2str(nth_take)];
movie_name=['2020.11.20,dt=10e-3 take ',num2str(nth_take),', T=',num2str(T_matrix(T_index)),', v_0=',num2str(v_0_matrix(v_0_index)),', delta_t=',num2str(delta_t_matrix(delta_t_index))];
% movie_name=['test3']
[movie_name,'.mat'];
load([movie_name,'.mat'],'num_transitions')
num_transitions_matrix=[num_transitions_matrix num_transitions];

            end
            nth_take=nth_take+1
        end
        nth_take=nth_take+1
    end
    nth_take=nth_take+1
end

%%

plot(v_0_matrix,num_transitions_matrix)