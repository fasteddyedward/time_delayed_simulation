clear
close all
%% Read diffusion.tdms if haven't done so
    % [ConvertedData,ConvertVer,ChanNames]=convertTDMS(true,'diffusion.tdms');


    A=ConvertedData.Data.MeasuredData;
    Data(1:48000,1:20)=0;
    Name(1:48000,1:20)=0;
    % B(1:48000,9)=A(9).Data
    for k=16:20
        Data(1:48000,k)=A(k).Data;
    %     Name(1:48000,k)=A(k).Name;
    end
    Data(1:48000,16:20);
    % x y in pixels: 0.06284 um/pixes
    pos_x=Data(1:48000,16)*0.06284; % mu m;
    pos_y=Data(1:48000,17)*0.06284; % mu m;
    Position_out_Nr=Data(1:48000,18);
    Position_out_Frame_nr=Data(1:48000,19);
    Position_out_function=Data(1:48000,20);

    x_matrix=[];
    y_matrix=[];
    for k=1:8
        x_matrix(k,:)=pos_x(k:8:48000);
        y_matrix(k,:)=pos_y(k:8:48000);
    end

    save('diffusion.mat')
%%
load('diffusion.mat')
%%
% What are the length scale of x? mu m or mm?
%% Analyzing data
D_matrix=[]
for k=1:8
close all

x=x_matrix(k,:);
y=y_matrix(k,:);

plot_hist_fit_T_eff='yes'    
dt=0.05 %ms % I see that in the tdms file: 0.05 ms for each frame, but then the each frame has 6 sets of position data
    %% Calculating sigma_x

%     figure(1)
    h=histogram(diff(x));
    Values=h.Values;
    Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(Bins,Values,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_x=fitresult.b1;
    sigma_x=fitresult.c1/sqrt(2);
    D_x=sigma_x^2/(2*dt);
    %% Calculating sigma_y
%     figure(2)
    h=histogram(diff(y));
    Values=h.Values;
    Bins=h.BinEdges(1:end-1)+0.5*h.BinWidth;
    % plot(x,y);
    [fitresult,gof]=fit_Gaussian(Bins,Values,plot_hist_fit_T_eff);
    % a1*exp(-((x-b1)/c1)^2)
    mu_y=fitresult.b1;
    sigma_y=fitresult.c1/sqrt(2);
    D_y=sigma_y^2/(2*dt);
    %%    
    D=(D_x+D_y)/2;
    D_matrix=[D_matrix D];
%     k
%     pause
    
end
%%
% D_matrix(4)=[]
hist(D_matrix)
% plot(D_matrix)

    
    
    