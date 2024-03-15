clear all; close all; clc;
eval( webread('http://b.link/rbeast') )  

%% FIELDS for your RAT
    format long g
    load FM_RAT_PT16
    n = 940245;   % Number of Signals per phase
    m = 4;        % Number of Channels

%% Structure called Anesthetized_rat_AW1 with fields of all phases
    AW16.PreCO2 = FM_RAT_PT16(:,1:n);
    AW16.PostCO2_Phase01 = FM_RAT_PT16(:,n+1:2*n);
    AW16.PostCO2_Phase02 = FM_RAT_PT16(:,2*n+1:3*n);
    AW16.PostCO2_Phase03 = FM_RAT_PT16(:,3*n+1:4*n);
    AW16.PostCO2_Phase04 = FM_RAT_PT16(:,4*n+1:5*n);
    AW16.PostCO2_Phase05 = FM_RAT_PT16(:,5*n+1:6*n);
    AW16.PostCO2_Phase06 = FM_RAT_PT16(:,6*n+1:end);

%% Create the data matrix
    x = AW16.PreCO2(1,:);  % change phase/ channel

%% Create time vector
    N = length(x);
    dt = 1/2999;
    time = dt:dt:N*dt;

%% Beast matrix: % x = fun_interp;
% 1. split my data in chunks of 1 second 
    num_points_per_second = 2999;
    num_of_seconds = floor(N/num_points_per_second);

%  2. initialize the variables that you want to save     
    all_optimal_number_cps = [];
    all_time_cps = [];
    all_values_cps = [];
    all_R2s = [];
        
 % 3. FOR loop for number_of_seconds   
    for i = 1:5*60 %num_of_seconds
       % 3.1 select_time and select_data 
           selected_time = time((i-1)*2999 + 1:(i)*2999);
           selected_data = x((i-1)*2999 + 1:(i)*2999);
           
       % (OPTIONAL) smooth the data
           selected_data = smoothdata(selected_data);
       
       %  3.2 find the optimal number of tcps in my 1 second chunk
        % select the range for the SEARCH
            optimal_number_search = beast(selected_data,'season','none','tcp.minmax',[5,100]);
            optimal_number_cps = optimal_number_search.trend.ncp_median;
        
       %  3.3 run beast this time with the optimal number of cps
            beast_result = beast(selected_data,'season','none','tcp.minmax',[1,optimal_number_cps]);
       
       %  3.4 get the cp information and R2 
            time_cps     = selected_time(sort(beast_result.trend.cp));
            values_cps   = beast_result.trend.Y(sort(beast_result.trend.cp));
            R2_cps        =  beast_result.R2;
       
            left_point    = beast_result.trend.Y(1);
            right_point  = beast_result.trend.Y(end);
       
       % 3.5 augment time_cps and values_cps with endpoints
            time_cps = [selected_time(1),time_cps,selected_time(end)];
            values_cps = [left_point; values_cps; right_point]';
                
       % 3.6 collect my results
             all_optimal_number_cps = [all_optimal_number_cps,optimal_number_cps];
             all_time_cps = [all_time_cps, time_cps];
             all_values_cps = [all_values_cps,values_cps];
             all_R2s = [all_R2s,R2_cps];           
    end
    
%% Plot BEAST vs DATA 
    figure
    scatter(time,x,'ok','filled','MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.05)
    xlim([0, all_time_cps(end)]); ylim([240,290]);
    xlabel('Time (seconds)'); set(gcf,'position',[10 10 1200 400]);
    set(gca,'FontSize',18); hold on;
    title([num2str(sum(all_optimal_number_cps)) ' Change Points via BEAST' ]);
    
    % augment time_cps and values_cps with endpoints
    plot(all_time_cps,all_values_cps,'Color','g','LineWidth',4) 
