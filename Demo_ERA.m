% function [froErr, mu, sigma] = ErrAnalysis(interp_points, m, ns)
close all; clear all; clc

%% Define the parameters
    format long g
    m = 4; 
    ns = 100;
    interp_points = 2^13;

%% loading the cps values and cps times 
    ch1_cps_values = load('ch01PreCO2_cps_values');
    ch1_cps_times = load('ch01PreCO2_cps_times');
   
    ch2_cps_values = load('ch02PreCO2_cps_values');
    ch2_cps_times = load('ch02PreCO2_cps_times');
    
    ch3_cps_values = load('ch03PreCO2_cps_values');
    ch3_cps_times = load('ch03PreCO2_cps_times');
    
    ch4_cps_values = load('ch04PreCO2_cps_values');
    ch4_cps_times = load('ch04PreCO2_cps_times');

%% Create the UNION of all CPS times
    union_cps_times = vertcat([ch1_cps_times ; ch2_cps_times; ch3_cps_times; ch4_cps_times]);
    union_cps_times = struct2cell(union_cps_times); %change from struct to double
    union_cps_times = cell2mat(union_cps_times);
    union_cps_times = unique(union_cps_times);
    union_cps_times = sort(union_cps_times);

    %% Changing data from struct to double
    ch1_cps_values = struct2cell(ch1_cps_values);
    ch1_cps_values = cell2mat(ch1_cps_values);
    ch1_cps_times = struct2cell(ch1_cps_times);
    ch1_cps_times = cell2mat(ch1_cps_times);

    ch2_cps_values = struct2cell(ch2_cps_values);
    ch2_cps_values = cell2mat(ch2_cps_values);
    ch2_cps_times = struct2cell(ch2_cps_times);
    ch2_cps_times = cell2mat(ch2_cps_times);

    ch3_cps_values = struct2cell(ch3_cps_values);
    ch3_cps_values = cell2mat(ch3_cps_values);
    ch3_cps_times = struct2cell(ch3_cps_times);
    ch3_cps_times = cell2mat(ch3_cps_times);

    ch4_cps_values = struct2cell(ch4_cps_values);
    ch4_cps_values = cell2mat(ch4_cps_values);
    ch4_cps_times = struct2cell(ch4_cps_times);
    ch4_cps_times = cell2mat(ch4_cps_times);

%% Computing the minimum distances between z and t_alphas 
    z = union_cps_times;
    D = distances(z);
    minDistance = min(D(:));

%% Create time vector
    dt = 0.01;
    t_end = 300;
    dt_desired = t_end/(interp_points-1);
    t_desired = 0:dt_desired:t_end;

%% Fill Missing Values
    [F1,TF1] = fillmissing(ch1_cps_values,'linear','SamplePoints',ch1_cps_times);
    [F2,TF2] = fillmissing(ch2_cps_values,'linear','SamplePoints',ch2_cps_times);
    [F3,TF3] = fillmissing(ch3_cps_values,'linear','SamplePoints',ch3_cps_times);
    [F4,TF4] = fillmissing(ch4_cps_values,'linear','SamplePoints',ch4_cps_times);

%% Create the data matrix
% Fill missing data in phases 1 and 6
% Interpolation matrix: Vq = interp1(X,V,Xq) interpolates to find Vq, the values of the underlying function V=F(X) at the query points Xq. 
    fun_interp_r1 = interp1(ch1_cps_times,F1,t_desired, 'linear',267);
    fun_interp_r2 = interp1(ch2_cps_times,F2,t_desired, 'linear',267);
    fun_interp_r3 = interp1(ch3_cps_times,F3,t_desired, 'linear',267);
    fun_interp_r4 = interp1(ch4_cps_times,F4,t_desired, 'linear',267);

    fun_interp =[fun_interp_r1; fun_interp_r2; fun_interp_r3; fun_interp_r4];
    x = fun_interp(:,1:end);

%% Create Hankel (data) matrices 
    H = form_H(x,m,interp_points-1,ns);
    H2 = form_H2(x,m,interp_points-1,ns);

%% SVD and rank -r truncation
    [U,S,V] = svd(H,'econ');
    r = rank(S);

%% Call ERA, ypred
    [Ar,Br,Cr] = run_ERA(H,H2,m,r);
    sys = ss(Ar,Br,Cr,0,dt_desired);
    u = zeros(length(t_desired),m);
    u(1) = 1;
    [ypred,t_desired] = lsim(sys,u,t_desired); 

%% Visualize the true/infe
    figure(1);
    for i = 1:m
        subplot(m,1,i);
        plot(t_desired(2:end),x(i,2:end),'k:.','LineWidth',1.5); hold on;
        plot(t_desired(2:end),ypred(2:end,i),'r--','LineWidth',1.5);hold on;
        ylabel(['CH' num2str(i, '%02d')]);
        sgtitle('Awake rat 16-PreCO2');
    end
    xlabel('Time');
    
%% Frobenius Error
    froErr = norm((x(:,2:end)-ypred(2:end,:)')/100,"fro");

%% Auxiliary Functions
% Creating the Hankel matrix H
    function H = form_H(x,m,n,ns)
        H = zeros(ns*m,n-ns);
        for j = 1:n-ns
            H_new = [];
            for i = j:j+ns-1
                H_new = [H_new; x(:,i)];
            end
            H(:,j) = H_new;
        end
    end

% Creating the shifted Hankel matrix H2
    function H2 = form_H2(x,m,n,ns)
       H2 =  zeros(ns*m,n-ns);
        for j = 2:n-ns+1
           H2_new = [];
           for i = j:j+ns-1
               H2_new = [H2_new; x(:,i)];
           end
           H2(:,j-1) = H2_new;
        end
    end

% ERA algorithm to infer A,B, and C
    function [Ar,Br,Cr] = run_ERA(H,H2,m,r)
        [U,S,V] = svd(H,'econ');
        Sigma = S(1:r,1:r);
        Ur = U(:,1:r);
        Vr = V(:,1:r);
        Ar = Sigma^(-.5)*Ur'*H2*Vr*Sigma^(-.5);
        Br = Sigma^(-.5)*Ur'*H(:,1:m);
        Cr = H(1:m,:)*Vr*Sigma^(-.5);
    end

    % Computing the distance between all pairs
    function D = distances(z)
    %  Input:
    %    real A(N,N), the one-step distance from time(i) to time(j).
    %
    %  Output:
    %    real D(N,N), the shortest total distance from time(i) to time(j).

      n = size (z, 1);
      D = z;

  for k = 1 : n
    for j = 1 : n
      for i = 1 : n
        D(i,j) = min ( D(i,j), D(i,k) + D(k,j) );
      end
    end
  end

  return
end
