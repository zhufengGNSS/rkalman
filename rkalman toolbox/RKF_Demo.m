%--------------------------------------------------------------------------
% Name:            RKF_Demo.m
%
% Description:     Application of the robust Kalman filter to the task of 
%                  tracking a sampled Wiener process on noisy measurements.
%
% Author:          Mattia Zorzi
%
% Date:            Agoust 20, 2015
%--------------------------------------------------------------------------

%%  True model
% Process noise variance
Q = 0.8;
% Measurement noise variance
R = 2;

%% Data generation

% Number of iterations  
N = 300;
% True State and measurement
x = zeros(N,1);
y = zeros(N,1);
% True initial state
x(1) = randn;
% First measurement
y(1) = x(1) + sqrt(Q)*randn;
% Update true state and measurements
for i=2:N
    x(i) = x(i-1) + sqrt(Q)*randn;
    y(i) = x(i) + sqrt(R)*randn;
end

%% Nominal noise/process covariance

% Process noise variance
Qi = 2.0;
% Measurement noise variance
Ri = 1.6;

%% Robust Kalman filtering
% Initial apriori state estimate
xp = rkalman(1,[sqrt(Qi) 0],1,[0 sqrt(Ri)],y,10^-1,1);

%% Plot Results
figure
% Plot true state
b = plot(1:N,x(1:N),'b');
hold on
% Plot estimates
r = plot(2:N,xp(1:N-1),'r');
% Plot measurements
g = plot(1:N,y(1:N),'g+');
title(['Robust Kalman Filtering: Tracking a sampled Wiener process']);
legend([b r g],'True Value','RKF Estimate','Measurements');
xlabel('Time')
ylabel('Value')
grid on
