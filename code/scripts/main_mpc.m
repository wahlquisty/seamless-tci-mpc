% Study disturbance rejection by MPC
% Simulation starts in stationarity, at BIS 50
% Two disturbances, one affecting DOH (surgical disturbance, affects the patient state), second one only
% measurement noise (does not affect the patient state)
% Using first patient in data set from Ionescu Clara M. et al. “Robust Predictive Control Strategy
% Applied for Propofol Dosing Using BIS as a Controlled Variable During Anesthesia"

% Date: 2024-10-23
clear all; clc; close all

%% Simulate disturbance - affecting DOH

% Small R and Q similar to having full measurement feedback
R = 1e-6;
Q = 1e-6*eye(4);

disturbance = 1; % 1 means disturbance affecting DOH, 2 means affecting measurement
[t, trueDOH_d1, ymeas_d1, ykalman_d1, u_d1] = simulate_mpc(R,Q,disturbance); % Simulate

% Plot results
figure(1)
subplot(4,1,1)
plot(t./60,u_d1)
xlabel('Time (min)')
ylabel('Dose (mg/s)')
title('u')

subplot(4,1,2)
plot(t./60,trueDOH_d1)
xlabel('Time (min)')
ylabel('DoH (BIS)')
title('True DoH in patient')
ylim([0,100])

subplot(4,1,3)
plot(t./60,ymeas_d1)
xlabel('Time (min)')
ylabel('Measurement (BIS)')
title('Measurement (BIS)')
ylim([0,100])

subplot(4,1,4)
plot(t./60,ykalman_d1)
xlabel('Time (min)')
ylabel('Kalman prediction (BIS)')
title('Kalman prediction (BIS)')
ylim([0,100])

%% Simulate disturbance - affecting measurement

R = 1e-6;
Q = 1e-6*eye(4);

disturbance = 2; % 1 means disturbance affecting DOH, 2 means affecting measurement
[t, trueDOH_d2, ymeas_d2, ykalman_d2, u_d2] = simulate_mpc(R,Q,disturbance);

% Plot results
figure(2)
subplot(4,1,1)
plot(t./60,u_d2)
xlabel('Time (min)')
ylabel('Dose (mg/s)')
title('u')

subplot(4,1,2)
plot(t./60,trueDOH_d2)
xlabel('Time (min)')
ylabel('DoH (BIS)')
title('True DoH in patient')
ylim([0,100])

subplot(4,1,3)
plot(t./60,ymeas_d2)
xlabel('Time (min)')
ylabel('Measurement (BIS)')
title('Measurement (BIS)')
ylim([0,100])

subplot(4,1,4)
plot(t./60,ykalman_d2)
xlabel('Time (min)')
ylabel('Kalman prediction (BIS)')
title('Kalman prediction (BIS)')
ylim([0,100])

%% Compute MSE over DOH and BIS 50

immse(trueDOH_d1,50*ones(length(trueDOH_d1),1))
immse(trueDOH_d2,50*ones(length(trueDOH_d2),1))

%% Save data to file

data = [t'./60 u_d1 trueDOH_d1 ymeas_d1 u_d2 trueDOH_d2 ymeas_d2];
% dlmwrite('csv/mpckalman.csv',data);
