% Study disturbance rejection by TCI
% Simulation starts in stationarity, at BIS 50
% Two disturbances, one affecting DOH (surgical disturbance, affects the patient state), second one only
% measurement noise (does not affect the patient state)
% Using first patient in data set from Ionescu Clara M. et al. “Robust Predictive Control Strategy
% Applied for Propofol Dosing Using BIS as a Controlled Variable During Anesthesia"

% Date: 2024-10-23
clear all; clc; close all

%% Simulate disturbance - affecting DOH

disturbance = 1; % Disturbance affecting DOH
[t, trueDOH_d1, ymeas_d1, u_d1] = simulate_tci(disturbance); % Simulate

% Plot results
figure(1)
subplot(3,1,1)
plot(t./60,u_d1)
xlabel('Time (min)')
ylabel('Dose (mg/s)')
title('u')

subplot(3,1,2)
plot(t./60,trueDOH_d1)
xlabel('Time (min)')
ylabel('DoH (BIS)')
title('True DoH in patient')
ylim([0,100])

subplot(3,1,3)
plot(t./60,ymeas_d1)
xlabel('Time (min)')
ylabel('Measurement (BIS)')
title('Measurement (BIS)')
ylim([0,100])

% Save data to file
data = [t'./60 u_d1 trueDOH_d1 ymeas_d1];
% dlmwrite('csv/tcidist1.csv',data);

%% Simulate disturbance - affecting measurement

disturbance = 2; % Disturbance affecting measurement
[t, trueDOH_d2, ymeas_d2, u_d2] = simulate_tci(disturbance); % Simulate

% Plot results
figure(2)
subplot(3,1,1)
plot(t./60,u_d2)
xlabel('Time (min)')
ylabel('Dose (mg/s)')
title('u')

subplot(3,1,2)
plot(t./60,trueDOH_d2)
xlabel('Time (min)')
ylabel('DoH (BIS)')
title('True DoH in patient')
ylim([0,100])

subplot(3,1,3)
plot(t./60,ymeas_d2)
xlabel('Time (min)')
ylabel('Measurement (BIS)')
title('Measurement (BIS)')
ylim([0,100])

% Save data to file
data = [t'./60 u_d2 trueDOH_d2 ymeas_d2];
% dlmwrite('csv/tcidist2.csv',data);

%% Compute MSE over DOH and BIS 50

immse(trueDOH_d1,50*ones(length(trueDOH_d1),1))
immse(trueDOH_d2,50*ones(length(trueDOH_d2),1))