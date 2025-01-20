% Full simulation scenario, induction, maintenance + two disturbances (one
% affecting DOH and one affecting measurement)
% Optimize R and diagonal Q to minimize MSE and simulate with the optimized
% parameters.
% We assume that the second disturbance is related to poor SQI (SQI
% = 50). Added noise from Andrzej Pawlowski et al. “MPC for Propofol Anesthesia:
% The Noise Issue"
% Plot result

% Date: 2024-10-23
clear all; clc; close all

%% First, simulate induction and disturbance scenario with non-optimal R and Q

Rmin = 0;
Rmax = 1;
Q = 1*eye(4);

addednoise = 1; % added noise

[t, trueDOH, ymeas, ykalman, sqi, u] = simulate_example(Rmin,Rmax,Q,addednoise);

% Plot simulation results
figure(1)
subplot(3,1,1)
plot(t./60,u)
xlabel('Time (min)')
ylabel('Dose (mg/s)')
title('u')

subplot(3,1,2)
plot(t./60,trueDOH)
xlabel('Time (min)')
ylabel('DoH (BIS)')
title('True DoH in patient')
ylim([0,100])

subplot(3,1,3)
plot(t./60,ymeas)
xlabel('Time (min)')
ylabel('Measurement (BIS)')
title('Measurement (BIS)')
ylim([0,100])

% %% Compute MSE between kalman prediction [BIS] and trueDOH [BIS]
% 
% v = [Rmin;Rmax;diag(Q)];
% 
% computeMSEexample(v) % sum

%% Optimize R and Q for simulation scenario

rng(1) % For reproducibility (other seed than for later simulations)

% Initial guesses from Kalman paper (these are not from the Kalman paper?!)
% Rminguess = 1.3277e+04;
% Rmaxguess = 9.8734e+06;
% Qvec = 1.0e+06*[0.0047,0.1149,5.8082,0.0034]';
Rminguess = 0.771;
Rmaxguess = 1.79;
Qvec = [0.0579, 0.0183, 0.027, 0.000212]';
v0 = [Rminguess;Rmaxguess;Qvec];

options = optimoptions('fmincon','Display','iter','MaxIterations',50);

A = -1*eye(6); % Ax <= b -> x >= 0
b = zeros(6,1);

tic
[optimalvec,fval,exitflag] = fmincon(@computeMSEexample,v0,A,b,[],[],[],[],[],options)
toc

Rminopt = optimalvec(1)
Rmaxopt = optimalvec(2)
Qopt = optimalvec(3:6).*eye(4)

%% Compute MSE for the optimal values obtained

% With noise
% Rminopt =
%    273.6357
% Rmaxopt =
%    19706
% Qopt =
%     0.8543         0       0         0
%          0    0.83         0         0
%          0         0    2.4601*10^4   0
%          0         0         0    0.142*10^4
% fval =
%    10.7689

v = [Rminopt;Rmaxopt;diag(Qopt)];
computeMSEexample(v) % sum mse

%% Simulate with optimal RQ 

[t, trueDOH, ymeas, ykalman, sqi, u] = simulate_example(Rminopt,Rmaxopt,Qopt,addednoise);

% Plot results
figure(1)
subplot(4,1,1)
plot(t./60,u)
xlabel('Time (min)')
ylabel('Dose (mg/s)')
title('u')

subplot(4,1,2)
plot(t./60,trueDOH)
xlabel('Time (min)')
ylabel('DoH (BIS)')
title('True DoH in patient')
ylim([0,100])

subplot(4,1,3)
plot(t./60,ymeas)
xlabel('Time (min)')
ylabel('Measurement (BIS)')
title('Measurement (BIS)')
ylim([0,100])

subplot(4,1,4)
plot(t./60,ykalman)
xlabel('Time (min)')
ylabel('Kalman prediction (BIS)')
title('Kalman prediction (BIS)')
ylim([0,100])

%% cut first five min for filter to converge
simulation_settings_example;

tcut = 5*60/h;
tnew = t(tcut:end) - 2*t(tcut); % time zero is at the time of when induction starts
trueDOHnew = trueDOH(tcut:end);
monitorBISnew = ymeas(tcut:end);
ykalmannew = ykalman(tcut:end);
unew = u(tcut:end);

%% Save data to file

data = [tnew'./60 unew trueDOHnew monitorBISnew ykalmannew];
% dlmwrite('csv/examplesqi.csv',data);

