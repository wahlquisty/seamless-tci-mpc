% Simulate induction with and without integral action to demonstrate the
% need for using integral action

% Date: 2024-10-23

clear all; clc; close all

%% Simulate induction and disturbance scenario - No integral action

Rmin = 0;
Rmax = 1;
Q = 1*eye(4);

hasintegralaction = 0; % no integral action

[t, trueDOH_noint, ymeas_noint, ykalman, u_noint] = simulate_integral(Rmin,Rmax,Q,hasintegralaction);

% Plot results
figure(1)
subplot(3,1,1)
plot(t./60,u_noint)
xlabel('Time (min)')
ylabel('Dose (mg/s)')
title('u')

subplot(3,1,2)
plot(t./60,trueDOH_noint)
xlabel('Time (min)')
ylabel('DoH (BIS)')
title('True DoH in patient')
ylim([0,100])

subplot(3,1,3)
plot(t./60,ymeas_noint)
xlabel('Time (min)')
ylabel('Measurement (BIS)')
title('Measurement (BIS)')
ylim([0,100])


%% With integral action

hasintegralaction = 1;

[t, trueDOH_int, ymeas_int, ykalman, u_int] = simulate_integral(Rmin,Rmax,Q,hasintegralaction);

% Plot results
figure(2)
subplot(3,1,1)
plot(t./60,u_int)
xlabel('Time (min)')
ylabel('Dose (mg/s)')
title('u')

subplot(3,1,2)
plot(t./60,trueDOH_int)
xlabel('Time (min)')
ylabel('DoH (BIS)')
title('True DoH in patient')
ylim([0,100])

subplot(3,1,3)
plot(t./60,ymeas_int)
xlabel('Time (min)')
ylabel('Measurement (BIS)')
title('Measurement (BIS)')
ylim([0,100])


%% cut first five min for filter to converge
simulation_settings_example;

tcut = 5*60/h;
tnew = t(tcut:end) - 2*t(tcut); % time zero is at the time of when induction starts
trueDOHnew_noint = trueDOH_noint(tcut:end);
unew_noint = u_noint(tcut:end);

trueDOHnew_int = trueDOH_int(tcut:end);
unew_int = u_int(tcut:end);

%% Save data to file

data = [tnew'./60 unew_noint trueDOHnew_noint unew_int trueDOHnew_int];
% dlmwrite('csv/example_integralaction.csv',data);

