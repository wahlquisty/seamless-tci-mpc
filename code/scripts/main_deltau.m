% Demonstrate how punishing delta u in cost function affects simulation to
% prevent ringing in the control signal
% Use MPC for control
% Date: 2024-10-23

clear all; clc; close all

%% Simulate with disturbance without delta u taken into consideration

disturbance = 1; % disturbance affecting DOH
hasdeltau = 0; % without adding cost function term

[t_nodeltau, trueDOH_nodeltau, ymeas_nodeltau, u_nodeltau] = simulate_deltau(disturbance,hasdeltau);

% Plot results
figure(1)
subplot(3,1,1)
plot(t_nodeltau./60,u_nodeltau)
xlabel('Time (min)')
ylabel('Dose (mg/s)')
title('u')

subplot(3,1,2)
plot(t_nodeltau./60,trueDOH_nodeltau)
xlabel('Time (min)')
ylabel('DoH (BIS)')
title('True DoH in patient')
ylim([0,100])

subplot(3,1,3)
plot(t_nodeltau./60,ymeas_nodeltau)
xlabel('Time (min)')
ylabel('Measurement (BIS)')
title('Measurement (BIS)')
ylim([0,100])

%% Simulate with disturbance with delta u taken into consideration

disturbance = 1; % disturbance affecting DOH
hasdeltau = 1; % with adding cost function term

[t_deltau, trueDOH_deltau, ymeas_deltau, u_deltau] = simulate_deltau(disturbance,hasdeltau);

% Plot results
figure(2)
subplot(3,1,1)
plot(t_deltau./60,u_deltau)
xlabel('Time (min)')
ylabel('Dose (mg/s)')
title('u')

subplot(3,1,2)
plot(t_deltau./60,trueDOH_deltau)
xlabel('Time (min)')
ylabel('DoH (BIS)')
title('True DoH in patient')
ylim([0,100])

subplot(3,1,3)
plot(t_deltau./60,ymeas_deltau)
xlabel('Time (min)')
ylabel('Measurement (BIS)')
title('Measurement (BIS)')
ylim([0,100])


%% cut first 8 min for visibility
simulation_settings;

tcut = 13*60/h;
tnew = t_deltau(tcut:end) - t_deltau(tcut); % time zero is at the time of when induction starts
trueDOHnew = trueDOH_deltau(tcut:end);
unew = u_deltau(tcut:end);

trueDOHnew_nodeltau = trueDOH_nodeltau(tcut:end);
unew_nodeltau = u_nodeltau(tcut:end);

%% Save data to file

data = [tnew'./60 unew_nodeltau trueDOHnew_nodeltau unew trueDOHnew];
% dlmwrite('csv/mpc_deltau.csv',data);


