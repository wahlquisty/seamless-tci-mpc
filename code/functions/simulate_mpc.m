function [t, trueDOH, ymeas, ykalman, u] = simulate_mpc(R,Q,dist)
%SIMULATE_MPC Simulate MPC under stationary conditions.
%
%   [t, trueDOH, ymeas, u] = SIMULATE_MPC(dist) in stationarity (BIS = 50)
%   Simulates two types of disturbances
%   - If dist = 1, the disturbance affects the Degree of Hypnosis (DOH),
%     which affects to the fourth state variable.
%   - If dist = 2, the disturbance affects the measurement, introducing 
%     noise into the observations.
%   - R and Q are covariances of the Kalman filter (R: scalar; Q: quadratic matrix)
%
%   Outputs:
%   - t: Time vector representing the simulation duration.
%   - trueDOH: True DoH values after simulation.
%   - ymeas: Measurement [BIS]
%   - ykalman: Kalman filter predictions [BIS]
%   - u: Control input values applied during the simulation, computed for
%   the MPC

%% Setup and initialize states, u and y

simulation_settings; % get simulation settings

r = get_refconc(N_horizon,E0,Ce50,gammahill); % get reference concentration for BIS 50

x = x0; % patient state
x_kalman = x; % State vector Kalman filter
x_kalman_m1 = x_kalman;
P0 = eye(n); % initial covariance matrix
P_m1 = P0; % P_m1 = Pn,n-1. Covariance matrix
I4 = eye(n); % identity of size 4

% Save simulation data in vectors
trueCe = zeros(N,1); % True Ce of patient [mg/L]
trueDOH = zeros(N,1); % True DOH of patient [BIS]
ymeas = zeros(N,1); % measurement [BIS]
ykalman = zeros(N,1); % predicted BIS by Kalman filter
u_all = zeros(N,1); % Dosing curves [mg/s]

%% Simulate
for i = 1:N

    % MPC - find control signal. Set up and solve the QP
    [Hqp, fT, A, b] = setupQP(F, G, alpha_u, umax, r, N_horizon, x_kalman_m1);

    [u_i,Jval] = quadprog(Hqp,fT,A,b,[],[],[],[],[],options);
    u = u_i(1); % Save first dose sample obtained from the optimization
    u_all(i) = u;

    % State update patient PKPD model
    x_ = F*x + G*u;
    Ce = H*x;

    trueCe(i) = Ce; % Save true Ce

    % Compute BIS
    doh_beforedist = computeBIS(Ce, E0, Ce50, gammahill); % Convert to BIS

    % Add disturbances
    if dist == 1 % disturbance affecting DOH
        if t(i) == td1
            doh = doh_beforedist + distsize; % add first disturbance to output from PD model
            Ce_dist = computeeffectconc(doh,Ce50,E0,gammahill);
            x_(4) = Ce_dist;  % patient new affected fourth state
        elseif t(i) == td2
            doh = doh_beforedist - distsize; % add first disturbance to output from PD model
            Ce_dist = computeeffectconc(doh,Ce50,E0,gammahill);
            x_(4) = Ce_dist;   % patient new affected fourth state
        else
            doh = doh_beforedist;
        end
        bis_meas = doh;
    else % disturbance affecting measurement
        doh = doh_beforedist;
        if t(i) >= td1 && t(i) <=td2
            bis_meas = doh + distsize; % add first disturbance to output from PD model
            % bis_meas = doh - dist; % add first disturbance to output from PD model
        else
            bis_meas = doh;
        end
    end

    % Save DoH and measurement
    trueDOH(i) = doh;
    ymeas(i) = bis_meas;

    % Kalman filter
    Ce_meas = computeeffectconc(bis_meas,Ce50,E0,gammahill); % feedback concentration inverted through perfect hill function
    % Ce_meas = Ce;

    % R = 0; % Assume no measurement uncertainty. (Rmax = 0.249799032336015)
    L = (P_m1*Ht)/(H*P_m1*Ht + R); % Kalman gain (replaced inv with division since 1dim)
    % Update current uncertainty:
    P = (I4 - L*H)*P_m1*(I4 - L*H)' + L*R*L';
    % Estimate the current state:
    x_kalman = x_kalman_m1 + L*(Ce_meas - H*x_kalman_m1);
    % Predict the next state:
    x_kalman_p1 = F*x_kalman + G*u;
    P_p1 = F*P*Ft + Q;
    % Memory update:
    x_kalman_m1 = x_kalman_p1;
    P_m1 = P_p1;
    % Predict output:
    Ce_kalman = H*x_kalman;
    ykalman(i) = computeBIS(Ce_kalman, E0, Ce50, gammahill); % Convert to BIS
    
    % update variables
    x = x_;

end

u = u_all;