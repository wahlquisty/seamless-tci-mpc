function [t, trueDOH, monitorBIS, ykalman, sqi, u] = simulate_example(Rmin,Rmax,Q,addednoise)
%SIMULATE_EXAMPLE Simulation example of indcution + maintenance with two
%disturbances
%
%   [t, trueDOH, ymeas, u] = SIMULATE_EXAMPLE(Rmin,Rmax,Q,addednoise)
%   Inputs:
%   - Rmin: smallest value of R, if SQI = 100, R = Rmin. Measurement
%   covariance in Kalman filter
%   - Rmax: largest value of R, if SQI = 0, R = Rmax
%   - Q: quadratic matrix, state covariance in Kalman filter
%   - addednoise: If addednoise = 0, no noise. addednoise= 1, noise profile
%   added to simulation.

%   Outputs:
%   - t: Time vector representing the simulation duration.
%   - trueDOH: True DoH values after simulation.
%   - ymeas: Measurement [BIS]
%   - ykalman: Kalman filter predictions [BIS]
%   - u: Control input values applied during the simulation, computed with
%   the MPC

simulation_settings_example; % get simulation settings

%% SQI drop for second disturbance (1 sample after disturbance)

SQIstep = 50 ; % Step size for SQI (100->SQIstep)
sqi = 100*ones(N,1);
sqi(t3/h+2:t4/h+1) = SQIstep; % index is off by one

%% Initialize states, u and y

x0 = zeros(n,1);
x = x0; % True state of system
xi = 0; % integral_state
xi_all = zeros(N,1);

x_kalman = x0; %zeros(4,1); % State vector Kalman filter
x_kalman_m1 = x_kalman;
P0 = eye(n); % initial covariance matrix. FIXME: how to initialize?
P_m1 = P0; % P_m1 = Pn,n-1. Covariance matrix
I4 = eye(n); % identity of size 4

x_filt = zeros(size(G_filt)); % filter state

trueCe = zeros(N,1);
trueDOH = zeros(N,1);
ymeas = zeros(N,1);
ykalman = zeros(N,1); % predicted BIS by Kalman filter
u_all = zeros(N,1);
monitorBIS = zeros(N,1);

options = optimoptions(@quadprog,'Display','off');

% Simulate
for i = 1:N

    % i

    if i < N-N_horizon
        r = reference(i:i+N_horizon-1);
    else
        r = reference(i)*ones(N_horizon,1);
    end

    % MPC - find control signal. Set up and solve the QP
    [Hqp, fT, A, b] = setupQP(F, G, alpha_u, umax, r, N_horizon, x_kalman_m1);

    [u_i,Jval] = quadprog(Hqp,fT,A,b,[],[],[],[],[],options);
    u_qp = u_i(1); % Save first dose sample obtained from the optimization

    u0 = beta_xi*xi;
    u = u_qp + u0; % add integral error

    % insert anti-windup
    if u > umax
        u = umax;
        xi = min(xi,(umax-u_qp)/beta_xi);
    elseif u < 0 
        u = 0;
        xi = max(xi,(0-u_qp)/beta_xi);
    end

    u = max(0,u);
    u = min(umax,u);

    u_all(i) = u;

    
    % State update patient PKPD model
    x_ = F_true*x + G_true*u;
    Ce = H_true*x;

    trueCe(i) = Ce; % Save true Ce

    % Compute BIS
    doh_beforedist = computeBIS(Ce, E0, Ce50, gammahill); % Convert to BIS

   % Add disturbances
     % disturbance affecting DOH
        if t(i) == t1
            doh = doh_beforedist + distsize; % add first disturbance to output from PD model
            Ce_dist = computeeffectconc(doh,Ce50,E0,gammahill);
            x_(4) = Ce_dist;  % patient new affected fourth state
        elseif t(i) == t2
            doh = doh_beforedist - distsize; % add first disturbance to output from PD model
            Ce_dist = computeeffectconc(doh,Ce50,E0,gammahill);
            x_(4) = Ce_dist;   % patient new affected fourth state
        else
            doh = doh_beforedist;
        end
     % disturbance affecting measurement
        % doh = doh_beforedist;
        if t(i) >= t3 && t(i) < t4
            bis_meas = doh + distsize; % add first disturbance to output from PD model
        else
            bis_meas = doh;
        end

    % Save DoH and measurement
    trueDOH(i) = doh;
    ymeas(i) = bis_meas;


        % delay measurement when signal quality is low
        delay = (1-sqi(i)/100)*maxdelay; % compute delay from SQI. Create
        if (t(i) >= t3 + 1) && (t(i) < t3+1+delay) % signal is fixed for some time that depends on the sqi
            ymonitor = monitorBIS(t3/h+1);
        else
            ymonitor = ymeas(i-delay/h); %trueDOH(i-delay); Ska det stå ymeas här?
        end

    if addednoise == 1
        ymonitor = ymonitor + noise(i); % add noise from noise profile

        % filter signal
        x_filt_ = F_filt*x_filt + G_filt*ymonitor;
        ymonitor = H_filt*x_filt_;
        x_filt = x_filt_;
    end

    ymonitor = max(0,ymonitor);
    ymonitor = min(ymonitor,E0);
    monitorBIS(i) = ymonitor;

    % Kalman filter
    Ce_meas = computeeffectconc(ymonitor,Ce50_vl,E0,gamma_vl); % feedback concentration inverted through vanluchene hill fct FIXME: Required to have nonlinear MPC? Ce50_vl not working

    trueCe(i) = Ce_meas;
    
    R = Rmin + (Rmax - Rmin)*(1-sqi(i)/100); % linear dependence between sqi and R from Rmax to Rmin
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
    Ce_kalman = H*x_kalman_m1;
    ykalman(i) = computeBIS(Ce_kalman, E0, Ce50_vl, gamma_vl); % E0 can be readily estimated

    % update variables
    x = x_;
    xi = xi + (r(1) - Ce); % integral error, FIXME, se över

    xi_all(i) = xi;
end

u = u_all;