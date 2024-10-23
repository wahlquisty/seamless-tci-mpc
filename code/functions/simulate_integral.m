function [t, trueDOH, monitorBIS, ykalman, u] = simulate_integral(Rmin,Rmax,Q,hasintegralaction)

simulation_settings_example; % get simulation settings

if hasintegralaction == 0 % no integral action
    beta_xi = 0;
end

%% SQI stable
sqi = 100*ones(N,1);

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

   doh = doh_beforedist;
   bis_meas = doh;

    % Save DoH and measurement
    trueDOH(i) = doh;
    ymeas(i) = bis_meas;

    ymonitor = ymeas(i);

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