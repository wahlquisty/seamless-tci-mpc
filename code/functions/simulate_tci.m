function [t, trueDOH, ymeas, u] = simulate_tci(dist)
%SIMULATE_TCI Simulate Target-Controlled Infusion (TCI) under stationary conditions.
%
%   [t, trueDOH, ymeas, u] = SIMULATE_TCI(dist) in stationarity (BIS = 50)
%   Simulates two types of disturbances
%   - If dist = 1, the disturbance affects the Degree of Hypnosis (DOH),
%     which affects to the fourth state variable.
%   - If dist = 2, the disturbance affects the measurement, introducing 
%     noise into the observations.
%
%   Outputs:
%   - t: Time vector representing the simulation duration.
%   - trueDOH: True DoH values after simulation.
%   - ymeas: Measurement [BIS]
%   - u: Control input values applied during the simulation, computed for
%   the TCI

%% Setup, initialize states, u and y

simulation_settings; % Simulation settings from separate file

r = get_refconc(N,E0,Ce50,gammahill); % get reference concentration for BIS 50

x = x0; % patient state

% Save simulation data in vectors
trueCe = zeros(N,1); % True Ce of patient [mg/L]
trueDOH = zeros(N,1); % True DOH of patient [BIS]
ymeas = zeros(N,1); % measurement [BIS]

% TCI: find control signal u. Set up and solve the QP
[Hqp, fT, A, b] = setupQP(F, G, alpha_u, umax, r, N, x0); % Setup QP
[u,Jval] = quadprog(Hqp,fT,A,b,[],[],[],[],[],options); % Solve the QP problem

% Simulate
for i = 1:N 

    % State update patient PKPD model
    x_ = F*x + G*u(i);
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
        else
            bis_meas = doh;
        end
    end

    % Save DoH and measurement
    trueDOH(i) = doh;
    ymeas(i) = bis_meas;

    % update variables
    x = x_;

end

end