function [t, trueDOH, ymeas, u] = simulate_deltau(whichdist,hasdeltau)
% See documentation for simulate_MPC. Only difference is that we add a term
% to the cost function that punished delta u = u(k) - u(k-1) to prevent
% ringing in the control signal.
% hasdeltau = 0; means no term added to the cost function
% hasdeltau = 1; means added term to cost function

simulation_settings % Get simulation settings

td1 = 5*60; % start time for disturbance
td2 = 15*60; % end time

if hasdeltau == 0 % if no delta u
    alpha_u = 0; % set term in cost function to zero to punished delta u
end

%% Initialize states, u and y

r = get_refconc(N_horizon,E0,Ce50,gammahill);

x = x0; % patient state
x_feedback = x; % state vector seen by mpc

% Save simulation data in vectors
trueCe = zeros(N,1); % True Ce of patient [mg/L]
trueDOH = zeros(N,1); % True DOH of patient [BIS]
ymeas = zeros(N,1); % measurement [BIS]
u_all = zeros(N,1); % Dosing curves [mg/s]

% Simulate
for i = 1:N

    % MPC - find control signal. Set up and solve the QP
    [Hqp, fT, A, b] = setupQP(F, G, alpha_u, umax, r, N_horizon, x_feedback);

    [u_i,Jval] = quadprog(Hqp,fT,A,b,[],[],[],[],[],options);
    u = u_i(1); % Save first dose sample obtained from the optimization
    u_all(i) = u;

    % State update patient PKPD model
    x_ = F*x + G*u;
    Ce = H*x;

    trueCe(i) = Ce; % Save true Ce

    % Compute BIS
    doh_beforedist = computeBIS(Ce, E0, Ce50, gammahill); % Convert to BIS

    x_feedback = x_;

    % Add disturbances
    if whichdist == 1 % disturbance affecting DOH
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
    else
        doh = doh_beforedist;
        if (t(i) >= td1 && t(i) < td2)
            bis_meas = doh + distsize; % add first disturbance to output from PD model
            Ce_dist2 = computeeffectconc(bis_meas,Ce50,E0,gammahill); 
            x_feedback(4) = Ce_dist2; % mpc sees a changed fourth state
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

u = u_all;