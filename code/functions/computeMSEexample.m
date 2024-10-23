function msetotal = computeMSEexample(v)
% computeMSEkalman computes MSE between kalman BIS prediction and trueDOH.
% Input v is [R;diag(Q)]

    Rmin = v(1);
    Rmax = v(2);
    Q = v(3:6).*eye(4);

    addednoise = 1; % with noise
    % addednoise = 0; % without noise

    [t, trueDOH, ~, ykalman, ~, ~] = simulate_example(Rmin,Rmax,Q,addednoise);

    msetotal = immse(trueDOH,ykalman);
end
