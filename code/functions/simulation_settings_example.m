%% Settings file for simulation

rng(12345) % seed for reproducability

%% Sampling

h = 10; % sampling time [s]
tmax = 70*60; % max simulation time [s]
N = tmax/h +1; % Nbr of sampling instances = horizon
umax = 400/60; % max infusion rate, [mg/s]
t = 0:h:tmax;
N_horizon = 10*60/h; % 5 minute horizon for mpc? for now

tstartinf = 10*60; % Time for start of infusion


%% Get patient data from file and Schnider model

dataPatients; % Load patient data

Id = 1;
[V1,V2,V3,Cl1,Cl2,Cl3,gammahill,E0,Emax,Ce50,C50gamma]=get_PKPDparams(Id,Patients);

[V1_v,V2_v,V3_v,Cl1_v,Cl2_v,Cl3_v]=create_perturbed_model(V1,V2,V3,Cl1,Cl2,Cl3); % create patient model

PKPDmodel = get_PKPDmodel(V1,V2,V3,Cl1,Cl2,Cl3); % model run by kalman filter

PKPDmodel_d = c2d(PKPDmodel,h); % discretization of model

% discrete model matrices
F = PKPDmodel_d.A;
Ft = F';
G = PKPDmodel_d.B;
H = PKPDmodel_d.C;
Ht = H';

% F_aug = [F, zeros(size(F,1), 1); -H, 1]; % Augmenting with integral action
% G_aug = [G; 0]; % Augmented input matrix
% H_aug = [H, 0]; % Output matrix stays the same

PKPDpatient = get_PKPDmodel(V1_v,V2_v,V3_v,Cl1_v,Cl2_v,Cl3_v); % true patient model
PKPDpatient_d = c2d(PKPDpatient,h); % discretization of model

% discrete patient system matrices
F_true = PKPDpatient_d.A;
G_true = PKPDpatient_d.B;
H_true = PKPDpatient_d.C;

% F_true_aug = [F_true, zeros(size(F_true,1), 1); -H_true, 1]; % Augmenting with integral action
% G_true_aug = [G_true; 0]; % Augmented input matrix
% H_true_aug = [H_true, 0]; % Output matrix stays the same

n = size(F,1);

% x0 = zeros(n,1); % Initial state

%% Disturbances time

distsize = 20; % optimization fails at 20 (15 works)

t1 = 30*60; % start time
t2 = 35*60; % end time

t3 = 50*60;
t4 = 55*60;

%% Average Hill function values (vanluchene)
gamma_vl = 2.69;
Ce50_vl = 4.92;
% Emax_vl = 87.5;
% Ce50gammamod = 72.6752;
E0_vl = 95.9;

%% reference values

% [BIS_baseline, CeMax, r50] = getlimitrefconc(N, E0, Ce50, gammahill);
Ce_BIS50 = Ce50*((50-E0)/(E0-50-E0))^(1/gamma_vl);

uref = Ce_BIS50/dcgain(PKPDmodel); % initial dose [mg/s] % continuous time?
x0 = -(PKPDmodel.A\(PKPDmodel.B*uref)); % initial state in stationarity

beta_xi = 0.2; % scaling factor for integral action

% Filter reference
reference_unfiltered = [zeros(tstartinf/h,1); Ce_BIS50*ones(tmax/h,1)];

s = tf('s');
T = 30; % time constant [s]
Gf = 1/(s*T + 1);
Hfilter = c2d(Gf,h);

[num, den] = tfdata(Hfilter, 'v');

reference = filter(num, den, reference_unfiltered); % filtered reference

% CeMax = computeeffectconc(40,Ce50,E0,gammahill);

%% quadprog settings

options = optimoptions(@quadprog,'Display','off');

alpha_u = 0.02; % in cost function to regulate punishment of u(k) - u(k-1)

%% Noise
structnoise = load('real_noise.mat');
noise_profile = repmat(structnoise.noise,1,3);

% noise
inoise = randi(3001,1,1);
noise = noise_profile(inoise:end);

maxdelay = 2*60; % If SQI = 0, how long delay?

%% Filter settings
s = tf('s');
T = 8; % time constant
Gf = 1/(s*T + 1)^2;
Hfilter = c2d(Gf,h);
sys_filter = ss(Hfilter);
F_filt = sys_filter.A;
G_filt = sys_filter.B;
H_filt = sys_filter.C;


