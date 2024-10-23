%% Settings file for simulations

rng(123456) % seed for reproducability

%% Sampling

h = 10; % sampling time [s]
tmax = 25*60; % max simulation time [s]
N = tmax/h +1; % Nbr of sampling instances, also horizon for TCI
umax = 400/60; % max infusion rate, [mg/s]
t = 0:h:tmax; % Time vector
N_horizon = 10*60/h; % 10 minute horizon for MPC

%% Get patient data from file and Schnider model

dataPatients; % Load patient data from file

Id = 1;
[V1,V2,V3,Cl1,Cl2,Cl3,gammahill,E0,Emax,Ce50,C50gamma]=get_PKPDparams(Id,Patients); % Get patient data
PKPDmodel = get_PKPDmodel(V1,V2,V3,Cl1,Cl2,Cl3);

PKPDmodel_d = c2d(PKPDmodel,h); % discretize patient model

% discrete model matrices
F = PKPDmodel_d.A;
Ft = F';
G = PKPDmodel_d.B;
H = PKPDmodel_d.C;
Ht = H';

n = size(F,1); % Number of states

%% Time for disturbance

distsize = 20; % Size of disturbance [BIS]

td1 = 5*60; % start time of disturbance
td2 = 10*60; % end time of disturbance

%% Reference values

r = get_refconc(N,E0,Ce50,gammahill);

uref = r(1)/dcgain(PKPDmodel); % initial dose [mg/s] for stationarity
x0 = -(PKPDmodel.A\(PKPDmodel.B*uref)); % initial state in stationarity

%% Optimization, quadprog settings

options = optimoptions(@quadprog,'Display','off');

alpha_u = 0.01; % in cost function to regulate punishment of u(k) - u(k-1)

