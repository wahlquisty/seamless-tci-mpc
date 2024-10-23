function [PKPDmodel]=get_PKPDmodel(V1,V2,V3,Cl1,Cl2,Cl3)
%% This function is required to set the parameters of simPropofol.slx
%Input: PK parameters V1,V2,V3,Cl1,Cl2,Cl3 [L] and [L/min]
%Output: PK and PD models for Propofol, gamma, E0, Emax, C50*gamma
% Four-compartment linear PKPD model

%% TRANSFER RATES

%Transfer Rate of Drug
%N.B. kij transfer rate from i to j
k10	= Cl1 / V1 / 60;                                                %[1/s]
k12	= Cl2 / V1 / 60;                                                %[1/s]
k13	= Cl3 / V1 / 60;                                                %[1/s]
k21	= Cl2 / V2 / 60;                                                %[1/s]
k31	= Cl3 / V3 / 60;                                                %[1/s]
k41	= 0.459 / 60;                                                   %[1/s]
ke0  = k41;                                                         %[1/s]

%% STATE SPACE LINEAR PHARMACOKINETIC-PHARMACODYNAMIK MODEL (PKPD)

%State Matrix
Assp =[-(k10+k12+k13)  k21   k31    0;
        k12           -k21     0    0;
        k13              0  -k31    0;
        ke0              0     0 -ke0];

%Input Matrix
% Bssp = [ 1; 0; 0; 0 ];
Bssp = [ 1/V1; 0; 0; 0 ];

%Output Matrix
% Cssp = [ 0 0 0 1] / V1;
Cssp = [ 0 0 0 1];

%Feedthrough Matrix
Dssp = 0;

PKPDmodel = ss( Assp , Bssp , Cssp , Dssp );



end