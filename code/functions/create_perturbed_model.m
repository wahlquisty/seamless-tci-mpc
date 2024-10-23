function [V1_v,V2_v,V3_v,CL1_v,CL2_v,CL3_v]=create_perturbed_model(V1,V2,V3,CL1,CL2,CL3)
% Creates perturbed 3-comp PK model with the Schnider random effects.

% Coefficient of variation dei parametri [%]
CV_V1 = 4.04;
CV_V2 = 1;
CV_V3 = 14.35;
CV_CL1 = 10.05;
CV_CL2 = 1;
CV_CL3 = 11.79;


% Calcolo le deviazioni standard dell'esponenziale
sigma_V1 = sqrt(log((CV_V1^2/100^2)+1));
sigma_V2 = sqrt(log((CV_V2^2/100^2)+1));
sigma_V3 = sqrt(log((CV_V3^2/100^2)+1));
sigma_CL1 = sqrt(log((CV_CL1^2/100^2)+1));
sigma_CL2 = sqrt(log((CV_CL2^2/100^2)+1));
sigma_CL3 = sqrt(log((CV_CL3^2/100^2)+1));


% Genero i parametri random
V1_v = V1*exp(normrnd(0,sigma_V1,1));
V2_v = V2*exp(normrnd(0,sigma_V2,1));
V3_v = V3*exp(normrnd(0,sigma_V3,1));
CL1_v = CL1*exp(normrnd(0,sigma_CL1,1));
CL2_v = CL2*exp(normrnd(0,sigma_CL2,1));
CL3_v = CL3*exp(normrnd(0,sigma_CL3,1));

end