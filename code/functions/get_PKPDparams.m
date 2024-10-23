function [V1,V2,V3,Cl1,Cl2,Cl3,gamma,E0,Emax,C50,C50gamma]=get_PKPDparams(Id,Patients)
%% This function returns PKPD parameters (three-compartment + hill function)
%Input: Id of the patient, Patients Database
%Output: PK model parameters for Propofol, V1, V2, V3, CL1, Cl2, Cl3 [L] and [L/min]
% Schnider model for PKPD model
% Three-compartment linear PK model

%% DATA

%loading data from external file
Age    = Patients( Id , 2 );
Height = Patients( Id , 3 );
Weight = Patients( Id , 4 );
Gender = Patients( Id , 5 );


%% MODEL PARAMETERS

if Gender == 1
    %lean body mass (James Formula for Men)
    lbm = 1.1*Weight - 128*( Weight / Height )^2;
else % 2 per donne
    %lean body mass (James Formula for Women)
    lbm = 1.07*Weight - 148*( Weight / Height )^2;
end

%Volume of the compartments (Schneider Model for Propofol)
V1 = 4.27;                                                          %[L]
V2 = 18.9 - 0.391*( Age - 53 );                                     %[L]
V3 = 238;                                                          %[L]

%Clearance of compartments
Cl1 = 1.89 + 0.0456*( Weight - 77 ) - 0.0681*( lbm - 59 ) +...      
      0.0264*( Height - 177 );                                      %[L/s]
Cl2 = 1.29 - 0.024*( Age - 53 );                                    %[L/s]
Cl3 = 0.836;                                                        %[L/s]

% Should it not be [L/min] for clearance?

%% Hill function parameters

C50 = Patients( Id, 6 );
gamma = Patients( Id, 7 );
E0 = Patients( Id, 8);
Emax = Patients( Id, 9);
C50gamma=C50^gamma;

end