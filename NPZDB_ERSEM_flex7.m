% NPZDB_ERSEM_flex7.m
% driver m-file for model ERSEM ecosystem model based on 
% Blackford et al. (2004) and Butenshon et al. 2016.
% Tatsuro Tanioka 08/29/17
% kernel program is NPZDB_ERSEM_flex_func7.m 
% 4 Phytoplankton types and 3 zooplankton types
% (nanoglagellates, microz, and meso zoo)
% C:N:P and max. assimilation rate of Mesozooplankton variable as a function of age
% Copepod has prey switching capability based on prey biomass and prey food
% quality
% ver 1: 08/29/17 Modified all the paramters to B16, added Fasham Preference Model
% ver 2: 08/30/17 C:N:P of mesozooplankton is a function of age
% ver 3: 08/30/17 max. growth of mesozooplankton is a function of age 
% ver 4: 08/30/17 Included switch to activate age-dependent copepod model
% ver 5: 08/31/17 Calibrated copepod age model with Meunier16 (actual data
% regression)

clear all
close all

% Basic setup
number_years = 5; % 3 years of experiment
run_name = "170831a"; % 

% Experimental Variables
global can_threshold % [1.] threshold value for cannibalism 1 = no threshold
global Pf_var % [1.] Switch to use flexible Pf (=1.0) for fixed Pf (=0.0)
global sigma_pf % [1.] standard deviation of the grazing spectrum matrix
global COPE_AGE_switch % [1.] Switch to use copepod age model (1 = true, 0 =false)
global COPE_AGE_cpslope % Slope of C:P against age (d-1)
global COPE_AGE_cpintercept % y-intercept of C:P against age (molC/molP)
global COPE_AGE_npslope % Slope of N:P against age (d-1)
global COPE_AGE_npintercept % y-intercept of N:P against age (molN/molP)
global COPE_AGE_rmaxslope % Slope of rmax against age (d-1)
global COPE_AGE_rmaxintercept % y-intercept of rmax against age (d-1 d-1)

can_threshold = 0.25 ; % threshold for feeding itself, 1 = no threshold
Pf_var = 1; % 0 = Fixed Pf (B16), 1 = Variable Pf (FQ+Fasham), 2= Variable Pf (Fasham)
sigma_pf = 0.05;
COPE_AGE_switch = 1; % Copepod age switch, 1= on, 0 = off
COPE_AGE_cpslope = 5.77; % From Meunier16
COPE_AGE_cpintercept = 52.6; % From Meunier16
COPE_AGE_npslope = 1.94; % From Meunier16
COPE_AGE_npintercept = 7.15; % From Meunier16
COPE_AGE_rmaxslope = -0.06; % From Meunier16, set to 0 if no age dependence
COPE_AGE_rmaxintercept = 1.36; % Calibrated, set to 1.0 if no age dependence

cope_age = [1:1:14]; % copepod age range
% -------------------
% Environment parameters
global A  %[/day]     CHEMOSTAT DILUTION RATE 
global Nin_N  %[mmol m^-3]     INPUT NO3 CONCENTRATION (N limited, N/P co-limitation, P limited)
global Nin_A  % [mmol m^-3]    INPUT NH4 CONCENTRATION 
global Nin_P  %[mmol m^-3]     INPUT P CONCENTRATION
global temperature %[C] Temperature
global PAR %[W m-2] Insolation
global O_rel % [mmol m^-3] % relative oxygen saturation

A = 0.1; % Chemostat dilution rate (Mixing rate) (1/day) 
Nin_N = 16; 
Nin_A = 1;
Nin_P = 1; 
temperature = 10.*ones(1,number_years.*365);
PAR = 200;
O_rel = 1.0; % Relative Oxygen Level
Bio_0 = 20.0; % Initial biomass [mg C m^-3]

% Techinal setup
option_safety = 1; %  0 = safety_off, 1 = safety_on 
safety = 1.0e-2; %[mgC m^-3] Default Negativity control factor = 0.01 (Butenschön et al., 2016)
num_ex = length(cope_age); % Number of experiments


% Ecosystem setting
global n_P % [unitless] Number of phytoplankton types
global n_Z % [unitless] Number of plankton types
global n_var % [unitless] Total Number of variables
n_P = 4; % # of phytoplankton types
n_Z = 3; % # of zootoplankton types
n_var = 48; % # of variables

% General Indenpendent Parameters

global Q10    %[unitless]        Q10 value 
global r_nitr  %[/day]       Max. Nitrification rate 
global h_N_nitr %[mmol m^-3] Cubic MMenten constant for nitrification 
global Q_N_RED %[mmol N /mgC] Redfield N:C Ratio
global Q_P_RED %[mmol P /mgC] Redfield P:C Ratio

Q10 = 2.0;   %[unitless]          Q10 value (Blackford, 2004, for diatoms and microzooplankton)
r_nitr  = 0.5; %[/day]       Max. Nitrification rate (Butenschön et al., 2016)
h_N_nitr = 0.5; %[mmol m^-3] Cubic MMenten constant for nitrification (Butenschön et al., 2016)
Q_N_RED = 0.0126;%[mmol N /mgC] Redfield N:C Ratio
Q_P_RED = 0.786e-3;%[mmol P /mgC] Redfield P:C Ratio

% Safety switch
% Strictly speaking, this violates conservation of mass but gives the
% system more stability. Set epsilon = 0.0 for mass balance check.
global epsilon %[mgC m^-3] Negativity control factor 
global epsilon_p %[mgC m^-3] Negativity control factor * (P:C)redfield (Butenschön et al., 2016)
global epsilon_n %[mgC m^-3] Negativity control factor * (N:C)redfield (Butenschön et al., 20162
global epsilon_chl %[mgC m^-3] Negativity control factor * (chl:C)min (Butenschön et al., 2016)
epsilon = double(option_safety>0).*safety; 
epsilon_p = epsilon*Q_P_RED; %[mmolP m^-3] Negativity control factor for P
epsilon_n = epsilon*Q_N_RED; %[mmolP m^-3] Negativity control factor for N
epsilon_chl = epsilon*0.0067; %[mgchl m^-3] Negativity control factor for Chl
% 
% Phytoplankton Parameters (diatom)
global R_ASS  %[/day]       Assimilation rate at 10C 
global R_RESTR % [/day]     Basal respiration at 10C 
global R_ACTR % [/day]      Activity respiration at 10C 
global R_NLUX % [/day]      rate of nutrient luxury uptake towards max. quota
global R_LYSIS %[/day]      Minimal lysis rate 
global P_EX   % [unitless] Exudation under nutrient stress 
global QMIN_N_PHYTO %[mmolN / mgC]    MIN. NQUOTA 
global QMAX_N_PHYTO  %[mmolN / mgC]    MAX. NQUOTA 
global QMIN_P_PHYTO  %[mmolP / mgC]    MIN. PQUOTA 
global QMAX_P_PHYTO %[mmolP / mgC]    MAX. PQUOTA 
global ALPHA  %[W/m^2/day] Initial slope of P-I curve 
global BETA   %[W/m^2/day] Light inhibition of photosynthesis 
global THETA_MAX %[mgChl / mg C] Max. Chl/C cell ratio  
global THETA_MIN %[mgChl / umg C] Min. Chl/C cell ratio  
global A_N  %[/ (mg C day)]    affinity for N uptake 
global A_A  %[/ (mg C day)]    affinity for amm uptake 
global A_P  %[/ (mg C day)]    affinity for P uptake 
global RSED_PHYTO %[day-1] max. subsiding velocity for phytoplankton under stress
global ST  % [unitless] Nutrient stress threshold 
global SIZE_PHYTO % [um] typical diameter of phytoplankton [Diatoms Flagellates Pico-plankton Dino-flagellates]

% Grazer Parameters
global R_ASS_ZOO      %[/day]  Assimilation rate at 10C
global H %[mgC m^-3] Food conc. where relative uptake is 0.5
global AE % [1] Assimilation efficiency
global P_DOM % [1] Fraction of excretion going to DOM
global R_RESTR_ZOO % [/day] Basal respiration at 10C 
global H_OXMORT % [mmol m^-3] Ox. saturation where respiration is 0.5;
global EU % [1] Excreted fraction of uptake
global R_MORTOX % [1/day] Mortality rate due to ox. limitation
global R_MORT % [1/day] Temp. independent mortality rate
global QMAX_N_ZOO  %[mmolN / mgC]  MAX N Quota for zooplankton
global QMAX_P_ZOO  %[mmolP / mgC]  MAX P Quota for zooplankton
global Z_MINFOOD  %[mgC m^-3]    Lower Foodthreshold for feeding
global cn %[unitless] relaxation rate of N-excretion 
global cp %[unitless] relaxation rate of P-excretion 
global Pf %[unitless] prey availability matrix
global q_Rexcr_meso % [unitless] Excreted fraction of POM uptake by mesozooplankton
global rho_TER % [1.] max.GGE(N)/max.GGE(P)
global COPE_AGE % Mean age of copepod (days)

% Detrital parameters
global Ra %[1/day] Fraction of detritus breaking down
global rfd_cn % [1] Redfield C/N ratio
global r_remin % [1/day] Remin rate
global V_SINK %[1/day] Sinking velocity [sPOM mPOM lPOM]
global q_dis % [1.] fraction going into labile DOM pool over semi-labile pool

% Bacteria Parameters
global r_ass_bac % [1/day] assimilation rate at 10C
global h_ox_bac % [1.] half saturation oxygen limitation
global r_resp_bac % [1.] respired fraction of uptake
global r_respox_bac % [1.] respired fraction of uptake at low oxygen conc.
global r_basal_bac % [1/day] basal respiration at 10C
global r_mort_bac % [1/day] mortality rate
global QMAX_N_BAC % [mmol N /mgC] Max cell-quotum N
global QMAX_P_BAC % [mmol P / mgC] Max cell-quotum P
global hp_bac % [mmol P /m^3] Michaelis constant for P uptake
global hn_bac % [mmol N / m^3] Michaelis constant for N uptake

% Set parameter values

% Phytoplankton parameters [diatoms flagellates pico-phytoplankton dinoflagellates]
R_ASS = [1.375 1.625 2.0 1.125]; %[/day]       Assimilation rate at 10C (B16)
R_RESTR = [0.04 0.04 0.045 0.035]; % [/day]     Basal respiration at 10C (B16)
R_ACTR = [0.1 0.25 0.25 0.25]; % [/day]      Activity respiration at 10C 
R_LYSIS = [0.05 0.05 0.05 0.055]; %[/day] Minimal lysis rate (B16)
R_NLUX = [1.0 1.0 1.0 1.0]; % [/day]      rate of nutrient luxury uptake towards max. quota (tuned)
P_EX = [0.05 0.2 0.2 0.05];  % [unitless] Exudation under nutrient stress 
 QMIN_N_PHYTO = [4.2e-3 5.0e-3 6.0e-3 4.2e-3]; %[mmolN / mgC]    MIN. NQUOTA (B16)
 QMAX_N_PHYTO = [4.2e-3.*1.075 5.0e-3.*1.075 6.0e-3.*1.05 4.2e-3.*1.05]; %[mmolN / mgC]    MAX. NQUOTA (B16)
 QMIN_P_PHYTO = [1.0e-4 2.25e-4 3.5e-4 1.0e-4]; %[mmolP / mgC]    MIN. PQUOTA (B16)
 QMAX_P_PHYTO = [1.0e-4.*2 2.25e-4.*2 3.5e-4.*1.5 1.0e-4.*2.7]; %[mmolP / mgC]    MAX. PQUOTA (B16)

ALPHA = [4.00 5.00 6.00 3.00];  %[W/m^2/day] Initial slope of P-I curve (B16)
BETA = [7.0e-2 1.0e-1 1.2e-1 6.0e-2];  %[W/m^2/day] Light inhibition of photosynthesis (B16)
THETA_MAX = [6.0e-2 2.5e-2 1.5e-2 4.5e-2]; %[mgChl / mgC] Max. Chl/C cell ratio (B16)
THETA_MIN = [0.0067 0.0067 0.0067 0.0067]; %[mgChl / mgC] Min. Chl/C cell ratio  
A_N = [2.5e-3 4.0e-3 6.0e-3 2.0e-3];  %[/ (mg C day)]    affinity for N uptake (B16)
A_A = [2.5e-3 4.0e-3 7.0e-3 2.0e-3]; %[/ (mg C day)]    affinity for amm uptake (B16)
A_P = [3.0e-3 4.0e-3 6.0e-3 2.0e-3];  %[/ (mg C day)]    affinity for P uptake (B16)
RSED_PHYTO = [5.0/100 0.0 0.0 5.0/100]; %[m day-1/m]] (assuming 100 m box)max. subsiding velocity 
ST = [0.70 0.7 0.7 0.7]; % [unitless] Nutrient stress threshold (B16)
SIZE_PHYTO = [100.0 10.0 1.0 50.0]; % [um] typical diameter of phytoplankton [Diatoms Flagellates Pico-plankton Dino-flagellates]

% Zooplankton parameters [heteroflagellates microzoo mesozoo] 
H_OXMORT = [7.81 7.81 7.81];% [mmol m^-3] Ox. saturation limitaion constant (B16);
H = [28.0 32.0 36.0];%[mgC m^-3] Food conc. where relative uptake is 0.5 (B16)
Z_MINFOOD = [12.0 12.0 12.0];%[mgC m^-3]    Lower Foodthreshold for feeding (B16)
P_DOM = [0.5 0.5 0.5]; % Fraction of excretion going to DOM (B16)
q_Rexcr_meso = 0.9; % [1.] Excreted fraction of POM uptake by mesozooplankton (B16)
EU = [0.5 0.5 0.5];% [1] Excreted fraction of uptake (B16)
AE = [0.4 0.5 0.6]; % [1] Assimilation efficiency (B16)
QMAX_N_ZOO = [0.0167 0.0167 NaN];%[mmolN / mgC]  MAX N Quota for zooplankton (B16)
QMAX_P_ZOO = [0.001 0.001 NaN]; %[mmolP / mgC]  MAX P Quota for zooplankton (B16)
R_MORTOX = [0.3 0.25 0.2];% [1/day] Mortality rate due to ox. limitation (B16)
R_MORT = [0.05 0.05 0.05];% [1/day] Temp. independent mortality rate (B16)
R_RESTR_ZOO = [0.025 0.02 0.015];% [/day] Basal respiration at 10C (B16)
cn = [0.5 0.5 NaN]; %[1/day] relaxation rate in N-excretion (B16)
cp = [0.5 0.5 NaN]; %[1/day] relaxation rate in P-excretion (B16)
R_ASS_ZOO = [1.5 1.25 1.0]; %[/day]  Assimilation rate at 10C (B16)

rho_TER = [1.0 1.0 1.0]; % GGE.N/GGE.P [sz uz mz]
% Prey availability matrix
    % From B16
Pf = [0.45 0.1 0.0; 0.25 0.15 0.0; 0.15 0.15 0.05; 0.0 0.15 0.15; 0.0 0.10 0.15; 0.15 0.2 0.05; 0.0 0.15 0.25; 0.0 0.0 0.25; 0.0 0.0 0.1];

% Detrital parameters
Ra = 0.01; %[1/day] Fraction of detritus breaking down
rfd_cn = 6.625;% [1] Redfield C/N ratio
r_remin = 0.05; % [1/day] Remin rate
V_SINK = [1./100 5.0./100 10.0/100];  %[1/day] Sinking velocity [sPOM mPOM lPOM] B16, assuming 100m box length
q_dis = 0.5; % [1.] fraction going into labile DOM pool over semi-labile pool (tuned)
% Bacteria Parameters
r_ass_bac = 4.0;% [1/day] assimilation rate at 10C
h_ox_bac = 0.3125;% [1.] half saturation oxygen limitation
r_resp_bac = 0.4;% [1.] respired fraction of uptake
r_respox_bac = 0.2;% [1.] respired fraction of uptake at low oxygen conc.
r_basal_bac = 0.1;% [1/day] basal respiration at 10C
r_mort_bac = 0.05; % [1/day] mortality rate
QMAX_N_BAC = 0.0126;% [mmol N /mgC] Max cell-quotum N
QMAX_P_BAC = 0.000786;% [mmol P / mgC] Max cell-quotum P
hp_bac = 0.1;% [mmol P /m^3] Michaelis constant for P uptake
hn_bac = 0.5;% [mmol N / m^3] Michaelis constant for N uptake

% Initial values 
Nn_0 = 1.0; % Initial Nn [mmol m^-3 N]
Na_0 = 0.1; % Initial Na [mmol m^-3 N]
Np_0 = 0.1; % Initial Np [mmol m^-3 P]
Pc_0   = Bio_0; % Initial Pc [mg C m^-3]
Pn_0   = Pc_0.*Q_N_RED; % Initial Pn [mmol N m^-3]
Pp_0   = Pc_0.*Q_P_RED; % Initial Pp 
Pchl_0   = Pc_0.*THETA_MIN(1); % Initial Pchl [mg C m^-3]
P2c_0   = Bio_0; % Initial P2c [mg C m^-3]
P2n_0   = P2c_0.*Q_N_RED; % Initial P2n [mmol N m^-3]
P2p_0   = P2c_0.*Q_P_RED; % Initial P2p 
P2chl_0   = P2c_0.*THETA_MIN(2); % Initial P2chl [mg C m^-3]
P3c_0   = Bio_0; % Initial P2c [mg C m^-3]
P3n_0   = P3c_0.*Q_N_RED; % Initial P2n [mmol N m^-3]
P3p_0   = P3c_0.*Q_P_RED; % Initial P2p 
P3chl_0   = P3c_0.*THETA_MIN(3); % Initial P2chl [mg C m^-3]
P4c_0   = Bio_0; % Initial P2c [mg C m^-3]
P4n_0   = P4c_0.*Q_N_RED; % Initial P2n [mmol N m^-3]
P4p_0   = P4c_0.*Q_P_RED; % Initial P2p 
P4chl_0   = P4c_0.*THETA_MIN(4); % Initial P2chl [mg C m^-3]
Zc_0  = Bio_0; % Initial Zc [mg C m^-3]
Zn_0   = Zc_0.*Q_N_RED; % Initial Zn [mmol N m^-3]
Zp_0   = Zc_0.*Q_P_RED; % Initial Zp [mmol P m^-3]
Z2c_0  = Bio_0; % Initial Z2c [mg C m^-3]
Z2n_0   = Z2c_0.*Q_N_RED; % Initial Z2n [mmol N m^-3]
Z2p_0   = Z2c_0.*Q_P_RED; % Initial Z2p [mmol P m^-3]
Z3c_0  = Bio_0; % Initial Z3c [mg C m^-3]
% Z3n_0   = Z3c_0.*Q_N_RED; % Initial Z3n [mmol N m^-3]
% Z3p_0   = Z3c_0.*Q_P_RED; % Initial Z3p [mmol P m^-3]
Z3n_0   = Z3c_0.*(COPE_AGE_npslope.*cope_age + COPE_AGE_npintercept).*1.0./(COPE_AGE_cpslope.*cope_age + COPE_AGE_cpintercept)./12; % Initial Z3n [mmol N m^-3]
Z3p_0   = Z3c_0.*1.0./(COPE_AGE_cpslope.*cope_age + COPE_AGE_cpintercept)./12; % Initial Z3p [mmol P m^-3]
Bc_0 = 20.0; % Initial Bc [mg C m^-3]
Bn_0 = Bc_0*Q_N_RED; % Initial Bn [mmol N m^-3]
Bp_0 = Bc_0*Q_P_RED; % Initial Bp [mmol  m^-3]
sPOC_0 = 1.0e-24; % Initial sPOC [mg C m^-3]
sPON_0 = sPOC_0*Q_N_RED; % Initial sPON [mmol N m^-3]
sPOP_0 = sPOC_0*Q_P_RED; % Initial sPOP [mmol  m^-3]
mPOC_0 = sPOC_0; % Initial mPOC [mg C m^-3]
mPON_0 = sPON_0; % Initial mPON [mmol N m^-3]
mPOP_0 = sPOP_0; % Initial mPOP [mmol  m^-3]
lPOC_0 = sPOC_0; % Initial lPOC [mg C m^-3]
lPON_0 = sPON_0; % Initial lPON [mmol N m^-3]
lPOP_0 = sPOP_0; % Initial lPOP [mmol  m^-3]
DOC_0 = 80.0; % Initial DOC [mg C m^-3] Polimene16
DON_0 = DOC_0*Q_N_RED; % Initial DON [mmol N m^-3] Polimene16
DOP_0 = DOC_0*Q_P_RED; % Initial DOP [mmol  m^-3] Polimene16
slDOC_0 = 1.0e-24; % Initial slDOC [mg C m^-3] Polimene16
slDON_0 = slDOC_0*Q_N_RED; % Initial slDON [mmol N m^-3] Polimene16
slDOP_0 = slDOC_0*Q_P_RED; % Initial slDOP [mmol  m^-3] Polimene16

% ------------- Calculating Numerical Solution -------------------
length_days = number_years.*365;
tspan = [1:length_days]; % Time (day), use a fixed time scale
n = length(tspan); % pointer to last time step
X = zeros(n,n_var); % Initialize array
result = [];
result_n = [];

% Running model 
for i = 1:num_ex
    COPE_AGE = cope_age(i)
    x0 = [Nin_N Nin_A Nin_P Pc_0 Pchl_0 Pn_0 Pp_0 P2c_0 P2chl_0 P2n_0 P2p_0 P3c_0 P3chl_0 P3n_0 P3p_0 P4c_0 P4chl_0 P4n_0 P4p_0 Zc_0 Zn_0 Zp_0 Z2c_0 Z2n_0 Z2p_0 Z3c_0 Bc_0 Bn_0 Bp_0 sPOC_0 sPON_0 sPOP_0 mPOC_0 mPON_0 mPOP_0 lPOC_0 lPON_0 lPOP_0 DOC_0 DON_0 DOP_0 Z3n_0(i) Z3p_0(i) slDOC_0 slDON_0 slDOP_0 0.0 0.0 0.0 0.0]'; % Starrtng values for dependent variables X
    [T,X] = ode15s('NPZDB_ERSEM_flex_func7',tspan,x0);
    result_n = cat(3,result_n,X);
end    
result = permute(result_n,[3,1,2]); % permutate result
X = result; % X = result matrix
X_ave = squeeze(nanmean(X(:,length_days-365:length_days,:),2)); % Avg. of the last year

% % Calculating Pf matrix for each mean age
Pf_result = [];
for i = 1:num_ex;
    Pf_n = Pf_flex_func(squeeze(X_ave(i,:)),sigma_pf,Pf_var);
    Pf_result = cat(3,Pf_result,Pf_n);
end
Pf_result = permute(Pf_result,[3,1,2]);
Pf1 = squeeze(Pf_result(1,:,:)); % Pf for mean age 1 
Pf10 = squeeze(Pf_result(10,:,:)); % Pf for mean age 6
Pf20 = squeeze(Pf_result(14,:,:)); % Pf for mean age 14
fname = ['output' sprintf('%s',run_name) '.mat']
save(fname,'result','COPE_AGE','Pf_var','sigma_pf','Pf_result') % saving output

% Calculating released C, N and P of mesozooplankton
release_result = [];
for i = 1:num_ex;
    loss_nut = mesoz_release_func(X_ave(i,:),cope_age(i));
    release_result = cat(3,release_result,loss_nut);
end
release_result = squeeze(release_result);
release_result = permute(release_result,[2,1]);
release_C = release_result(:,1);
release_N = release_result(:,2);
release_P = release_result(:,3);

% --------------------------------------
%-----------------------------------------
% Plotting results 
Phyto_C = X_ave(:,4)+X_ave(:,8)+X_ave(:,12)+X_ave(:,16); % mgC/m^3
Phyto_N = X_ave(:,6)+X_ave(:,10)+X_ave(:,14)+X_ave(:,18); % mmol N/m^3
Phyto_P = X_ave(:,7)+X_ave(:,11)+X_ave(:,15)+X_ave(:,19); % mmol P/m^3
Zoo_C = X_ave(:,20)+X_ave(:,23)+X_ave(:,26); % mgC/m^3
Zoo_N = X_ave(:,21)+X_ave(:,24)+X_ave(:,42); % mmol N/m^3
Zoo_P = X_ave(:,22)+X_ave(:,25)+X_ave(:,43); % mmol P/m^3
POC = X_ave(:,30)+X_ave(:,33)+X_ave(:,36); % mgC/m^3
PON = X_ave(:,31)+X_ave(:,34)+X_ave(:,37); % mmol N/m^3
POP = X_ave(:,32)+X_ave(:,35)+X_ave(:,38); % mmol P/m^3

figure(1)
suptitle(fname)
subplot(2,2,1) % C:P against cope_age
    plot(cope_age,(X_ave(:,12)./12)./X_ave(:,15),'LineWidth',0.5)  % C:P small plankton (pico-p)
    hold on
    plot(cope_age,(X_ave(:,8)./12)./X_ave(:,11),'LineWidth',0.5)  % C:P small plankton (flagellates)
    hold on
    plot(cope_age,(X_ave(:,16)./12)./X_ave(:,19),'LineWidth',0.5)  % C:P dinoflagellates
    hold on
    plot(cope_age,(X_ave(:,4)./12)./X_ave(:,7),'LineWidth',0.5)  % C:P diatoms
    hold on
    plot(cope_age,Phyto_C./12./Phyto_P,'LineWidth',2.0); % C:P plankton
    hold on
    plot(cope_age,(X_ave(:,20)./12)./X_ave(:,22),'--','LineWidth',0.5)  % C:P heterotrophs
    hold on
    plot(cope_age,(X_ave(:,23)./12)./X_ave(:,25),'--','LineWidth',0.5)  % C:P microz
    hold on
    plot(cope_age,(X_ave(:,26)./12)./X_ave(:,43),'--','LineWidth',0.5)  % C:P mesoz
    hold on
    plot(cope_age,Zoo_C./12./Zoo_P,'--','LineWidth',2.0); % C:P Zoo (with meso z)
    hold on
    plot(cope_age,POC./12./POP,'k','LineWidth',2.0); % POM
    hold on
    plot(cope_age,((X_ave(:,27))./12)./(X_ave(:,29)),'-.','LineWidth',1.5); % Bacteria
    xlim([0 14])
    ylim([0 450])
    ylabel('C:P')
    hold on
    xlabel('Mean age (days)')
    hold on
    xticks([0:1:14])
    hold on
   legend('pico-p','flagellates','dino-flagellates','diatoms','Bulk P','pico z', 'micro z','mesoz','Bulk Z','POM','Bacteria','Location','northeast')
    title('Mean age against C:P')
subplot(2,2,2) % Plotting C:N against cope_age
    plot(cope_age,(X_ave(:,12)./12)./X_ave(:,14),'LineWidth',0.5)  % small plankton (pico-p)
    hold on
    plot(cope_age,(X_ave(:,8)./12)./X_ave(:,10),'LineWidth',0.5)  % small plankton (flagellates)
    hold on
    plot(cope_age,(X_ave(:,16)./12)./X_ave(:,18),'LineWidth',0.5)  % dinoflagellates
    hold on
    plot(cope_age,(X_ave(:,4)./12)./X_ave(:,6),'LineWidth',0.5)  % diatoms
    hold on
    plot(cope_age,Phyto_C./12./Phyto_N,'LineWidth',2.0); %phytoplankton
    hold on
    plot(cope_age,(X_ave(:,20)./12)./X_ave(:,21),'--','LineWidth',0.5)  % heterotrophs
    hold on
    plot(cope_age,(X_ave(:,23)./12)./X_ave(:,24),'--','LineWidth',0.5)  % microz
    hold on
    plot(cope_age,(X_ave(:,26)./12)./X_ave(:,42),'--','LineWidth',0.5)  % mesoz
    hold on
    plot(cope_age,Zoo_C./12./Zoo_N,'--','LineWidth',2.0); % Zoo (with meso z)
    hold on
    plot(cope_age,POC./12./PON,'k','LineWidth',2.0); % POM
    hold on
    plot(cope_age,((X_ave(:,27))./12)./(X_ave(:,28)),'-.','LineWidth',1.5); % Bacteria
    xlim([0 14])
    ylabel('C:N')
    hold on
    xlabel('Mean age (days)')
    hold on
    xticks([0:1:14])
    hold on
    xlim([0 14])
    ylim([0 30])
     legend('pico-p','flagellates','dino-flagellates','diatoms','Bulk P','pico z', 'micro z','mesoz','Bulk Z','POM','Bacteria','Location','northeast')
    title('Mean age against C:N')
subplot(2,2,3) % N:P
    plot(cope_age,(X_ave(:,14))./X_ave(:,15),'LineWidth',0.5)  % small plankton (pico-p)
    hold on
    plot(cope_age,(X_ave(:,10))./X_ave(:,11),'LineWidth',0.5)  % small plankton (flagellates)
    hold on
    plot(cope_age,(X_ave(:,18))./X_ave(:,19),'LineWidth',0.5)  % dinoflagellates
    hold on
    plot(cope_age,(X_ave(:,6))./X_ave(:,7),'LineWidth',0.5)  % diatoms
    hold on
    plot(cope_age,Phyto_N./Phyto_P,'LineWidth',2.0); % phytoplankton
    hold on
    plot(cope_age,(X_ave(:,21))./X_ave(:,22),'--','LineWidth',0.5)  % heterotrophs
    hold on
    plot(cope_age,(X_ave(:,24))./X_ave(:,25),'--','LineWidth',0.5)  % microz
    hold on
    plot(cope_age,(X_ave(:,42))./X_ave(:,43),'--','LineWidth',0.5)  % mesoz
    hold on
    plot(cope_age,Zoo_N./Zoo_P,'--','LineWidth',2.0); % Zoo (with meso z)
    hold on
    plot(cope_age,PON./POP,'k','LineWidth',2.0); % POM
    hold on
    plot(cope_age,((X_ave(:,28)))./(X_ave(:,29)),'-.','LineWidth',1.5); % Bacteria
    hold on
    plot(cope_age,((X_ave(:,1)+X_ave(:,2)))./(X_ave(:,3)),':','LineWidth',1.5); % DIN:DIP
    hold on
    xlim([0 14])
    ylim([0 64])
    ylabel('N:P')
    hold on
    xlabel('Mean age (days)')
    hold on
    xticks([0:1:14])
    hold on
    legend('pico-p','flagellates','dino-flagellates','diatoms','Bulk P','pico z', 'micro z','mesoz','Bulk Z','POM','Bacteria','DIN:DIP','Location','northeast')
    title('Mean age against N:P')
subplot(2,2,4) % Plotting Biomass against cope_age
    plot(cope_age,X_ave(:,12),'LineWidth',0.5)  % small plankton (pico-p)
    hold on
    plot(cope_age,X_ave(:,8),'LineWidth',0.5)  % small plankton (flagellates)
    hold on
    plot(cope_age,X_ave(:,16),'LineWidth',0.5)  % dinoflagellates
    hold on
    plot(cope_age,X_ave(:,4),'LineWidth',0.5)  % diatoms
    hold on
    plot(cope_age,Phyto_C,'LineWidth',2.0); %phytoplankton
    hold on
    plot(cope_age,X_ave(:,20),'--','LineWidth',0.5)  % heterotrophs
    hold on
    plot(cope_age,X_ave(:,23),'--','LineWidth',0.5)  % microz
    hold on
    plot(cope_age,X_ave(:,26),'--','LineWidth',0.5)  % mesoz
    hold on
    plot(cope_age,Zoo_C,'--','LineWidth',2.0); % Zoo (with meso z)
    hold on
    plot(cope_age,X_ave(:,27),'-.','LineWidth',1.5); % Bacteria
    xlim([0.5 10])
    ylabel('Biomass (mgC/m^{3})')
    hold on
    xlabel('Mean age (days)')
    hold on
    xticks([0:1:14])
    hold on
    xlim([0 14])
     legend('pico-p','flagellates','dino-flagellates','diatoms','Bulk P','pico z', 'micro z','mesoz','Bulk Z','Bacteria','Location','northeast')
    title('Mean age against biomass')
figure(3)
clf
suptitle('POM and DOM')
subplot(2,2,1) % POC:POP and DOC:DOP
    plot(cope_age,(X_ave(:,30)./12)./X_ave(:,32),'LineWidth',0.5)  % SPOC:SPOP
    hold on
    plot(cope_age,(X_ave(:,33)./12)./X_ave(:,35),'LineWidth',0.5)  % MPOC:MPOP
    hold on
    plot(cope_age,(X_ave(:,36)./12)./X_ave(:,38),'LineWidth',0.5)  % LPOC:LPOP
    hold on    
    plot(cope_age,POC./12./POP,'LineWidth',2.0)  % POC:POP
    hold on    
    plot(cope_age,X_ave(:,39)./12./X_ave(:,41),'LineWidth',2.0)  % DOC:DOP
    ylabel('C:P')
    hold on
    xlabel('Mean age (days)')
    hold on
    xticks([0:1:14])
    hold on
    xlim([0 14])
    ylim([0 700])
    title('POC:POP and DOC:DOP')
    legend('SPOM','MPOM','LPOM','Bulk POM','DOM')
subplot(2,2,2) % POC:PON and DOC:DON
    plot(cope_age,(X_ave(:,30)./12)./X_ave(:,31),'LineWidth',0.5)  % SPOC:SPON
    hold on
    plot(cope_age,(X_ave(:,33)./12)./X_ave(:,34),'LineWidth',0.5)  % MPOC:MPON
    hold on
    plot(cope_age,(X_ave(:,36)./12)./X_ave(:,37),'LineWidth',0.5)  % LPOC:LPON
    hold on    
    plot(cope_age,POC./12./PON,'LineWidth',2.0)  % POC:LPON
    hold on    
    plot(cope_age,X_ave(:,39)./12./X_ave(:,40),'LineWidth',2.0)  % DOC:DON
    ylabel('C:N')
    hold on
    xlabel('Mean age (days)')
    hold on
    xticks([0:1:14])
    hold on
    xlim([0 14])
    ylim([0 60])
    title('POC:PON and DOC:DON')
    legend('SPOM','MPOM','LPOM','Bulk POM','DOM')
subplot(2,2,3) % PON:POP
    plot(cope_age,(X_ave(:,31))./X_ave(:,32),'LineWidth',0.5)  % SPON:SPOP
    hold on
    plot(cope_age,(X_ave(:,34))./X_ave(:,35),'LineWidth',0.5)  % MPON:MPOP
    hold on
    plot(cope_age,(X_ave(:,37))./X_ave(:,38),'LineWidth',0.5)  % LPON:LPOP
    hold on    
    plot(cope_age,PON./POP,'LineWidth',2.0)  % PON:POP
    hold on    
    plot(cope_age,X_ave(:,40)./X_ave(:,41),'LineWidth',2.0)  % DON:DOP
    ylabel('N:P')
    hold on
    xlabel('Mean age (days)')
    hold on
    xticks([0:1:14])
    yticks([0:10:60])
    hold on
    xlim([0 14])
    ylim([0 64])
    title('PON:POP and DON:DOP')
    legend('SPOM','MPOM','LPOM','Bulk POM','DOM')
subplot(2,2,4) % POC and DOC
    plot(cope_age,X_ave(:,30),'LineWidth',0.5)  % SPOC
    hold on
    plot(cope_age,X_ave(:,33),'LineWidth',0.5)  % MPOC
    hold on
    plot(cope_age,X_ave(:,36),'LineWidth',0.5)  % LPOC
    hold on    
    plot(cope_age,POC,'LineWidth',2.0)  % Bulk POC
    hold on    
    plot(cope_age,X_ave(:,39),'LineWidth',2.0)  % DOC
    ylabel('Carbon(mgC m^{-3})')
    hold on
    xlabel('Mean age (days)')
    hold on
    xticks([0:1:14])
    hold on
    xlim([0 14])
    %ylim([0 700])
    title('POC and DOC')
    legend('SPOM','MPOM','LPOM','Total POM','DOC')
figure(4) % nutrient conc.
clf
suptitle(fname)
    subplot(2,2,1)
        plot(cope_age,X_ave(:,1),'LineWidth',0.5)  % NO3
        ylabel('NO3 (mmol m^{-3})')
        xlabel('Mean age (days)')
        hold on
        xticks([0:1:14])
        hold on
        title('Nitrate')
        hold on
        xlim([0 14])
    subplot(2,2,2)
        plot(cope_age,X_ave(:,2),'LineWidth',0.5)  % NH4
        ylabel('NH4 (mmol m^{-3})')
        xlabel('Mean age (days)')
        hold on
        xticks([0:1:14])
        hold on
        title('Ammonium')
        hold on
        xlim([0 14])
    subplot(2,2,3)
        plot(cope_age,X_ave(:,3),'LineWidth',0.5)  % PO4
        ylabel('PO4 (mmol m^{-3})')
        xlabel('Mean age (days)')
        hold on
        xticks([0:1:14])
        hold on
        title('Phosphate')
        hold on
        xlim([0 14])
    subplot(2,2,4)
        plot(cope_age,(X_ave(:,1)+X_ave(:,2))./X_ave(:,3),'LineWidth',0.5)  % NH4/PO4
        hold on
        plot(cope_age,(X_ave(:,37))./X_ave(:,38),'LineWidth',0.5)  % lPON/lPOP
        hold on
        plot(cope_age,(X_ave(:,40))./X_ave(:,41),'LineWidth',0.5)  % lDON/lDOP
        hold on
        plot(cope_age,(X_ave(:,45))./X_ave(:,46),'LineWidth',0.5)  % slDON/slDOP
        hold on
        plot(cope_age,(X_ave(:,42))./X_ave(:,43),'LineWidth',0.5)  % mesoz
        ylabel('N/P')
        xlabel('Mean age (days)')
        hold on
        plot(cope_age,release_N./release_P) % mesoz release
        xticks([0:1:14])
        hold on
        legend('DIN/DIP','lPON/lPOP','lDON/lDOP','slDON/slDOP','mesoz','release N/P')
        hold on
        xlim([0 14])
figure(5) % Body N:P for mesozvs release N:P and lPON:lPOP 
clf
suptitle(fname)
    body_np = X_ave(:,42)./X_ave(:,43);
    plot(body_np,release_N./release_P)
    hold on
    plot(body_np,X_ave(:,37)./X_ave(:,38))
    hold on
    plot(body_np,body_np,'--')
    ylabel('N:P')
    xlabel('body N:P')
    xlim([5 35])
    ylim([5 35])
    legend('release N/P','lPON/lPOP','1:1 line')
    axis square