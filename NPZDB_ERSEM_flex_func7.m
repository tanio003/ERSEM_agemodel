% NPZDB_ERSEM_flex_func7.m
% Function file for ERSEM NPZDB model based on Blackford et al. 2004.
% version 7 Tatsuro Tanioka 08/29/17
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
function dXdt=NPZDB_ERSEM_flex_func7(t,X)
% Variables
% 1 = Nn Nitrate
% 2 = Na Ammonium 
% 3 = Np Phosphates
% 4 = P1c = Diatoms
% 5 = P1chl 
% 6 = P1n
% 7 = P1p
% 8 = P2c = Flagellates
% 9 = P2chl
% 10 = P2n
% 11 = P2p
% 12 = P3c = pico-plankton
% 13 = P3chl
% 14 = P3n
% 15 = P3p
% 16 = P4c = dino-flagellates
% 17 = P4chl
% 18 = P4n
% 19 = P4p
% 20 = z1c = Heterotrophic flagellates
% 21 = z1n
% 22 = Z1p
% 23 = z2c = Micro-zooplankton
% 24 = z2n
% 25 = z2p
% 26 = z3c = Meso-zooplankton
% 27 = Bc = Bacteria
% 28 = Bn
% 29 = Bp
% 30 = sPOC = small POM
% 31 = sPON
% 32 = sPOP
% 33 = mPOC = medium POM
% 34 = mPON
% 35 = mPOP
% 36 = lPOC = large POM 
% 37 = lPON
% 38 = lPOP
% 39 = DOC = labile DOM
% 40 = DON
% 41 = DOP
% 42 = z3n = Meso-zooplankton N
% 43 = z3p = Meso-zooplankton P
% 44 = slDOC
% 45 = slDON 
% 46 = slDOP
% Mass Balance Check Variables
% 47 = Total N input
% 48 = Total N loss
% 49 = Total P input
% 50 = Total P loss

% declare some global variables 
% Ecosystem setting
global n_P % [unitless] Number of phytoplankton types
global n_Z % [unitless] Number of plankton types
global n_var
global A  %[/day]     CHEMOSTAT DILUTION RATE 
global Nin_N  %[mmol m^-3]     INPUT NO3 CONCENTRATION (N limited, N/P co-limitation, P limited)
global Nin_A  % [mmol m^-3]    INPUT NH4 CONCENTRATION 
global Nin_P  %[mmol m^-3]     INPUT P CONCENTRATION
global Q10    %[unitless]        Q10 value (Butenschön et al., 2016)
global r_nitr  %[/day]       Max. Nitrification rate (Butenschön et al., 2016)
global h_N_nitr %[mmol m^-3] Cubis MMenten constant for nitrification (Butenschön et al., 2016)
global temperature %[C] Temperature
global PAR %[W m-2] Insolation
global O_rel % [mmol m^-3] % relative oxygen saturation
global epsilon %[mgC m^-3] Negativity control factor = 0.01 (Butenschön et al., 2016)
global epsilon_p %[mgC m^-3] Negativity control factor for P = 0.01 * (P:C)redfield (Butenschön et al., 2016)
global epsilon_n %[mgC m^-3] Negativity control factor for N = 0.01 * (N:C)redfield (Butenschön et al., 2016)
global epsilon_chl %[mgC m^-3] Negativity control factor for Chl = 0.01 * (chl:C)min (Butenschön et al., 2016)
global Q_N_RED %[mmol N /mgC] Redfield N:C Ratio
global Q_P_RED %[mmol P /mgC] Redfield P:C Ratio

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
global ST %[unitless] level of nutrient stress below which phytoplankton subsides
global SIZE_PHYTO % [um] typical diameter of phytoplankton [Diatoms Flagellates Pico-plankton Dino-flagellates]

% Grazer Parameters
global R_ASS_ZOO      %[/day]  Assimilation rate at 10C
global H %[mgC m^-3] Food conc. where relative uptake is 0.5
global AE % [1] Assimilation efficiency
global P_DOM % Fraction of excretion going to DOM
global R_RESTR_ZOO % [/day] Basal respiration at 10C 
global H_OXMORT % [mmol m^-3] Ox. saturation where respiration is 0.5;
global EU % [1] Excreted fraction of uptake
global R_MORTOX % [1/day] Mortality rate due to ox. limitation
global R_MORT % [1/day] Temp. independent mortality rate
global QMAX_N_ZOO  %[mmolN / mgC]  MAX N Quota for zooplankton
global QMAX_P_ZOO  %[mmolP / mgC]  MAX P Quota for zooplankton
global Z_MINFOOD  %[mgC m^-3]    Lower Foodthreshold for feeding
global cn %[1/day] relaxation rate of N-excretion 
global cp %[1/day] relaxation rate of P-excretion 
global q_Rexcr_meso % Excreted fraction of POM uptake by mesozooplankton

% Newly added global parameters
global Pf %[unitless] prey availability matrix
global sigma_pf % [1.] standard deviation of the grazing spectrum matrix
global can_threshold % [1.] threshold value for cannibalism 1 = no threshold
global Pf_var % [1.] Switch to use flexible Pf (=1) for fixed Pf (=0)
global rho_TER % [1.] max.GGE(N)/max.GGE(P)
global COPE_AGE % Mean age of copepod (days)
global COPE_AGE_switch % [1.] Switch to use copepod age model (1 = true, 0 =false)
global COPE_AGE_cpslope % Slope of C:P against age (d-1)
global COPE_AGE_cpintercept % y-intercept of C:P against age (molC/molP)
global COPE_AGE_npslope % Slope of N:P against age (d-1)
global COPE_AGE_npintercept % y-intercept of N:P against age (molN/molP)
global COPE_AGE_rmaxslope % Slope of rmax against age (d-1 d-1)
global COPE_AGE_rmaxintercept % y-intercept of rmax against age (d-1)

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

% Initialize ouput column vector
dXdt = zeros(n_var,1);

% temperature dependence 
t_water = temperature;
n = round(t); %round t to the next integer(t~=n);
temp = Q10.^((t_water(n)-10)./10) - Q10.^((t_water(n)-32)./3); % temperature function eq. (239), B16 

% Phytoplankton functions
P_C_phyto = [X(4) X(8) X(12) X(16)];
P_Chl_phyto = [X(5) X(9) X(13) X(17)];
P_N_phyto = [X(6) X(10) X(14) X(18)];
P_P_phyto = [X(7) X(11) X(15) X(19)];

P_C_phyto_prime = P_C_phyto-epsilon.*ones(1,n_P);
P_N_phyto_prime = P_N_phyto-epsilon_n.*ones(1,n_P);
P_P_phyto_prime = P_P_phyto-epsilon_p.*ones(1,n_P);
P_Chl_phyto_prime = P_Chl_phyto - epsilon_chl.*ones(1,n_P);

Q_N_phyto = P_N_phyto./P_C_phyto;
Q_P_phyto = P_P_phyto./P_C_phyto;
Q_N_phyto_prime = P_N_phyto_prime./P_C_phyto;
Q_P_phyto_prime = P_P_phyto_prime./P_C_phyto;
theta_phyto = P_Chl_phyto./P_C_phyto; 
theta_phyto_prime = P_Chl_phyto_prime./P_C_phyto; 

fi = (1.0 - exp((-1.*ALPHA.*PAR.*theta_phyto)./(R_ASS.*temp))).*exp((-1.*BETA*PAR.*theta_phyto)./(R_ASS.*temp)); % Light limitation eq. (5)
N_stress = (Q_N_phyto - QMIN_N_PHYTO)./(QMAX_N_PHYTO - QMIN_N_PHYTO); % N stress for eq. (8)
P_stress = (Q_P_phyto - QMIN_P_PHYTO)./(QMAX_P_PHYTO - QMIN_P_PHYTO); % P stress for eq. (8)
n_stress = nanmax(0,nanmin(N_stress,P_stress)); % Nutreint stress  eq. (242) of B16
n_stress = nanmin(1,(double(n_stress>0).*n_stress)); % Use this value only when a real number and constrain between 0-1 
photosynthesis = R_ASS.* temp.* fi.* P_C_phyto; % Photosynthesis eq. (4)
excretion_C_phyto = photosynthesis .* (P_EX + (ones(1,n_P)-n_stress).*(ones(1,n_P)-P_EX)); % excretion of excess C eq. (9)
respiration_phyto = R_RESTR.*temp.*(P_C_phyto-epsilon*ones(1,n_P)) + (photosynthesis-excretion_C_phyto).*R_ACTR; % respiration eq. (10)
lysis = (1.0./(n_stress+0.1)) .* R_LYSIS.* (P_C_phyto-epsilon); % lysis eq. (7)
lysis_N = lysis.*Q_N_phyto_prime; % loss of N by lysis
lysis_P = lysis.*Q_P_phyto_prime; % loss of P by lysis
phi_phyto = n_stress.*((THETA_MAX-THETA_MIN).*((R_ASS.*temp.*fi)./(ALPHA.*PAR.*theta_phyto_prime))); % proportion of photosynthate directed to Chl eq. (10) B16 
uptake_N_required = (photosynthesis-excretion_C_phyto-respiration_phyto).*QMAX_N_PHYTO + QMAX_N_PHYTO.*P_C_phyto-(P_N_phyto-epsilon_n.*ones(1,n_P)); % eq. (13)
uptake_N_available = (A_N.*(X(1)-epsilon_n) + A_A.*(X(2)-epsilon_n)).*P_C_phyto;

uptake_P_required = (photosynthesis-excretion_C_phyto-respiration_phyto).*QMAX_P_PHYTO + QMAX_P_PHYTO.*P_C_phyto-(P_P_phyto-epsilon_p.*ones(1,n_P)); % eq. (11)
uptake_P_available = A_P.*(X(3)-epsilon_p).*P_C_phyto; 

excretion_N_phyto = -1.*uptake_N_required.*double((uptake_N_required<=0)); % eq. (16) B16
uptake_N = double((uptake_N_required>0)).*nanmin(uptake_N_required,uptake_N_available); % eq. (16) B16
excretion_P_phyto = -1.*uptake_P_required.*double((uptake_P_required<=0)); % eq. (16) B16
uptake_P = double((uptake_P_required>0)).*nanmin(uptake_P_required,uptake_P_available); % eq. (16) B16

sinking_phyto = RSED_PHYTO.*nanmax(0.0,(ST-n_stress)).*(P_C_phyto-epsilon.*ones(1,n_P)); % sinking under stress eq. (16)

% Zooplankton Functions
    % Define Prey array
      Prey = [X(27) X(12) X(8) X(4) X(16) X(20) X(23) X(26) X(33)]; %[bacteria picop flage diatoms dinoflag sz mz lz mPOC]
      Prey_N = [X(28) X(14) X(10) X(6) X(18) X(21) X(24) X(42) X(34)];%[bacteria picop flage diatoms dinoflag sz mz lz mPON]
      Prey_P = [X(29) X(15) (X(11)) (X(7)) (X(19)) (X(22)) (X(25)) X(43) X(35)]; %[bacteria picop flage diatoms dinoflag sz mz lz lPON]
      Q_N_Prey = Prey_N./Prey;  
      Q_P_Prey = Prey_P./Prey;
    % Define Predator array
      Predator = [X(20) X(23) X(26)]; % Predator C [small_z micro_z meso_z]
      Predator_N = [X(21) X(24) X(42)]; % Predator N [small_z micro_z meso_z]
      Predator_P = [X(22) X(25) X(43)]; % Predator N [small_z micro_z meso_z]
      Q_N_Predator = Predator_N./Predator;
      Q_P_Predator = Predator_P./Predator;
    % Construct Relative prey availability matrix Pf_flex
      N_mat = Q_N_Predator./Q_N_Prey';
      P_mat = Q_P_Predator./Q_P_Prey';
      food_z1 = [Prey(1) Prey(2) Prey(3) 0 0 Prey(6) 0 0 0]'; % Food for Z1
      food_z2 = [Prey(1) Prey(2) Prey(3) Prey(4) Prey(5) Prey(6) Prey(7) 0 0]'; % Food for Z2
      food_z3 = [0 0 Prey(3) Prey(4) Prey(5) Prey(6) Prey(7) Prey(8) Prey(9)]'; % Food for Z3
      food = [food_z1 food_z2 food_z3];
    % Multiply prey biomass by density independent prey preference 
      food_frac(:,1) = food(:,1).*Pf(:,1);
      food_frac(:,2) = food(:,2).*Pf(:,2);
      food_frac(:,3) = food(:,3).*Pf(:,3);
    % Food quality
      Q_NP_optfood = Predator_N./Predator_P.*rho_TER; % optimal food N:P
      %x = ((Q_N_Prey'./Q_P_Prey') - ones(9,3).*(Q_N_RED./Q_P_RED))./(Q_N_RED./Q_P_RED); % (N:P(prey) - N:P(optimal))/N:P(optimal)
      x = ((Q_N_Prey'./Q_P_Prey') - ones(9,3).*Q_NP_optfood)./(ones(9,3).*Q_NP_optfood); % (N:P(prey) - N:P(optimal))/N:P(optimal)
      fq = normpdf(x,0.0,sigma_pf); % food quality
      
      % Assign 0 to constrained relationships
      fq(1,3) = 0; % meso z can't eat bacteria
      fq(2,3) = 0; % meso z can't eat picophytoplankton
      fq(4,1) = 0; % small z can't eat diatoms
      fq(5,1) = 0; % small z can't eat dinoflagellates
      fq(7,1) = 0; % small z can't eat micro z
      fq(8,1) = 0; % small z can't eat meso z
      fq(8,2) = 0; % micro z can't eat meso z
      fq(9,1) = 0; % small z can't mPOM
      fq(9,2) = 0; % micro z can't mPOM
      
      if Pf_var == 0;
        Pf_flex = Pf;
      elseif Pf_var == 1;
        Pf_flex = fq.*food_frac; % Food preference = food quality + Fasham
      elseif Pf_var == 2;
        Pf_flex = food_frac; % Food preference = Fasham formulation
      else
        Pf_flex = Pf; % Fixed Preference
      end

      % Normalize each sum of each preadtor's total preference to 1
      Pf_flex(:,1) = Pf_flex(:,1)./(sum(Pf_flex(:,1)));
      Pf_flex(:,2) = Pf_flex(:,2)./(sum(Pf_flex(:,2)));
      Pf_flex(:,3) = Pf_flex(:,3)./(sum(Pf_flex(:,3)));
      % Cannabalsim control for small z and micro z when exceeds threshold
      Pf_flex(6,1) = Pf_flex(6,1).*double(Pf_flex(6,1)<=can_threshold);
      Pf_flex(7,2) = Pf_flex(7,2).*double(Pf_flex(7,2)<=can_threshold);
      Pf_flex(8,3) = Pf_flex(8,3).*double(Pf_flex(8,3)<=can_threshold);
      % Re-Normalize each sum of each preadtor's total preference to 1 - can_threshold 
      Pf_flex(:,1) = (Pf_flex(:,1).*(1-can_threshold.*double(Pf_flex(6,1)<=0)))/(sum(Pf_flex(:,1)));
      Pf_flex(:,2) = (Pf_flex(:,2).*(1-can_threshold.*double(Pf_flex(7,2)<=0)))/(sum(Pf_flex(:,2)));
      Pf_flex(:,3) = (Pf_flex(:,3).*(1-can_threshold.*double(Pf_flex(8,3)<=0)))/(sum(Pf_flex(:,3)));
      % Insert the Canibalism threshold term
      Pf_flex(6,1) = Pf_flex(6,1) + can_threshold.*double(Pf_flex(6,1)<=0); % small z like can_const
      Pf_flex(7,2) = Pf_flex(7,2) + can_threshold.*double(Pf_flex(7,2)<=0); % micro z like can_const
      Pf_flex(8,3) = Pf_flex(8,3) + can_threshold.*double(Pf_flex(8,3)<=0); % micro z like can_const
      % Final Pf matrix
      Pf_fix = Pf;
      Pf_flex = Pf_flex.*double(Pf_var>0) + Pf_fix.*double(Pf_var<=0);  
      Pf_flex(:,1) =  Pf_fix(:,1); % small z has fixed preference
      Pf_flex(:,2) =  Pf_fix(:,2); % micro z has fixed preference
      %Pf_flex(:,3) =  Pf_fix(:,3);

      Pf_flex;
      
      % 1. Heterotrophic flagellates functions (X(20)-X(22))
        Q_N_sz = X(21)./X(20); % sz N quota
        Q_N_sz_prime = (X(21)-epsilon_n)./X(20); 
        Q_P_sz = X(22)./X(20); % sz P quota
        Q_P_sz_prime = (X(22)-epsilon_p)./X(20); % sz P quota
      % Total prey available in C,N,P eq. (29)
        FCtot_sz = nansum(Pf_flex(:,1)'.*(Prey.*Prey./(Prey+Z_MINFOOD(1)))); % Total food available for sz 
        FCtot_N_sz = nansum(Pf_flex(:,1)'.*(Prey_N.*Prey_N./(Prey_N+Z_MINFOOD(1).*Q_N_Prey))); % Total food N available for sz 
        FCtot_P_sz = nansum(Pf_flex(:,1)'.*(Prey_P.*Prey_P./(Prey_P+Z_MINFOOD(1).*Q_P_Prey))); % Total food P available for sz 
      % Uptake of C, N, P  eq. (30-32)
        uptake_sz = R_ASS_ZOO(1).*FCtot_sz./(FCtot_sz + H(1)).*temp.*X(20); % uptake of carbon
        uptake_N_sz = R_ASS_ZOO(1).*FCtot_N_sz./(FCtot_N_sz + H(1).*(FCtot_N_sz./FCtot_sz)).*temp.*X(21); % uptake of N
        uptake_P_sz = R_ASS_ZOO(1).*FCtot_P_sz./(FCtot_P_sz + H(1).*(FCtot_P_sz./FCtot_sz)).*temp.*X(22); % uptake of P 
      % Flux from each individidual food [bacteria picop flage diatoms dinoflag sz mz lz mPOM] % eq. (31)
        flux_sz = Pf_flex(:,1)'.*(Prey.*Prey./(Prey+Z_MINFOOD(1))).*uptake_sz./FCtot_sz; % flux C
        flux_N_sz = Pf_flex(:,1)'.*(Prey_N.*Prey_N./(Prey_N+Z_MINFOOD(1).*Q_N_Prey)).*uptake_N_sz./FCtot_N_sz;  % flux N
        flux_P_sz = Pf_flex(:,1)'.*(Prey_P.*Prey_P./(Prey_P+Z_MINFOOD(1).*Q_P_Prey)).*uptake_P_sz./FCtot_P_sz;  % flux P
      % Respiration of C, eq. (35)
        respiration_sz = R_RESTR_ZOO(1).*temp.*(X(20)-epsilon) + uptake_sz.*(1.0-AE(1)).*(1.0-EU(1));
      % Excretion of C, N, P, eq. (34)
        excretion_sz = uptake_sz.*(1.0-AE(1)).*EU(1); % excretion C
        excretion_N_sz = uptake_N_sz.*(1.0-AE(1)).*EU(1); % excretion N
        excretion_P_sz = uptake_P_sz.*(1.0-AE(1)).*EU(1); % excretion P
      % Loss of C, N, P due to Mortality, % eq. (37)
        eO2_sz = (1.0+H_OXMORT(1)).*O_rel./(O_rel + H_OXMORT(1)); % oxygen limitation eq. (249) 
        mortality_sz = ((1.0-eO2_sz).*R_MORTOX(1) + R_MORT(1)).*(X(20)-epsilon); % C loss
        mortality_N_sz = mortality_sz.*Q_N_sz_prime; % N loss
        mortality_P_sz = mortality_sz.*Q_P_sz_prime; % P loss
      % Release of excess N and P, eq. (36)
        loss_N_sz = nanmax(0,((Q_N_sz_prime-QMAX_N_ZOO(1)).*(X(20)-epsilon).*cn(1))); 
        loss_P_sz = nanmax(0,((Q_P_sz_prime-QMAX_P_ZOO(1)).*(X(20)-epsilon).*cp(1)));  

    % 2. Micro-Zooplankton functions (X(23)-X(25))
        Q_N_mz = X(24)./X(23); % mz N quota
        Q_N_mz_prime = (X(24)-epsilon_n)./X(23);
        Q_P_mz = X(25)./X(23); % mz P quota
        Q_P_mz_prime = (X(25)-epsilon_p)./X(23); % mz P quota
      % Total prey available in C,N,P eq. (29)
        FCtot_mz = nansum(Pf_flex(:,2)'.*(Prey.*Prey./(Prey+Z_MINFOOD(2)))); % Total food available 
        FCtot_N_mz = nansum(Pf_flex(:,2)'.*(Prey_N.*Prey_N./(Prey_N+Z_MINFOOD(2).*Q_N_Prey))); % Total food N available 
        FCtot_P_mz = nansum(Pf_flex(:,2)'.*(Prey_P.*Prey_P./(Prey_P+Z_MINFOOD(2).*Q_P_Prey))); % Total food P available 
     %  Uptake of C, N, P  eq. (30-32)
        uptake_mz = R_ASS_ZOO(2).*FCtot_mz./(FCtot_mz + H(2)).*temp.*X(23); % uptake of carbon 
        uptake_N_mz = R_ASS_ZOO(2).*FCtot_N_mz./(FCtot_N_mz + H(2).*(FCtot_N_mz./FCtot_mz)).*temp.*X(24); % uptake of N
        uptake_P_mz = R_ASS_ZOO(2).*FCtot_P_mz./(FCtot_P_mz + H(2).*(FCtot_P_mz./FCtot_mz)).*temp.*X(25); % uptake of P
     % Flux from each individidual food [bacteria picop flage diatoms
      % dinoflag sz mz lz mPOM] % eq. (31)
        flux_mz = Pf_flex(:,2)'.*(Prey.*Prey./(Prey+Z_MINFOOD(2))).*uptake_mz./FCtot_mz; 
        flux_N_mz = Pf_flex(:,2)'.*(Prey_N.*Prey_N./(Prey_N+Z_MINFOOD(2).*Q_N_Prey)).*uptake_N_mz./FCtot_N_mz;
        flux_P_mz = Pf_flex(:,2)'.*(Prey_P.*Prey_P./(Prey_P+Z_MINFOOD(2).*Q_P_Prey)).*uptake_P_mz./FCtot_P_mz;
      % Respiration of C, eq. (35)
        respiration_mz = R_RESTR_ZOO(2).*temp.*(X(23)-epsilon) + uptake_mz.*(1.0-AE(2)).*(1.0-EU(2)); % eq. (21) 
      % Excretion of C, N, P, eq. (34)
        excretion_mz = uptake_mz.*(1.0-AE(2)).*EU(2); % 
        excretion_N_mz = uptake_N_mz.*(1.0-AE(2)).*EU(2); % 
        excretion_P_mz = uptake_P_mz.*(1.0-AE(2)).*EU(2); % 
      % Loss of C, N, P due to Mortality ,eq. (37)
        eO2_mz = (1.0+H_OXMORT(2)).*O_rel./(O_rel + H_OXMORT(2)); % 
        mortality_mz = ((1.0-eO2_mz).*R_MORTOX(2) + R_MORT(2)).*(X(23)-epsilon); % 
        mortality_N_mz = mortality_mz.*Q_N_mz_prime;
        mortality_P_mz = mortality_mz.*Q_P_mz_prime;
      
      % Release of excess N and P, eq. (36)
        loss_N_mz = nanmax(0,((Q_N_mz_prime-QMAX_N_ZOO(2)).*(X(23)-epsilon).*cn(2))); % 
        loss_P_mz = nanmax(0,((Q_P_mz_prime-QMAX_P_ZOO(2)).*(X(23)-epsilon).*cp(2))); % 
       
     % 3. Meso-Zooplankton functions (X(26), X(42-43))
     % ---- Age dependent functions ------------
     % Calculating C:N:P as a function of age
        Q_P_meso = 1.0/(COPE_AGE_cpslope.*COPE_AGE + COPE_AGE_cpintercept)./12.*double(COPE_AGE_switch > 0) + Q_P_RED.*double(COPE_AGE_switch < 1) ; % P:C [mol P/gC]according to Meunier16
        Q_NP_meso = (COPE_AGE_npslope.*COPE_AGE + COPE_AGE_npintercept).*double(COPE_AGE_switch > 0) + (Q_N_RED./Q_P_RED).*double(COPE_AGE_switch < 1); % N:P [mol N/ mol P] according to Meunier16
        Q_N_meso = Q_NP_meso.*Q_P_meso; % Calcualte N:C [mol N/gC]
     % Calculating max. assimilation rate as a function of age
        R_ASS_meso = (COPE_AGE_rmaxslope*COPE_AGE + COPE_AGE_rmaxintercept).*double(COPE_AGE_switch > 0) + R_ASS_ZOO(3).*double(COPE_AGE_switch < 1); % max.assimilation rate according to Meunier16
     % ----------------------------------------
     % Total prey available in C,N,P eq. (29)
        FCtot_lz = nansum(Pf_flex(:,3)'.*(Prey.*Prey./(Prey+Z_MINFOOD(3)))); % Total food available 
        FCtot_N_lz = nansum(Pf_flex(:,3)'.*(Prey_N.*Prey_N./(Prey_N+Z_MINFOOD(3).*Q_N_Prey))); % Total food N available 
        FCtot_P_lz = nansum(Pf_flex(:,3)'.*(Prey_P.*Prey_P./(Prey_P+Z_MINFOOD(3).*Q_P_Prey))); % Total food P available 
     %  Uptake of C, N, P  eq. (30-32)   
        %uptake_lz = R_ASS_ZOO(3).*FCtot_lz./(FCtot_lz + H(3)).*temp.*X(26); % uptake of carbon 
        uptake_lz = R_ASS_meso.*FCtot_lz./(FCtot_lz + H(3)).*temp.*X(26); % uptake of carbon 
        %uptake_N_lz = R_ASS_ZOO(3).*FCtot_N_lz./(FCtot_N_lz + H(3).*(FCtot_N_lz./FCtot_lz)).*temp.*X(26).*Q_N_RED; % uptake of N 
        %uptake_P_lz = R_ASS_ZOO(3).*FCtot_P_lz./(FCtot_P_lz + H(3).*(FCtot_P_lz./FCtot_lz)).*temp.*X(26).*Q_P_RED; % uptake of P 
        uptake_N_lz = R_ASS_meso.*FCtot_N_lz./(FCtot_N_lz + H(3).*(FCtot_N_lz./FCtot_lz)).*temp.*X(26).*Q_N_meso; % uptake of N 
        uptake_P_lz = R_ASS_meso.*FCtot_P_lz./(FCtot_P_lz + H(3).*(FCtot_P_lz./FCtot_lz)).*temp.*X(26).*Q_P_meso; % uptake of P 
     % Flux from each individidual food [bacteria picop flage diatoms dinoflag sz mz lz mPOM], eq. (31)
        flux_lz = Pf_flex(:,3)'.*(Prey.*Prey./(Prey+Z_MINFOOD(3))).*uptake_lz./FCtot_lz; 
        flux_N_lz = Pf_flex(:,3)'.*(Prey_N.*Prey_N./(Prey_N+Z_MINFOOD(3).*Q_N_Prey)).*uptake_N_lz./FCtot_N_lz;
        flux_P_lz = Pf_flex(:,3)'.*(Prey_P.*Prey_P./(Prey_P+Z_MINFOOD(3).*Q_P_Prey)).*uptake_P_lz./FCtot_P_lz;
     % Respiration of C, eq. (35)
        respiration_lz = R_RESTR_ZOO(3).*temp.*(X(26)-epsilon) + uptake_lz.*(1.0-AE(3)).*(1.0-EU(3)); % eq. (21) 
      % Excretion of C,N,P with enhanced excretion for scavenging on POM, eq. (41)
        excretion_lz = (uptake_lz-flux_lz(9)).*(1.0-AE(3)).*EU(3) + q_Rexcr_meso.*flux_lz(9); % 
        excretion_N_lz = (uptake_N_lz-flux_N_lz(9)).*(1.0-AE(3)).*EU(3) + q_Rexcr_meso.*flux_N_lz(9); % 
        excretion_P_lz = (uptake_P_lz-flux_P_lz(9)).*(1.0-AE(3)).*EU(3) + q_Rexcr_meso.*flux_P_lz(9); % 
      % Loss of C, N, P due to Mortality,eq. (37)
        eO2_lz = (1.0+H_OXMORT(3)).*O_rel./(O_rel + H_OXMORT(3)); % eq. (24) 
        mortality_lz = ((1.0-eO2_lz).*R_MORTOX(3) + R_MORT(3)).*(X(26)-epsilon); % eq. (23)
        %mortality_N_lz = mortality_lz.*Q_N_RED;
        %mortality_P_lz = mortality_lz.*Q_P_RED;
        mortality_N_lz = mortality_lz.*Q_N_meso;
        mortality_P_lz = mortality_lz.*Q_P_meso;
        
% Loss due to predation for each prey
   flux_predator = [flux_sz' flux_mz' flux_lz'];
   flux_N_predator = [flux_N_sz' flux_N_mz' flux_N_lz'];
   flux_P_predator = [flux_P_sz' flux_P_mz' flux_P_lz'];
   predation_bac = flux_sz(1) + flux_mz(1) + flux_lz(1);
   predation_N_bac = flux_N_sz(1) + flux_N_mz(1) + flux_N_lz(1);
   predation_P_bac = flux_P_sz(1) + flux_P_mz(1) + flux_P_lz(1);
   predation_p1 = flux_sz(4) + flux_mz(4) + flux_lz(4);
   predation_N_p1 = flux_N_sz(4) + flux_N_mz(4) + flux_N_lz(4);
   predation_P_p1 = flux_P_sz(4) + flux_P_mz(4) + flux_P_lz(4);
   predation_p2 = flux_sz(3) + flux_mz(3) + flux_lz(3);
   predation_N_p2 = flux_N_sz(3) + flux_N_mz(3) + flux_N_lz(3); 
   predation_P_p2 = flux_P_sz(3) + flux_P_mz(3) + flux_P_lz(3); 
   predation_p3 = flux_sz(2) + flux_mz(2) + flux_lz(2);
   predation_N_p3 = flux_N_sz(2) + flux_N_mz(2) + flux_N_lz(2); 
   predation_P_p3 = flux_P_sz(2) + flux_P_mz(2) + flux_P_lz(2);
   predation_p4 = flux_sz(5) + flux_mz(5) + flux_lz(5);
   predation_N_p4 = flux_N_sz(5) + flux_N_mz(5) + flux_N_lz(5); 
   predation_P_p4 = flux_P_sz(5) + flux_P_mz(5) + flux_P_lz(5);
   predation_sz = flux_sz(6) + flux_mz(6) + flux_lz(6);
   predation_N_sz = flux_N_sz(6) + flux_N_mz(6) + flux_N_lz(6); 
   predation_P_sz = flux_P_sz(6) + flux_P_mz(6) + flux_P_lz(6);
   predation_mz = flux_sz(7) + flux_mz(7) + flux_lz(7);
   predation_N_mz = flux_N_sz(7) + flux_N_mz(7) + flux_N_lz(7); 
   predation_P_mz = flux_P_sz(7) + flux_P_mz(7) + flux_P_lz(7);
   predation_lz = flux_sz(8) + flux_mz(8) + flux_lz(8);
   predation_N_lz = flux_N_sz(8) + flux_N_mz(8) + flux_N_lz(8); 
   predation_P_lz = flux_P_sz(8) + flux_P_mz(8) + flux_P_lz(8);         
   predation_mPOM = flux_sz(9) + flux_mz(9) + flux_lz(9);
   predation_N_mPOM = flux_N_sz(9) + flux_N_mz(9) + flux_N_lz(9);  
   predation_P_mPOM = flux_P_sz(9) + flux_P_mz(9) + flux_P_lz(9); 
   
% Flux adjustment for condervation of mass of meso-zooplankton for strict
% regulation of C:N:P, eq. (268-270)
% All of the released C, N, P goes to lPOM pool as in B16.
    netC_lz = uptake_lz - respiration_lz - excretion_lz - mortality_lz - predation_lz; % net C change
    netN_lz = uptake_N_lz - excretion_N_lz - mortality_N_lz - predation_N_lz; % net N change
    netP_lz = uptake_P_lz - excretion_P_lz - mortality_P_lz - predation_P_lz; % net P change
    
    % First adjust N and P (which ever is in more excess) relative to C when excess
%     relN_lz1 = nanmax(0,netN_lz-Q_N_RED.*netC_lz); % N released into PON pool if excess over C
%     relP_lz1 = nanmax(0,netP_lz-Q_P_RED.*netC_lz); % P released into POP pool if excess over C
    relN_lz1 = nanmax(0,netN_lz-Q_N_meso.*netC_lz); % N released into PON pool if excess over C
    relP_lz1 = nanmax(0,netP_lz-Q_P_meso.*netC_lz); % P released into POP pool if excess over C
%     nlim_plim_lz = netN_lz - (Q_N_RED./Q_P_RED).*netP_lz; % Positive = excess N, Negative = excess P
%     plim_nlim_lz = netP_lz - (Q_P_RED./Q_N_RED).*netN_lz; % Positive = excess P, Negative = excess N
    nlim_plim_lz = netN_lz - (Q_N_meso./Q_P_meso).*netP_lz; % Positive = excess N, Negative = excess P
    plim_nlim_lz = netP_lz - (Q_P_meso./Q_N_meso).*netN_lz; % Positive = excess P, Negative = excess N
    relN_lz1 = relN_lz1.*double(nlim_plim_lz > 0); % Release only when excess over P
    relP_lz1 = relP_lz1.*double(plim_nlim_lz > 0); % Release only when excess over N
    netN_temp_lz =  netN_lz - relN_lz1; % Net N change after first adjustment
    netP_temp_lz =  netP_lz - relP_lz1; % Net P change after first adjustment
    % Second adjust C to get right C:X when C:X > Redfield C:X 
%     clim_nutlim1 = nanmax(netC_lz-netP_temp_lz./Q_P_RED,netC_lz-netN_temp_lz./Q_N_RED); % Positive = excess C, Negative = excess nutrient
    clim_nutlim1 = nanmax(netC_lz-netP_temp_lz./Q_P_meso,netC_lz-netN_temp_lz./Q_N_meso); % Positive = excess C, Negative = excess nutrient
    relC_lz = clim_nutlim1.*double(clim_nutlim1>0); % first release of carbon relative to more limiting nutrient
    % Next adjust N:P balance
%     nlim_plim_lz2 = netN_temp_lz - (Q_N_RED./Q_P_RED).*netP_temp_lz; % Positive = excess N, Negative = excess P
%     plim_nlim_lz2 = netP_temp_lz - (Q_P_RED./Q_N_RED).*netN_temp_lz; % Positive = excess P, Negative = excess N
    nlim_plim_lz2 = netN_temp_lz - (Q_N_meso./Q_P_meso).*netP_temp_lz; % Positive = excess N, Negative = excess P
    plim_nlim_lz2 = netP_temp_lz - (Q_P_meso./Q_N_meso).*netN_temp_lz; % Positive = excess P, Negative = excess N
    relN_lz2 = nlim_plim_lz2.*double(nlim_plim_lz2 > 0); % release N if excess over P
    relP_lz2 = plim_nlim_lz2.*double(plim_nlim_lz2 > 0); % release P if excess over N
    netN_temp_lz2 =  netN_temp_lz - relN_lz2; % Net N change after second adjustment
    netP_temp_lz2 =  netP_temp_lz - relP_lz2; % Net P change after second adjustment

    % Total Release of N and P after adjustments 
    relN_lz = relN_lz1 + relN_lz2;
    relP_lz = relP_lz1 + relP_lz2;
    % Check : Should Equal to zero -> Checked
%      (netC_lz - relC_lz)-(netN_lz - relN_lz)./Q_N_RED;
%      (netC_lz - relC_lz)-(netP_lz - relP_lz)./Q_P_RED;
%      (netN_lz - relN_lz) - (Q_N_RED./Q_P_RED).*(netP_lz - relP_lz);
     (netC_lz - relC_lz)-(netN_lz - relN_lz)./Q_N_meso;
     (netC_lz - relC_lz)-(netP_lz - relP_lz)./Q_P_meso;
     (netN_lz - relN_lz) - (Q_N_meso./Q_P_meso).*(netP_lz - relP_lz);

% Bacteria functions X(27)-X(29)
    % Bacteria quota
    Q_N_bac = X(28)./X(27);
    Q_N_bac_prime = (X(28)-epsilon_n)./X(27);
    Q_P_bac = X(29)./X(27);
    Q_P_bac_prime = (X(29)-epsilon_p)./X(27);
    % Uptake of DOC
    n_bac = nanmin((X(2)+X(40))./(X(2)+X(40)+hn_bac),(X(3)+X(41))./(X(3)+X(41)+hp_bac)); % nutrient scarcity eq. (29)
    eO2_bac = O_rel./(O_rel + h_ox_bac); % oxygen limitation eq. (28)
    uptake_DOC_bac = nanmin(r_ass_bac.*temp.*eO2_bac.*n_bac.*X(27),(X(39)-epsilon)); % uptake of DOC eq. (27)
    uptake_DON_bac = uptake_DOC_bac.*(X(40)-epsilon_n)./X(39);
    uptake_DOP_bac = uptake_DOC_bac.*(X(41)-epsilon_n)./X(39);
    
    % Uptake and loss of N, P from inorgnic pool based on internal
    % nutrient 
    excessN_bac = Q_N_bac - QMAX_N_BAC; 
    excessP_bac = Q_P_bac - QMAX_P_BAC; 
    
    uptake_DIN_bac = double(excessN_bac < 0).*-1.*excessN_bac.*X(27).*(X(2)-epsilon_n)./((X(2)-epsilon_n) + hn_bac);
    uptake_DIP_bac = double(excessP_bac < 0).*-1.*excessP_bac.*X(27).*(X(3)-epsilon_p)./((X(3)-epsilon_p) + hp_bac);
    loss_N_bac = double(excessN_bac >= 0).*excessN_bac.*X(27);
    loss_P_bac = double(excessP_bac >= 0).*excessP_bac.*X(27);
    
    % Loss of C due to respiration = activity respiration + rest
    % respiration
    respiration_bac = uptake_DOC_bac.*(r_resp_bac.*O_rel + r_respox_bac.*(1.0-O_rel)) + ...
        r_basal_bac.*temp.*(X(27)-epsilon); % eq. (30)
    % Loss due to Mortality
    mortality_bac = r_mort_bac.*temp.*(X(27)-epsilon); % eq. (31)
    mortality_N_bac = mortality_bac.*Q_N_bac_prime;
    mortality_P_bac = mortality_bac.*Q_P_bac_prime;
    
% recycling

% small POM processes (X(30)-X(32))
% 1. small POM: derives from flagellates (p2), pico-plankton (p3), and heterotrophic
% flagellates (sz)
  % Gain due to lysis of phytoplankton p2 and p3  
  % qdloss = fraction allocated to POM pool over DOM pool
  qdloss_p2 = nanmin(QMIN_N_PHYTO(2)./Q_N_phyto(2),QMIN_P_PHYTO(2)./Q_P_phyto(2));
  qdloss_p3 = nanmin(QMIN_N_PHYTO(3)./Q_N_phyto(3),QMIN_P_PHYTO(3)./Q_P_phyto(3));
  lysis_sPOC = qdloss_p2.*lysis(2) + qdloss_p3.*lysis(3);  
  lysis_sPON = qdloss_p2.*lysis_N(2) + qdloss_p3.*lysis_N(3);
  lysis_sPOP = qdloss_p2.*lysis_P(2) + qdloss_p3.*lysis_P(3);

  % Gain due to mortality of sz
  mortality_sz_sPOC = (1.0-P_DOM(1)).*mortality_sz; 
  mortality_sz_sPON = (1.0-P_DOM(1)).*mortality_N_sz; 
  mortality_sz_sPOP = (1.0-P_DOM(1)).*mortality_P_sz;
  
  mortality_sPOC = mortality_sz_sPOC;
  mortality_sPON = mortality_sz_sPON;
  mortality_sPOP = mortality_sz_sPOP;

  % Gain due to excretion of sz
  excretion_sPOC = (1.0-P_DOM(1)).*excretion_sz; 
  excretion_sPON = (1.0-P_DOM(1)).*excretion_N_sz; 
  excretion_sPOP = (1.0-P_DOM(1)).*excretion_P_sz;
 
  % Loss due to breakdown into DOM
  breakdown_sPOC = X(31).*(rfd_cn.*12).*Ra; % break down of sPOC eq. (35)
  breakdown_sPON = breakdown_sPOC.*(X(31)-epsilon_n)./X(30);
  breakdown_sPOP = breakdown_sPOC.*(X(32)-epsilon_p)./X(30);
  
  % Loss due to sinking
  sinking_sPOC = V_SINK(1).*(X(30)-epsilon); 
  sinking_sPON = V_SINK(1).*(X(31)-epsilon_n); 
  sinking_sPOP = V_SINK(1).*(X(32)-epsilon_p); 
  
% medium POM processes (X(33)-X(35))
% 2. medium POM: derives from diatoms (p1), dino-flagellates (p4), and microzooplankton (mz)
  % Gain due to lysis of phytoplankton p1 and p4  
  % qdloss = fraction allocated to POM pool over DOM pool
  qdloss_p1 = nanmin(QMIN_N_PHYTO(1)./Q_N_phyto(1),QMIN_P_PHYTO(1)./Q_P_phyto(1));
  qdloss_p4 = nanmin(QMIN_N_PHYTO(4)./Q_N_phyto(4),QMIN_P_PHYTO(4)./Q_P_phyto(4));
  lysis_mPOC = qdloss_p1.*lysis(1) + qdloss_p4.*lysis(4);  
  lysis_mPON = qdloss_p1.*lysis_N(1) + qdloss_p4.*lysis_N(4);
  lysis_mPOP = qdloss_p1.*lysis_P(1) + qdloss_p4.*lysis_P(4);
  % Gain due to mortality of mz
  mortality_mz_mPOC = (1.0-P_DOM(2)).*mortality_mz; 
  mortality_mz_mPON = (1.0-P_DOM(2)).*mortality_N_mz; 
  mortality_mz_mPOP = (1.0-P_DOM(2)).*mortality_P_mz;
  
  mortality_mPOC = mortality_mz_mPOC;
  mortality_mPON = mortality_mz_mPON;
  mortality_mPOP = mortality_mz_mPOP;
  
  % Gain due to excretion of mz
  excretion_mPOC = (1.0-P_DOM(2)).*excretion_mz; 
  excretion_mPON = (1.0-P_DOM(2)).*excretion_N_mz; 
  excretion_mPOP = (1.0-P_DOM(2)).*excretion_P_mz;
  
  % Loss due to breakdown into DOM
  breakdown_mPOC = X(34).*(rfd_cn.*12).*Ra; % break down of mPOC eq. (35)
  breakdown_mPON = breakdown_mPOC.*(X(34)-epsilon_n)./X(33);
  breakdown_mPOP = breakdown_mPOC.*(X(35)-epsilon_p)./X(33);
  
  % Loss due to sinking
  sinking_mPOC = V_SINK(2).*(X(33)-epsilon); 
  sinking_mPON = V_SINK(2).*(X(34)-epsilon_n); 
  sinking_mPOP = V_SINK(2).*(X(35)-epsilon_p); 
% large POM processes (X(36)-X(38))
% 3. large POM: derives from  mesozooplankton (lz)
  % Gain due to lysis of phytoplankton p1 and p4  
  % Gain due to mortality of lz
  mortality_lz_lPOC = (1.0-P_DOM(3)).*mortality_lz; 
  mortality_lz_lPON = (1.0-P_DOM(3)).*mortality_N_lz; 
  mortality_lz_lPOP = (1.0-P_DOM(3)).*mortality_P_lz;
  
  mortality_lPOC = mortality_lz_lPOC;
  mortality_lPON = mortality_lz_lPON;
  mortality_lPOP = mortality_lz_lPOP;
  % Gain due to excretion of lz and release of excess C,N,P adjustment
  excretion_lPOC = (1.0-P_DOM(3)).*excretion_lz + relC_lz; 
  excretion_lPON = (1.0-P_DOM(3)).*excretion_N_lz + relN_lz; 
  excretion_lPOP = (1.0-P_DOM(3)).*excretion_P_lz + relP_lz;

  % Loss due to breakdown into DOM
  breakdown_lPOC = X(37).*(rfd_cn.*12).*Ra; % break down of lPOC eq. (35)
  breakdown_lPON = breakdown_lPOC.*(X(37)-epsilon_n)./X(36);
  breakdown_lPOP = breakdown_lPOC.*(X(38)-epsilon_p)./X(36);

  % Loss due to sinking
  sinking_lPOC = V_SINK(3).*(X(36)-epsilon); 
  sinking_lPON = V_SINK(3).*(X(37)-epsilon_n); 
  sinking_lPOP = V_SINK(3).*(X(38)-epsilon_p); 

% Labile DOM Functions
    % Gain due to lysis of phytoplankton
    lysis_OC = (1.0-qdloss_p1).*lysis(1) + (1.0-qdloss_p2).*lysis(2) + ...
        (1.0-qdloss_p3).*lysis(3) + (1.0-qdloss_p4).*lysis(4); 
    lysis_DOC = q_dis.*lysis_OC;
    lysis_ON = (1.0-qdloss_p1).*lysis_N(1) + (1.0-qdloss_p2).*lysis_N(2) + ...
        (1.0-qdloss_p3).*lysis_N(3) + (1.0-qdloss_p4).*lysis_N(4); 
    lysis_DON = q_dis.*lysis_ON;
    lysis_OP = (1.0-qdloss_p1).*lysis_P(1) + (1.0-qdloss_p2).*lysis_P(2) + ...
        (1.0-qdloss_p3).*lysis_P(3) + (1.0-qdloss_p4).*lysis_P(4);
    lysis_DOP = q_dis.*lysis_OP;
    % Gain due to excretion of phytoplankton and zooplankton 
    excretion_OC = nansum(excretion_C_phyto) + P_DOM(1).*excretion_sz + P_DOM(2).*excretion_mz + P_DOM(3).*excretion_lz;
    excretion_ON = nansum(excretion_N_phyto) + P_DOM(1).*excretion_N_sz + P_DOM(2).*excretion_N_mz + P_DOM(3).*excretion_N_lz;
    excretion_OP = nansum(excretion_P_phyto) + P_DOM(1).*excretion_P_sz + P_DOM(2).*excretion_P_mz + P_DOM(3).*excretion_P_lz;
    excretion_DOC = q_dis.*excretion_OC;
    excretion_DON = q_dis.*excretion_ON;
    excretion_DOP = q_dis.*excretion_OP;
    % Gain due to mortality of zooplankton and bacteria
    mortality_OC = P_DOM(1).*mortality_sz + P_DOM(2).*mortality_mz + P_DOM(3).*mortality_lz + mortality_bac;
    mortality_ON = P_DOM(1).*mortality_N_sz + P_DOM(2).*mortality_N_mz + P_DOM(3).*mortality_N_lz + mortality_N_bac;
    mortality_OP = P_DOM(1).*mortality_P_sz + P_DOM(2).*mortality_P_mz + P_DOM(3).*mortality_P_lz + mortality_P_bac;
    mortality_DOC = q_dis.*mortality_OC;
    mortality_DON = q_dis.*mortality_ON;
    mortality_DOP = q_dis.*mortality_OP;
    % Gain due to break down of POM
    gain_POC = breakdown_sPOC + breakdown_mPOC + breakdown_lPOC;
    gain_PON = breakdown_sPON + breakdown_mPON + breakdown_lPON;
    gain_POP = breakdown_sPOP + breakdown_mPOP + breakdown_lPOP;
    % Loss due to break down of semi-labile DOM
    decomp_slDOC = Ra.*(X(44)-epsilon);
    decomp_slDON = Ra.*(X(45)-epsilon_n);
    decomp_slDOP = Ra.*(X(46)-epsilon_p);
    % Loss of N and P due to remineralization
    remin_DON = r_remin.*(X(40)-epsilon_n); 
    remin_DOP = r_remin.*(X(41)-epsilon_p); 
    
% semi-Labile DOM Functions
    % Gain due to lysis of phytoplankton
    lysis_slDOC = (1-q_dis).*lysis_OC;
    lysis_slDON = (1-q_dis).*lysis_ON;
    lysis_slDOP = (1-q_dis).*lysis_OP;
    % Gain due to excretion of phytoplankton and zooplankton 
    excretion_slDOC = (1-q_dis).*excretion_OC;
    excretion_slDON = (1-q_dis).*excretion_ON;
    excretion_slDOP = (1-q_dis).*excretion_OP;
    % Gain due to mortality of zooplankton and bacteria
    mortality_slDOC = (1-q_dis).*mortality_OC;
    mortality_slDON = (1-q_dis).*mortality_ON;
    mortality_slDOP = (1-q_dis).*mortality_OP;

    
% Inorganic Pool Soource  terms     
% Right now, assume everything goes back to inorganic pool

recycling_N = loss_N_sz + loss_N_mz + loss_N_bac + remin_DON; %
recycling_P = loss_P_sz + loss_P_mz + loss_P_bac + remin_DOP; %
% 
% Nitrification
nitrification = r_nitr.*temp.*(X(2)-epsilon_n); % nitrification eq. (no sediments)
% 
% N uptake fractionation
fno3 = (A_N.*(X(1)-epsilon_n))./(A_N.*(X(1)-epsilon_n) + A_A.*(X(2)-epsilon_n)); % fractionation for no3

% ----------------------------------------------------------------------
% differential equations
dXdt(1) = A.*(Nin_N-X(1)-epsilon_n)-nansum(fno3.*uptake_N) + nitrification;   % dNn/dt eq. (153)
dXdt(2) = A.*(Nin_A-X(2)-epsilon_n)-nansum((1-fno3).*uptake_N) - uptake_DIN_bac - nitrification + recycling_N;   % dNa/dt eq. (154)
dXdt(3) = A.*(Nin_P-X(3)-epsilon_p)-nansum(uptake_P) - uptake_DIP_bac + recycling_P;   % dNp/dt eq. (155)
dXdt(4) = photosynthesis(1) - respiration_phyto(1) - excretion_C_phyto(1) - lysis(1) - sinking_phyto(1) - predation_p1; % dPc/dt eq. (129)
dXdt(5) = photosynthesis(1).*(phi_phyto(1)) - (respiration_phyto(1) + excretion_C_phyto(1) + lysis(1) + sinking_phyto(1) + predation_p1).*theta_phyto_prime(1); % dPchl/dt eq. (130)
dXdt(6) = uptake_N(1) - excretion_N_phyto(1) -lysis_N(1) - sinking_phyto(1).*Q_N_phyto_prime(1) - predation_N_p1; % dPn/dt eq. (131)
dXdt(7) = uptake_P(1) - excretion_P_phyto(1) -lysis_P(1) - sinking_phyto(1).*Q_P_phyto_prime(1) - predation_P_p1; % dPp/dt eq. (132)
dXdt(8) = photosynthesis(2) - respiration_phyto(2) - excretion_C_phyto(2) - lysis(2) - sinking_phyto(2) - predation_p2; % dPc/dt eq. (129)
dXdt(9) = photosynthesis(2).*(phi_phyto(2)) - (respiration_phyto(2) + excretion_C_phyto(2) + lysis(2) + sinking_phyto(2) + predation_p2).*theta_phyto_prime(2); % dPchl/dt eq. (130)
dXdt(10) = uptake_N(2) - excretion_N_phyto(2) -lysis_N(2) - sinking_phyto(2).*Q_N_phyto_prime(2) - predation_N_p2; % dPn/dt eq. (131)
dXdt(11) = uptake_P(2) - excretion_P_phyto(2) -lysis_P(2) - sinking_phyto(2).*Q_P_phyto_prime(2) - predation_P_p2; % dPp/dt eq. (132)
dXdt(12) = photosynthesis(3) - respiration_phyto(3) - excretion_C_phyto(3) - lysis(3) - sinking_phyto(3) - predation_p3; % dPc/dt eq. (129)
dXdt(13) = photosynthesis(3).*(phi_phyto(3)) - (respiration_phyto(3) + excretion_C_phyto(3) + lysis(3) + sinking_phyto(3) + predation_p3).*theta_phyto_prime(3); % dPchl/dt eq. (130)
dXdt(14) = uptake_N(3) - excretion_N_phyto(3) -lysis_N(3) - sinking_phyto(3).*Q_N_phyto_prime(3) - predation_N_p3; % dPn/dt eq. (131)
dXdt(15) = uptake_P(3) - excretion_P_phyto(3) -lysis_P(3) - sinking_phyto(3).*Q_P_phyto_prime(3) - predation_P_p3; % dPp/dt eq. (132)
dXdt(16) = photosynthesis(4) - respiration_phyto(4) - excretion_C_phyto(4) - lysis(4) - sinking_phyto(4) - predation_p4; % dPc/dt eq. (129)
dXdt(17) = photosynthesis(4).*(phi_phyto(4)) - (respiration_phyto(4) + excretion_C_phyto(4) + lysis(4) + sinking_phyto(4) + predation_p4).*theta_phyto_prime(4); % dPchl/dt eq. (130)
dXdt(18) = uptake_N(4) - excretion_N_phyto(4) -lysis_N(4) - sinking_phyto(4).*Q_N_phyto_prime(4) - predation_N_p4; % dPn/dt eq. (131)
dXdt(19) = uptake_P(4) - excretion_P_phyto(4) -lysis_P(4) - sinking_phyto(4).*Q_P_phyto_prime(4) - predation_P_p4; % dPp/dt eq. (132)
dXdt(20) = uptake_sz - respiration_sz - excretion_sz - mortality_sz - predation_sz; %dZc1/dt
dXdt(21) = uptake_N_sz - excretion_N_sz - mortality_N_sz - predation_N_sz - loss_N_sz; %dZn1/dt 
dXdt(22) = uptake_P_sz - excretion_P_sz - mortality_P_sz - predation_P_sz - loss_P_sz; %dZp1/dt 
dXdt(23) = uptake_mz - respiration_mz - excretion_mz - mortality_mz - predation_mz; %dZc2/dt
dXdt(24) = uptake_N_mz - excretion_N_mz - mortality_N_mz - predation_N_mz - loss_N_mz; %dZn2/dt 
dXdt(25) = uptake_P_mz - excretion_P_mz - mortality_P_mz - predation_P_mz - loss_P_mz; %dZp2/dt 
dXdt(26) = netC_lz - relC_lz; %dZc3/dt
dXdt(27) = uptake_DOC_bac - respiration_bac - mortality_bac - predation_bac; % dBc/dt
dXdt(28) = uptake_DON_bac + uptake_DIN_bac - mortality_N_bac - predation_N_bac - loss_N_bac; % dBn/dt 
dXdt(29) = uptake_DOP_bac + uptake_DIP_bac - mortality_P_bac - predation_P_bac - loss_P_bac; % dBp/dt
dXdt(30) = lysis_sPOC + mortality_sPOC + excretion_sPOC - breakdown_sPOC -sinking_sPOC; % dsPOC/dt
dXdt(31) = lysis_sPON + mortality_sPON + excretion_sPON - breakdown_sPON -sinking_sPON; % dsPON/dt
dXdt(32) = lysis_sPOP + mortality_sPOP + excretion_sPOP - breakdown_sPOP -sinking_sPOP; % dsPOP/dt
dXdt(33) = lysis_mPOC  + mortality_mPOC + excretion_mPOC - breakdown_mPOC -sinking_mPOC - predation_mPOM; % dmPOC/dt
dXdt(34) = lysis_mPON  + mortality_mPON + excretion_mPON - breakdown_mPON -sinking_mPON - predation_N_mPOM; % dmPON/dt
dXdt(35) = lysis_mPOP  + mortality_mPOP + excretion_mPOP - breakdown_mPOP -sinking_mPOP - predation_P_mPOM; % dmPOP/dt
dXdt(36) = mortality_lPOC + excretion_lPOC - breakdown_lPOC -sinking_lPOC; % dlPOC/dt
dXdt(37) = mortality_lPON + excretion_lPON - breakdown_lPON -sinking_lPON; % dlPON/dt
dXdt(38) = mortality_lPOP + excretion_lPOP - breakdown_lPOP -sinking_lPOP; % dlPOP/dt
dXdt(39) = lysis_DOC + excretion_DOC + mortality_DOC + gain_POC + decomp_slDOC- uptake_DOC_bac; % dDOC/dt
dXdt(40) = lysis_DON + excretion_DON + mortality_DON + gain_PON + decomp_slDON- uptake_DON_bac - remin_DON; % dDON/dt
dXdt(41) = lysis_DOP + excretion_DOP + mortality_DOP + gain_POP + decomp_slDOP- uptake_DOP_bac - remin_DOP; % dDOP/dt
% dXdt(42) = netN_temp_lz2; %dZn3/dt
% dXdt(43) = netP_temp_lz2; %dZp3/dt
dXdt(42) = netN_lz-relN_lz; %dZn3/dt
dXdt(43) = netP_lz-relP_lz; %dZp3/dt
dXdt(44) = lysis_slDOC + excretion_slDOC + mortality_slDOC - decomp_slDOC; %dslDOC/dt
dXdt(45) = lysis_slDON + excretion_slDON + mortality_slDON - decomp_slDON; %dslDON/dt
dXdt(46) = lysis_slDOP + excretion_slDOP + mortality_slDOP - decomp_slDOP; %dslDOP/dt
dXdt(47) = A.*(Nin_N+Nin_A);
dXdt(48) = nansum(sinking_phyto.*Q_N_phyto) + sinking_sPON + sinking_mPON + sinking_lPON + A.*(X(1)+X(2));
dXdt(49) = A.*(Nin_P);
dXdt(50) = nansum(sinking_phyto.*Q_P_phyto) + sinking_sPOP + sinking_mPOP + sinking_lPOP + A.*(X(3));