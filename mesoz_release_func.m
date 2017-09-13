% mesoz_release_func.m
% Function for calculating released C, N, and P (in mmol m^-3) for Mesozooplankton
% Released as a result of internal stoichiometric adjustment
% Tatsuro Tanioka 17/08/30
function loss_nut = mesoz_release_func(X,age)
global epsilon %[mgC m^-3] Negativity control factor = 0.01 (Butenschön et al., 2016)
global epsilon_p %[mgC m^-3] Negativity control factor for P = 0.01 * (P:C)redfield (Butenschön et al., 2016)
global epsilon_n %[mgC m^-3] Negativity control factor for N = 0.01 * (N:C)redfield (Butenschön et al., 2016)
global temperature % [degC]
global Q10 % 
global O_rel %
global Q_P_RED
global Q_N_RED
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
% -------------------------------------------------------------------
% Fixed temperature dependence 
t_water = temperature(1);
% n = round(t); %round t to the next integer(t~=n);
temp = Q10.^((t_water-10)./10) - Q10.^((t_water-32)./3); % temperature function eq. (239), B16 
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
% Calculate Pf_flex
Pf_flex = Pf_flex_func(X,sigma_pf,Pf_var);

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
        Q_P_meso = 1.0/(COPE_AGE_cpslope.*age + COPE_AGE_cpintercept)./12.*double(COPE_AGE_switch > 0) + Q_P_RED.*double(COPE_AGE_switch < 1) ; % P:C [mol P/gC]according to Meunier16
        Q_NP_meso = (COPE_AGE_npslope.*age + COPE_AGE_npintercept).*double(COPE_AGE_switch > 0) + (Q_N_RED./Q_P_RED).*double(COPE_AGE_switch < 1); % N:P [mol N/ mol P] according to Meunier16
        Q_N_meso = Q_NP_meso.*Q_P_meso; % Calcualte N:C [mol N/gC]
     % Calculating max. assimilation rate as a function of age
        R_ASS_meso = (COPE_AGE_rmaxslope*age + COPE_AGE_rmaxintercept).*double(COPE_AGE_switch > 0) + R_ASS_ZOO(3).*double(COPE_AGE_switch < 1); % max.assimilation rate according to Meunier16
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
    
    % Total released C, N, P:
    loss_C = relC_lz./12;
    %loss_N = excretion_N_lz + relN_lz;
    %loss_P = excretion_P_lz + relP_lz;
    loss_N = relN_lz;
    loss_P = relP_lz;
    loss_nut = [loss_C loss_N loss_P];
    
    % Check: these should equal zero, if calculate correctly -> checked
     (netC_lz - relC_lz)-(netN_lz - relN_lz)./Q_N_meso;
     (netC_lz - relC_lz)-(netP_lz - relP_lz)./Q_P_meso;
     (netN_lz - relN_lz) - (Q_N_meso./Q_P_meso).*(netP_lz - relP_lz);


