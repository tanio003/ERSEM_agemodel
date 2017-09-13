% Function for calculating Pf matrix as a function of sigma
% updated for flex_func7
function Pf_flex = Pf_flex_func(X,sigma,pf_var)
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
      x = ((Q_N_Prey'./Q_P_Prey') - ones(9,3).*Q_NP_optfood)./(ones(9,3).*Q_NP_optfood); % (N:P(prey) - N:P(optimal))/N:P(optimal)
      fq = normpdf(x,0.0,sigma); % food quality
      
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
  
      if pf_var == 0;
        Pf_flex = Pf;
      elseif pf_var == 1;
        Pf_flex = fq.*food_frac; % Food preference = food quality + Fasham
      elseif pf_var == 2;
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
      Pf_flex = Pf_flex.*double(pf_var>0) + Pf_fix.*double(pf_var<=0);  
      Pf_flex(:,1) =  Pf_fix(:,1); % small z has fixed preference
      Pf_flex(:,2) =  Pf_fix(:,2); % micro z has fixed preference
      %Pf_flex(:,3) =  Pf_fix(:,3);

      Pf_flex;

