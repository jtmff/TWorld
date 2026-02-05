%     Cardiac model T-World
%     Copyright (C) 2026 Jakub Tomek. Contact: jakub.tomek.mff@gmail.com
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function [output] = model_TWorld(t,y,flag_ode, cellType, ICaL_Multiplier, ...
    INa_Multiplier, Itos_Multiplier, Itof_Multiplier, INaL_Multiplier, IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKb_Multiplier,INaCa_Multiplier,...
    INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, IClCa_Multiplier, IClb_Multiplier, Jrel_Multiplier,Jup_Multiplier, ions_na_o,ions_ca_o,ions_k_o,ICaL_fraction_junc,INaCa_fraction_junc, PKA_P, Whole_cell_PP1, runSignallingPathway, constantsSig, sex, stimAmp, stimDur, vcParameters, apClamp, extraParams)


%% Handling sex-specific models:
% the parameter sex can be used to set male/female model
CMDN_Multiplier = 1;
if strcmp(sex, 'M')

    INaCa_Multiplier = INaCa_Multiplier / 1.0724;
    IKr_Multiplier = IKr_Multiplier * 1.1;
    IKs_Multiplier = IKs_Multiplier * 1.09;
    IK1_Multiplier = IK1_Multiplier * 1.07;
    IpCa_Multiplier = IpCa_Multiplier * 0.8; 
    CMDN_Multiplier = 1/1.1;

elseif strcmp(sex, 'F')
    INaCa_Multiplier = INaCa_Multiplier * 1.0724;
    IKr_Multiplier = IKr_Multiplier * 0.87;
    IKs_Multiplier = IKs_Multiplier * 0.905;
    IK1_Multiplier = IK1_Multiplier * 0.92;
    IpCa_Multiplier = IpCa_Multiplier * 1.28;
    CMDN_Multiplier = 1.1;
end


%% PKA signalling:
% simulating or using direct phosphorylation values
if runSignallingPathway == 1 % if signalling is simulated
    ydot_signalling = getPKASignalling(y(93+(1:57)), constantsSig);
    [fICaL_PKA,fIKs_PKA,fPLB_PKA,fTnI_PKA,fINa_PKA,fINaK_PKA,fRyR_PKA,fIKb_PKA]=  getEffectiveFraction(y(93+(1:57)), constantsSig);
    % RyR and I_Kb phosphorylation are currently not used due to lack of
    % convincing evidence, but are kept for future use if there is a good
    % argument for their effects.

    fMyBPC_PKA = fTnI_PKA;
    % Concentration of uninhibited PP1 in the cytosolic compartment
    pp1_PP1f_cyt_sum = constantsSig(38) - constantsSig(37) + y(93+39);
    % Concentration of uninhibited PP1 in the cytosolic compartment
    PP1f_cyt = 0.5 * (sqrt(pp1_PP1f_cyt_sum ^ 2.0 + 4.0 * constantsSig(38) * constantsSig(37)) - pp1_PP1f_cyt_sum);
    Whole_cell_PP1  = constantsSig(36) / constantsSig(5) + constantsSig(35) / constantsSig(6) + PP1f_cyt / constantsSig(7);
else % or phosphorylation levels can be directly set
    fINa_PKA = PKA_P(1);
    fICaL_PKA = PKA_P(2);
    fINaK_PKA = PKA_P(3);
    fIKs_PKA = PKA_P(4);
    fPLB_PKA = PKA_P(5);
    fTnI_PKA = PKA_P(6);
    fMyBPC_PKA = PKA_P(7);
    % Whole_cell_PP1 is provided as an input in this case.
end

%% State variables are extracted and named
membrane_v = y(1);

ions_na_junc = y(2);
ions_na_sl = y(3);
ions_na_i = y(4);
ions_k_i = y(5);
ions_ca_junc = y(6);
ions_ca_sl = y(7);
ions_ca_i = y(8);
ions_cl_i = y(9);
ions_ca_SR = y(10);

buffers_NaBj = y(11);
buffers_NaBsl = y(12);
buffers_TnClow = y(13); 
buffers_TnCHc = y(14);
buffers_TnCHm = y(15);
buffers_CaM = y(16);
buffers_Myosin_ca = y(17);
buffers_Myosin_mg = y(18);
buffers_SRB = y(19);
buffers_SLLj = y(20);
buffers_SLLsl = y(21);
buffers_SLHj = y(22);
buffers_SLHsl = y(23);
buffers_Csqn = y(24);

ina_m = y(25);
ina_h = y(26);
ina_j = y(27);
inal_ml = y(28);
inal_hl = y(29);
inal_hl_p = y(30);

ito_xtos = y(31);
ito_ytos = y(32);
ito_xtof = y(33);
ito_ytof = y(34);
ito_xtos_p = y(35);
ito_ytos_p = y(36);
ito_xtof_p = y(37);
ito_ytof_p = y(38);

iks_xs_junc = y(39);
iks_xs_sl = y(40);

ikr_c0 = y(41);
ikr_c1 = y(42);
ikr_c2 = y(43);
ikr_o = y(44);
ikr_i = y(45);

ical_d=y(46);
ical_ff=y(47);
ical_fs=y(48);
ical_fcaf=y(49);
ical_fcas=y(50);
ical_jca=y(51);
ical_nca=y(52);
ical_nca_i=y(53);
ical_ffp=y(54);
ical_fcafp=y(55);
ical_pureCDI_junc = y(56);
ical_pureCDI_sl = y(57);

ryr_R = y(58);
ryr_O = y(59);
ryr_I = y(60);
ryr_CaRI = y(61);
ryr_R_p = y(62);
ryr_O_p = y(63);
ryr_I_p = y(64);
ryr_CaRI_p = y(65);
jrel_icaldep_act = y(66);
jrel_icaldep_f1 = y(67);
jrel_icaldep_f2 = y(68);

contraction_XS=max(0,y(69));
contraction_XW=max(0,y(70));
contraction_Ca_TRPN=max(0,y(71));
contraction_TmBlocked=y(72);
contraction_ZETAS=y(73);
contraction_ZETAW=y(74);

camk_trap = y(75);
camk_f_ICaL = y(76);
camk_f_RyR = y(77);
camk_f_PLB = y(78);
casig_serca_trap = y(79);

ina_hp = y(80);
ina_jp = y(81);
ina_m_PKA = y(82);
ina_h_PKA = y(83);
ina_j_PKA = y(84);
ina_h_dualP = y(85);
ina_j_dualP = y(86);
ical_d_P = y(87);
ical_ff_P = y(88);
ical_fs_P = y(89);
ical_fcaf_P = y(90);
ical_fcas_P = y(91);
ical_fBPf = y(92);
ical_fcaBPf = y(93);

% 94-150 are signalling constants

%% Model Parameters
R = 8314;       
Frdy = 96485;   
Temp = 310;     
FoRT = Frdy/R/Temp;
Cmem = 1.3810e-10;   
Qpow = (Temp - 310)/10;
F = Frdy;
T = Temp;
vfrt = membrane_v * F/(R*T);

% Cell dimensions and compartment volumes
cellLength = 100;     % cell length in um
cellRadius = 10.25;   
Vcell = pi*cellRadius^2*cellLength*1e-15;  
Vmyo = 0.65*Vcell; 
Vsr = 0.035*Vcell; 
Vsl = 0.02*Vcell; 
Vjunc = 0.0539*0.01*Vcell;

% Diffusion coefficients between compartments
J_ca_juncsl = 1/3.06685e12;
J_ca_slmyo =  1/0.74556e11;
J_na_juncsl = 1/(1.6382e12/3*100);
J_na_slmyo = 1/(1.8308e10/3*100);

% Fractional distribution of ionic currents between compartments
Fjunc = 0.11; Fsl = 1-Fjunc;
% different localizations are used for ICaL (ICaL_fractionSS), NCX (INaCa_fractionSS), and I(Ca)Cl (50:50)

% Fixed ion concentrations (extracellular Cl and intracellular Mg)
ions_cl_o = 150;  
ions_mg_i = 0.5;  

% Nernst Potentials
ena_junc = (1/FoRT)*log(ions_na_o/ions_na_junc);    
ena_sl = (1/FoRT)*log(ions_na_o/ions_na_sl);       
ek = (1/FoRT)*log(ions_k_o/ions_k_i);	        
eca_junc = (1/FoRT/2)*log(ions_ca_o/ions_ca_junc); 
eca_sl = (1/FoRT/2)*log(ions_ca_o/ions_ca_sl);     
ecl = (1/FoRT)*log(ions_cl_i/ions_cl_o);           


%% CaMK and Ca signalling
[d_camk_trap, d_camk_f_ICaL, d_camk_f_RyR, d_camk_f_PLB, d_casig_serca_trap, casig_SERCA_act] = getCaSignalling(ions_ca_junc, camk_trap, camk_f_ICaL, camk_f_RyR, camk_f_PLB, casig_serca_trap, Whole_cell_PP1);

%% I_Na and I_NaL (fast and late Na current)
fINap = camk_f_RyR; % CaMKII phosphorylation is currently assumed to be that of RyR

[I_NaFast_junc, I_NaFast_sl, d_ina_m, d_ina_h, d_ina_hp, d_ina_j, d_ina_jp, d_ina_m_PKA, d_ina_h_PKA, d_ina_j_PKA, d_ina_h_dualP, d_ina_j_dualP] = getINa_Grandi(membrane_v, ina_m, ina_h, ina_hp, ina_j, ina_jp, fINap, ena_junc, ena_sl, INa_Multiplier,fINa_PKA, Fjunc,ina_m_PKA, ina_h_PKA, ina_j_PKA, ina_h_dualP, ina_j_dualP);

[I_NaL_junc, I_NaL_sl,d_inal_ml,d_inal_hl, d_inal_hl_p] = getINaL_ORd2011(membrane_v, inal_ml, inal_hl, inal_hl_p, fINap, ena_junc, ena_sl, cellType, Fjunc, Fsl, INaL_Multiplier);

% Na current aggregation variables:
% 1) for ionic update purposes, total sodium junc and sl are calculated
I_Na_junc = I_NaFast_junc + I_NaL_junc ;
I_Na_sl = I_NaFast_sl + I_NaL_sl ;

% 2) for plotting purposes, storing total fast and late sodium current
I_NaFast = I_NaFast_junc + I_NaFast_sl;
I_NaL = I_NaL_junc + I_NaL_sl;

%% I_Nabk (Background sodium current)
GNaB = 0.594e-3 * INab_Multiplier;    % [mS/uF]
I_Nabk_junc = Fjunc*GNaB*(membrane_v-ena_junc);
I_Nabk_sl = Fsl*GNaB*(membrane_v-ena_sl);
I_Nabk = I_Nabk_junc+I_Nabk_sl;

%% I_NaK: Na/K Pump Current
IbarNaK = 2.10774 * INaK_Multiplier;     % [uA/uF]
KmNaip = 11;         
KmNaip_PKA = 8.4615;
KmKo = 1.5;         

fnak = 0.75 + (0.00375- ((140 - ions_na_o)/50) * 0.001)*membrane_v ; % Voltage and nao-sensitivity
I_NaK_junc_noPKA = Fjunc*IbarNaK*fnak*ions_k_o /(1+(KmNaip/ions_na_junc)^4) /(ions_k_o+KmKo);
I_NaK_sl_noPKA = Fsl*IbarNaK*fnak*ions_k_o /(1+(KmNaip/ions_na_sl)^4) /(ions_k_o+KmKo);
I_NaK_junc_PKA = Fjunc*IbarNaK*fnak*ions_k_o /(1+(KmNaip_PKA/ions_na_junc)^4) /(ions_k_o+KmKo);
I_NaK_sl_PKA = Fsl*IbarNaK*fnak*ions_k_o /(1+(KmNaip_PKA/ions_na_sl)^4) /(ions_k_o+KmKo);

I_NaK_junc = (1 - fINaK_PKA) * I_NaK_junc_noPKA + fINaK_PKA * I_NaK_junc_PKA;
I_NaK_sl = (1 - fINaK_PKA) * I_NaK_sl_noPKA + fINaK_PKA * I_NaK_sl_PKA;
I_NaK = I_NaK_junc+I_NaK_sl;

%% I_Kr (Rapidly Activating K Current, hERG current)
[I_Kr, d_ikr_c0, d_ikr_c1, d_ikr_c2, d_ikr_o, d_ikr_i ] = getIKr(membrane_v,ikr_c0,ikr_c1, ikr_c2, ikr_o, ikr_i,...
    ions_k_o, ek, vfrt, cellType, IKr_Multiplier);

%% I_Ks (Slowly Activating K Current)
pNaK = 0.01833;
eks = (1/FoRT)*log((ions_k_o+pNaK*ions_na_o)/(ions_k_i+pNaK*ions_na_i));
[I_Ks, d_iks_xs_junc, d_iks_xs_sl ] = getIKs_Bartos(membrane_v, iks_xs_junc, iks_xs_sl, ions_ca_junc, ions_ca_sl, eks, Fjunc, Fsl, fIKs_PKA, cellType, IKs_Multiplier);

%% IKb (Background potassium current)
xkb=1.0/(1.0+exp(-(membrane_v-10.8968)/(23.9871)));
GKb= 0.010879* IKb_Multiplier;
I_Kb=GKb*xkb*(membrane_v-ek);


%% I_to: Transient Outward K Current (slow and fast components)
fItop = camk_f_RyR;

[I_to, I_tof, I_tos, d_ito_xtos, d_ito_ytos, d_ito_xtof, d_ito_ytof, d_ito_xtos_p, d_ito_ytos_p, d_ito_xtof_p, d_ito_ytof_p] = getIto(membrane_v, ito_xtos, ito_ytos, ito_xtof, ito_ytof, ito_xtos_p, ito_ytos_p, ito_xtof_p, ito_ytof_p, ek, fItop, cellType, Itos_Multiplier, Itof_Multiplier);

%% I_K1: Time-Independent K Current
I_K1 = getIK1_CRLP(membrane_v, ions_k_o, ek, cellType, IK1_Multiplier);

%% I_ClCa, I_Clbk (Ca-activated Cl current, Background Cl current)
GClCa = IClCa_Multiplier * 0.01615;   % [mS/uF]
GClB = IClb_Multiplier * 0.00241;        
KdClCa = 100e-3;    

I_ClCa_junc = 0.5 * GClCa/(1+KdClCa/ions_ca_junc)*(membrane_v-ecl);
I_ClCa_sl = 0.5 * GClCa/(1+KdClCa/ions_ca_sl)*(membrane_v-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;
I_Clbk = GClB*(membrane_v-ecl);

%% I_Ca (L-type Calcium Current)
fICaLp = camk_f_ICaL;

if (fICaL_PKA > 0) % If PKA is active
    [I_Ca_junc,I_CaNa_junc,ICaK_junc,I_Ca_sl,I_CaNa_sl,ICaK_sl,d_ical_d, d_ical_ff,d_ical_fs,d_ical_fcaf,d_ical_fcas,d_ical_jca,d_ical_nca,d_ical_nca_i,...
        d_ical_ffp, d_ical_fcafp, d_ical_d_P,d_ical_ff_P,d_ical_fs_P,d_ical_fcaf_P,d_ical_fcas_P,d_ical_fBPf,d_ical_fcaBPf, d_ical_pureCDI_junc, d_ical_pureCDI_sl, PhiCaL_junc, PhiCaL_i, gammaCaoMyo, gammaCaiMyo] = getICaL(membrane_v, ical_d,ical_ff,ical_fs,ical_fcaf,ical_fcas,ical_jca,ical_nca,ical_nca_i,ical_ffp,ical_fcafp,...
        ical_d_P,ical_ff_P,ical_fs_P,ical_fcaf_P,ical_fcas_P,ical_fBPf,ical_fcaBPf, ...
        ical_pureCDI_junc, ical_pureCDI_sl, fICaLp, fICaL_PKA, ions_ca_i, ions_ca_sl, ions_ca_junc, ions_ca_o, ions_na_i, ions_na_sl, ions_na_o, ions_k_i, ions_k_o, ions_cl_i, ions_cl_o, cellType, ICaL_fraction_junc, ICaL_Multiplier, extraParams );
else
    [I_Ca_junc,I_CaNa_junc,ICaK_junc,I_Ca_sl,I_CaNa_sl,ICaK_sl,d_ical_d, d_ical_ff,d_ical_fs,d_ical_fcaf,d_ical_fcas,d_ical_jca,d_ical_nca,d_ical_nca_i,...
        d_ical_ffp, d_ical_fcafp, d_ical_d_P,d_ical_ff_P,d_ical_fs_P,d_ical_fcaf_P,d_ical_fcas_P,d_ical_fBPf,d_ical_fcaBPf, d_ical_pureCDI_junc, d_ical_pureCDI_sl, PhiCaL_junc, PhiCaL_i, gammaCaoMyo, gammaCaiMyo] = getICaL_noPKA(membrane_v, ical_d,ical_ff,ical_fs,ical_fcaf,ical_fcas,ical_jca,ical_nca,ical_nca_i,ical_ffp,ical_fcafp,...
        ical_d_P,ical_ff_P,ical_fs_P,ical_fcaf_P,ical_fcas_P,ical_fBPf,ical_fcaBPf, ...
        ical_pureCDI_junc, ical_pureCDI_sl, fICaLp, fICaL_PKA, ions_ca_i, ions_ca_sl, ions_ca_junc, ions_ca_o, ions_na_i, ions_na_sl, ions_na_o, ions_k_i, ions_k_o, ions_cl_i, ions_cl_o, cellType, ICaL_fraction_junc, ICaL_Multiplier, extraParams );
end


I_CaK = ICaK_junc + ICaK_sl;
I_Ca = I_Ca_junc+I_Ca_sl;
I_CaNa = I_CaNa_junc+I_CaNa_sl;
I_Catot = I_Ca+I_CaK+I_CaNa;

%% I_ncx (Na/Ca Exchanger flux)
[ I_NCX_sl, I_NCX_junc] = getINaCa_ORd2011(membrane_v,F,R,T, ions_na_junc, ions_na_sl, ions_na_i, ions_na_o, ions_ca_junc, ions_ca_sl, ions_ca_i, ions_ca_o, cellType, INaCa_Multiplier, INaCa_fraction_junc);

%% I_pca (Sarcolemmal Ca Pump Current)
IbarSLCaP = 0.02064  * IpCa_Multiplier;
KmPCa = 0.5e-3;     % [mM]

Q10SLCaP = 2.35;    % [none]
I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*ions_ca_junc^1.6/(KmPCa^1.6+ions_ca_junc^1.6);
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*ions_ca_sl^1.6/(KmPCa^1.6+ions_ca_sl^1.6);
I_pca = I_pca_junc+I_pca_sl;

%% I_Cabk (Background calcium current)
GCaB = ICab_Multiplier * 5.15575e-04;    % [uA/uF]
I_Cabk_junc = Fjunc*GCaB*(membrane_v-eca_junc);
I_Cabk_sl = Fsl*GCaB*(membrane_v-eca_sl);
I_Cabk = I_Cabk_junc+I_Cabk_sl;

%% SR release and leak
[J_SRCarel, J_SRleak, Jrel_ICaLdep, d_ryr_R, d_ryr_O, d_ryr_I, d_ryr_CaRI, d_ryr_R_p, d_ryr_O_p, d_ryr_I_p, d_ryr_CaRI_p, d_jrel_icaldep_act, d_jrel_icaldep_f1, d_jrel_icaldep_f2] = getJrel(ions_ca_junc, ions_ca_SR, I_Ca_junc, ryr_R, ryr_O, ryr_I, ryr_CaRI, ryr_R_p, ryr_O_p, ryr_I_p, jrel_icaldep_act, jrel_icaldep_f1, jrel_icaldep_f2, ryr_CaRI_p, camk_f_RyR, Jrel_Multiplier, extraParams);

%% SR reuptake via SERCA
J_serca = getJup(ions_ca_SR, ions_ca_i, Temp, camk_f_PLB, fPLB_PKA, casig_SERCA_act, cellType, Jup_Multiplier, extraParams);

%% Stimulation 
amp=stimAmp;
duration=stimDur;
if t<=duration
    Istim=amp;
else
    Istim=0.0;
end

%% Contraction
[Ta, d_contraction_XS, d_contraction_XW, d_contraction_Ca_TRPN, d_contraction_TmBlocked, d_contraction_ZETAS, d_contraction_ZETAW] = getContractionLand(ions_ca_i, contraction_XS, contraction_XW, contraction_Ca_TRPN, contraction_TmBlocked, contraction_ZETAS, contraction_ZETAW,  fTnI_PKA, fMyBPC_PKA, extraParams);

%% Sodium and Calcium Buffering
Bmax_Naj = 7.561;       % [mM] 
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      
% koff_tncl = 19.6e-3;    % [1/ms]
% kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      
koff_tnchca = 0.032e-3; % [1/ms]
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms]
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_CaM = 24e-3 * CMDN_Multiplier;       % [mM] 
koff_cam = 238e-3;      % [1/ms]
kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 17.85854e-3;     % [mM]  
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 33.923e-3*Vmyo/Vsl;      % [mM]    
Bmax_SLlowj = 4.89983e-4*Vmyo/Vjunc;    % [mM]   
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 12.15423e-3*Vmyo/Vsl;  % [mM]
Bmax_SLhighj = 1.75755e-4*Vmyo/Vjunc;  % [mM]
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 136.55214e-3*Vmyo/Vsr;    % [mM]   
koff_csqn = 65;         % [1/ms]
kon_csqn = 100;         % [1/mM/ms]
d_buffers_NaBj = kon_na*ions_na_junc*(Bmax_Naj-buffers_NaBj)-koff_na*buffers_NaBj;        % [mM/ms]
d_buffers_NaBsl = kon_na*ions_na_sl*(Bmax_Nasl-buffers_NaBsl)-koff_na*buffers_NaBsl;       % [mM/ms]

d_buffers_TnClow = Bmax_TnClow * d_contraction_Ca_TRPN; % myofilament - based on the contraction model
d_buffers_TnCHc = kon_tnchca*ions_ca_i*(Bmax_TnChigh-buffers_TnCHc-buffers_TnCHm)-koff_tnchca*buffers_TnCHc; %  [mM/ms]
d_buffers_TnCHm = kon_tnchmg*ions_mg_i*(Bmax_TnChigh-buffers_TnCHc-buffers_TnCHm)-koff_tnchmg*buffers_TnCHm;   %  [mM/ms]
d_buffers_CaM = kon_cam*ions_ca_i*(Bmax_CaM-buffers_CaM)-koff_cam*buffers_CaM;                 %  [mM/ms]
d_buffers_Myosin_ca = kon_myoca*ions_ca_i*(Bmax_myosin-buffers_Myosin_ca-buffers_Myosin_mg)-koff_myoca*buffers_Myosin_ca;    % [mM/ms]
d_buffers_Myosin_mg = kon_myomg*ions_mg_i*(Bmax_myosin-buffers_Myosin_ca-buffers_Myosin_mg)-koff_myomg*buffers_Myosin_mg;      % [mM/ms]
d_buffers_SRB = kon_sr*ions_ca_i*(Bmax_SR-buffers_SRB)-koff_sr*buffers_SRB;                    %  [mM/ms]

J_CaB_cytosol = d_buffers_TnClow + d_buffers_TnCHc + d_buffers_CaM + d_buffers_Myosin_ca + d_buffers_SRB;

% Junctional and SL Ca Buffers
d_buffers_SLLj = kon_sll*ions_ca_junc*(Bmax_SLlowj-buffers_SLLj)-koff_sll*buffers_SLLj;       %  [mM/ms]
d_buffers_SLLsl = kon_sll*ions_ca_sl*(Bmax_SLlowsl-buffers_SLLsl)-koff_sll*buffers_SLLsl;      %  [mM/ms]
d_buffers_SLHj = kon_slh*ions_ca_junc*(Bmax_SLhighj-buffers_SLHj)-koff_slh*buffers_SLHj;      %  [mM/ms]
d_buffers_SLHsl = kon_slh*ions_ca_sl*(Bmax_SLhighsl-buffers_SLHsl)-koff_slh*buffers_SLHsl;     %  [mM/ms]
J_CaB_junction = d_buffers_SLLj + d_buffers_SLHj;
J_CaB_sl = d_buffers_SLLsl + d_buffers_SLHsl;

%% Ion concentrations
% SR Ca Concentrations
d_buffers_Csqn = kon_csqn*ions_ca_SR*(Bmax_Csqn-buffers_Csqn)-koff_csqn*buffers_Csqn;       % Csqn      [mM/ms]
d_ions_ca_SR = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel) - d_buffers_Csqn;         % Ca_sr     [mM/ms] %Ratio 3 leak current

% Sodium Concentrations
I_Na_tot_junc = I_Na_junc+I_Nabk_junc+3*I_NCX_junc+3*I_NaK_junc+I_CaNa_junc;   % [uA/uF]
I_Na_tot_sl = I_Na_sl+I_Nabk_sl+3*I_NCX_sl+3*I_NaK_sl+I_CaNa_sl;   % [uA/uF]

d_ions_na_junc = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(ions_na_sl-ions_na_junc)-d_buffers_NaBj;
d_ions_na_sl = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(ions_na_junc-ions_na_sl)...
    +J_na_slmyo/Vsl*(ions_na_i-ions_na_sl)-d_buffers_NaBsl;
d_ions_na_i = J_na_slmyo/Vmyo*(ions_na_sl-ions_na_i);             % [mM/msec]

% Potassium Concentration
I_K_tot = I_to+I_Kr+I_Ks+I_K1-2*I_NaK+I_CaK+I_Kb+Istim;     % [uA/uF]
d_ions_k_i = -I_K_tot*Cmem/(Vmyo*Frdy);           % [mM/msec]

% Cl concentrations
I_Cl_tot = I_ClCa+I_Clbk;

d_ions_cl_i = -I_Cl_tot*Cmem/(-1*Vmyo*Frdy);           % [mM/msec]

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_Cabk_junc+I_pca_junc-2*I_NCX_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_Cabk_sl+I_pca_sl-2*I_NCX_sl;            % [uA/uF]
d_ions_ca_junc = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(ions_ca_sl-ions_ca_junc)...
    -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  % Ca_j
d_ions_ca_sl = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(ions_ca_junc-ions_ca_sl)...
    + J_ca_slmyo/Vsl*(ions_ca_i-ions_ca_sl)-J_CaB_sl;   % Ca_sl
d_ions_ca_i = -J_serca*Vsr/Vmyo-J_CaB_cytosol+J_ca_slmyo/Vmyo*(ions_ca_sl-ions_ca_i);


%% Membrane Potential
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk;                        % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
d_membrane_v = -(I_tot);

%% Derivative
ydot = [    d_membrane_v;    d_ions_na_junc;    d_ions_na_sl;    d_ions_na_i;    d_ions_k_i;    d_ions_ca_junc;    d_ions_ca_sl;    d_ions_ca_i;    d_ions_cl_i;    d_ions_ca_SR;  ...
    d_buffers_NaBj;    d_buffers_NaBsl;    d_buffers_TnClow; d_buffers_TnCHc;   d_buffers_TnCHm;    d_buffers_CaM;    d_buffers_Myosin_ca;    d_buffers_Myosin_mg;    d_buffers_SRB;    d_buffers_SLLj;
    d_buffers_SLLsl; d_buffers_SLHj;    d_buffers_SLHsl;    d_buffers_Csqn;    d_ina_m;    d_ina_h;    d_ina_j;    d_inal_ml;    d_inal_hl;    d_inal_hl_p;...
    d_ito_xtos; d_ito_ytos;    d_ito_xtof;    d_ito_ytof;    d_ito_xtos_p;    d_ito_ytos_p;    d_ito_xtof_p;    d_ito_ytof_p;    d_iks_xs_junc;    d_iks_xs_sl;  ...
    d_ikr_c0; d_ikr_c1;    d_ikr_c2;    d_ikr_o;    d_ikr_i;    d_ical_d;    d_ical_ff;    d_ical_fs;    d_ical_fcaf;    d_ical_fcas;   ...
    d_ical_jca; d_ical_nca;    d_ical_nca_i;    d_ical_ffp;    d_ical_fcafp;    d_ical_pureCDI_junc;    d_ical_pureCDI_sl;    d_ryr_R;    d_ryr_O;    d_ryr_I;    ...
    d_ryr_CaRI; d_ryr_R_p;    d_ryr_O_p;    d_ryr_I_p;    d_ryr_CaRI_p;    d_jrel_icaldep_act;    d_jrel_icaldep_f1;    d_jrel_icaldep_f2;    d_contraction_XS;    d_contraction_XW;    ...
    d_contraction_Ca_TRPN;d_contraction_TmBlocked;    d_contraction_ZETAS;    d_contraction_ZETAW;    d_camk_trap;    d_camk_f_ICaL;    d_camk_f_RyR;    d_camk_f_PLB;    d_casig_serca_trap; d_ina_hp; ...
    d_ina_jp; d_ina_m_PKA;  d_ina_h_PKA; d_ina_j_PKA; d_ina_h_dualP; d_ina_j_dualP; d_ical_d_P; d_ical_ff_P; d_ical_fs_P; d_ical_fcaf_P;...
    d_ical_fcas_P; d_ical_fBPf; d_ical_fcaBPf];

if (runSignallingPathway)
    ydot = [ydot; ydot_signalling];
end

if flag_ode==1
    output = ydot;
else
    output = [I_NaFast, I_NaL, I_to, I_Ca, I_Kr, I_Ks, I_K1, I_NCX_junc, I_NCX_sl, I_NaK, I_Kb, I_Nabk, I_Cabk, I_pca, J_SRCarel,...
        J_SRleak, I_ClCa, I_ClCa_junc, I_ClCa_sl, I_Clbk, I_Catot, J_serca, I_Ca_junc, I_Ca_sl,Jrel_ICaLdep, I_tof, I_tos, Ta, ...
        fINa_PKA, fICaL_PKA, fINaK_PKA, fIKs_PKA, fPLB_PKA, fTnI_PKA, fMyBPC_PKA];
end

end


%% CaMKII and other Ca signalling
function [d_camk_trap, d_camk_f_ICaL, d_camk_f_RyR, d_camk_f_PLB, d_casig_serca_trap, casig_SERCA_act] = getCaSignalling(ions_ca_junc, camk_trap, camk_f_ICaL, camk_f_RyR, camk_f_PLB, casig_serca_trap, Whole_cell_PP1)
PP1_tot = Whole_cell_PP1;%0.13698; %

CaMK0  = 2*0.05; % Equilibrium fraction of active CaMKII binding sites
Km_CaMK_Ca = 5*0.0015; %[mmol/L] CaMKII affinity for Ca2+/CaM activation %Adjusted because of smaller cleft space

bound = CaMK0 * (1 - camk_trap) / (1 + Km_CaMK_Ca / ions_ca_junc);
CaMK_active = bound + camk_trap; % Fraction of active CaMKII

alpha = 0.05;
beta  = 6.8e-4;
d_camk_trap = alpha * bound * CaMK_active - beta * camk_trap * (0.1 + 0.9 * PP1_tot / 0.1371); %dCaMK_Trap/dt

tau_plb = 100000; % [ms] Time constant of CaMKII PLB phosphorylation
tau_ryr = 10000; % [ms] Time constant of CaMKII RyR phosphorylation
tau_cal = tau_ryr; % Time constant of CaMKII ICaL phosphorylation


K_Phos_CaMK = 0.35;  % Affinity of PLB, ICaL etc for CaMKII
CaMK_Phos_junc_ICaL =  CaMK_active / (CaMK_active + K_Phos_CaMK);
CaMK_Phos_junc_RyR =  CaMK_active / (CaMK_active + 1);
CaMK_Phos_junc_PLB =  CaMK_active / (CaMK_active + 10);

d_camk_f_ICaL = (CaMK_Phos_junc_ICaL - camk_f_ICaL) / tau_cal;
d_camk_f_RyR = (CaMK_Phos_junc_RyR - camk_f_RyR)  / tau_ryr;
d_camk_f_PLB = (CaMK_Phos_junc_PLB - camk_f_PLB)  / tau_plb;

alpha_serca = 0.05;
bound_serca = CaMK0 * (1 - casig_serca_trap) / (1 + Km_CaMK_Ca / ions_ca_junc);
casig_SERCA_act = bound_serca + casig_serca_trap; % Fraction of active CaMKII
d_casig_serca_trap = alpha_serca * bound_serca * casig_SERCA_act - beta * casig_serca_trap * (0.1 + 0.9 * PP1_tot / 0.1371); %dCaMK_Trap/dt
end

%% INa
function [I_NaFast_junc, I_NaFast_sl, dm, dh, dhp, dj, djp, dm_P, dh_P, dj_P, dhp_P, djp_P] = getINa_Grandi(v, m, h, hp, j, jp, fINap, ena_junc, ena_sl, INa_Multiplier,fINaP, Fjunc, m_P, h_P,j_P,hp_P,jp_P)
% m gate
mss = 1 / ((1 + exp( -(56.86 + v) / 9.03 ))^2);
taum = 0.1292 * exp(-((v+45.79)/15.54)^2) + 0.06487 * exp(-((v-4.823)/51.12)^2);
dm = (mss - m) / taum;

% h gate
ah = (v >= -40) * (0)...
    + (v < -40) * (0.057 * exp( -(v + 80) / 6.8 ));
bh = (v >= -40) * (0.77 / (0.13*(1 + exp( -(v + 10.66) / 11.1 )))) ...
    + (v < -40) * ((2.7 * exp( 0.079 * v) + 3.1*10^5 * exp(0.3485 * v)));
tauh = 1 / (ah + bh);
hss = 1 / ((1 + exp( (v + 71.55)/7.43 ))^2);

dh = (hss - h) / tauh;

% j gate
aj = (v >= -40) * (0) ...
    +(v < -40) * (((-2.5428 * 10^4*exp(0.2444*v) - 6.948*10^-6 * exp(-0.04391*v)) * (v + 37.78)) / ...
    (1 + exp( 0.311 * (v + 79.23) )));
bj = (v >= -40) * ((0.6 * exp( 0.057 * v)) / (1 + exp( -0.1 * (v + 32) ))) ...
    + (v < -40) * ((0.02424 * exp( -0.01052 * v )) / (1 + exp( -0.1378 * (v + 40.14) )));
tauj = 1 / (aj + bj);
jss = 1 / ((1 + exp( (v + 71.55)/7.43 ))^2);
dj = (jss - j) / tauj;

% gating CaMK-P
hssp = 1 / ((1 + exp( (v + 71.55 + 6 )/7.43 ))^2);
dhp = (hssp - hp) / tauh;
taujp = 1.46 * tauj;
djp = (jss - jp) / taujp;

% gating PKA
mss_P =  1 / ((1 + exp( -(56.86 + v + 5) / 9.03 ))^2);
taum = 0.1292 * exp(-((v+45.79)/15.54)^2) + 0.06487 * exp(-((v-4.823)/51.12)^2);
dm_P = (mss_P - m_P) / taum;

hss_P = 1 / ((1 + exp( (v + 71.55+5.0)/7.43 ))^2);  %BetaAdrenergic;
dh_P=(hss_P-h_P)/tauh;   %BetaAdrenergic;
jss_P=hss_P;
dj_P = (jss_P - j_P) / tauj;

% Both phosphorylated
hssp_P = 1 / ((1 + exp( (v + 71.55 + 6 + 5.0)/7.43 ))^2); %BetaAdrenergic;
dhp_P = (hssp_P - hp_P) / tauh;
jssp_P=hssp_P;
djp_P = (jssp_P - jp_P) / taujp;

% Putting together the channels behavior and fraction
GNa = 22.08788 * INa_Multiplier;
GNa_P= GNa * 1.25; 

fINa_P = fINaP; 
fINa_BP = fINap*fINa_P ;
fINa_CaMKonly = fINap-fINa_BP ;
fINa_PKAonly = fINa_P-fINa_BP ;

INaBase_NP =  GNa*m^3.0*h*j ; % Non-Phosphorylated
INaBase_CaMK = GNa*m^3.0*hp*jp ;
INaBase_PKA = GNa_P*m_P^3.0*h_P*j_P ;
INaBase_BP = GNa_P*m_P^3.0*hp_P*jp_P ;

% Combining populations
INaBase =  ((1-fINa_CaMKonly-fINa_PKAonly-fINa_BP)*INaBase_NP + fINa_CaMKonly*INaBase_CaMK + fINa_PKAonly*INaBase_PKA + fINa_BP*INaBase_BP) ;
I_NaFast_junc = Fjunc * INaBase*(v-ena_junc);
I_NaFast_sl = (1 - Fjunc) * INaBase * (v-ena_sl);
end


%% INaL
function [I_NaL_junc, I_NaL_sl, dmL, dhL, dhLp] = getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ena_junc, ena_sl, celltype, Fjunc, Fsl, INaL_Multiplier)
mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
tm = 0.1292 * exp(-((v+45.79)/15.54)^2) + 0.06487 * exp(-((v-4.823)/51.12)^2);
tmL=tm;
dmL=(mLss-mL)/tmL;
hLss=1.0/(1.0+exp((v+87.61)/7.488));
thL=145;
dhL=(hLss-hL)/thL;
hLssp=1.0/(1.0+exp((v+93.81)/7.488));
thLp=3.0*thL;
dhLp=(hLssp-hLp)/thLp;
GNaL=0.04229 * INaL_Multiplier * (1+fINaLp); 
if celltype==1
    GNaL=GNaL*0.7;
end

I_NaL_junc = Fjunc * GNaL*(v-ena_junc)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);
I_NaL_sl = Fsl * GNaL*(v-ena_sl)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);
end



function [ICaL_junc,ICaNa_junc,ICaK_junc,ICaL_sl,ICaNa_sl,ICaK_sl,dd,dff,dfs,dfcaf,dfcas,...
    djca,dnca, dnca_i, dffp,dfcafp, dd_P,dff_P,dfs_P,dfcaf_P,dfcas_P,dfBPf,dfcaBPf, d_ical_pureCDI_junc, d_ical_pureCDI_sl, PhiCaL_junc, PhiCaL_i, gammaCaoMyo, gammaCaiMyo] = getICaL(v, d,ff,fs,fcaf,fcas,jca,nca, nca_i,ffp,fcafp,...
    d_P,ff_P,fs_P,fcaf_P,fcas_P,fBPf,fcaBPf,...
    ical_pureCDI_junc, ical_pureCDI_sl, fICaLp, fICaLP, ions_ca_i, ions_ca_sl, ions_ca_junc, ions_ca_o, ions_na_i, ions_na_sl, ions_na_o, ions_k_i, ko, ions_cl_i, clo, celltype, ICaL_fractionSS, ICaL_PCaMultiplier, extraParams)


%physical constants
R=8314.0;
T=310.0;
F=96485.0;
vffrt=v*F*F/(R*T);
vfrt=v*F/(R*T);

%calculate ICaL, ICaNa, ICaK
dss=min(1.0763*exp(-1.0070*exp(-0.0829*(v+3.62483))), 1);
td= 1.5+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
dd=(dss-d)/td;
fss=1.0/(1.0+exp((v+19.58)/3.696));

tff=6.17111+1.0/(0.00126*exp(-(v+26.63596)/(9.69961))+0.00126*exp((v+26.63596)/(9.69961)));
tfs=2719.22489+1.0/(7.19411e-05*exp(-(v+5.74631)/(10.87690))+7.19411e-05*exp((v+5.74631)/(16.31535)));

Aff=0.52477;
Afs=1.0-Aff;
dff=(fss-ff)/tff;
dfs=(fss-fs)/tfs;
f=Aff*ff+Afs*fs;
fcass=fss;
tfcaf=13.50673+1.0/(0.15420*exp(-(v-1.31611)/(11.33960))+0.15420*exp((v-1.31611)/(11.33960)));
tfcas=177.95813+1.0/(4.73955e-04*exp((-v+ 0.79049)/(0.81777)) + 4.73955e-04*exp((v+2.40474)/(1.90812)));

Afcaf=0.3+0.6/(1.0+exp((v-9.24247)/(27.96201)));

Afcas=1.0-Afcaf;
dfcaf=(fcass-fcaf)/tfcaf;
dfcas=(fcass-fcas)/tfcas;
fca=Afcaf*fcaf+Afcas*fcas;

tjca = 66;
jcass = 1.0/(1.0+exp((v+17.66945)/(3.21501)));
djca=(jcass-jca)/tjca;
tffp=2.5*tff;
dffp=(fss-ffp)/tffp;
fp=Aff*ffp+Afs*fs;
tfcafp=2.5*tfcaf;
dfcafp=(fcass-fcafp)/tfcafp;
fcap=Afcaf*fcafp+Afcas*fcas;

%% SS nca
Kmn=0.00222;
k2n=957.85903;
km2n=jca*0.84191;
anca=1.0/(k2n/km2n+(1.0+Kmn/ions_ca_junc)^(3.80763));
dnca=anca*k2n-nca*km2n;

%% myoplasmic nca
anca_i = 1.0/(k2n/km2n+(1.0+Kmn/ions_ca_sl)^(3.80763));
dnca_i = anca_i*k2n-nca_i*km2n;

%% junctional driving force
Io = 0.5*(ions_na_o + ko + clo + 4*ions_ca_o)/1000 ; % ionic strength outside. /1000 is for things being in micromolar
Ii = 0.5*(ions_na_sl + ions_k_i + ions_cl_i + 4*ions_ca_sl)/1000 ; % ionic strength inside. /1000 is for things being in micromolar % USING NONLOCAL CA
% The ionic strength is too high for basic DebHuc. We'll use Davies
dielConstant = 74; % water at 37°.
temp = 310; % body temp in kelvins.
constA = 1.82*10^6*(dielConstant*temp)^(-1.5);

gamma_cai = 10^(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_cao = 10^(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_nai = 10^(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_nao = 10^(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_ki = 10^(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_ko = 10^(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));

PhiCaL_junc =  4.0*vffrt*(gamma_cai*ions_ca_sl*exp(2.0*vfrt)-gamma_cao*ions_ca_o)/(exp(2.0*vfrt)-1.0); % USING NONLOCAL CA
PhiCaNa_junc =  1.0*vffrt*(gamma_nai*ions_na_sl*exp(1.0*vfrt)-gamma_nao*ions_na_o)/(exp(1.0*vfrt)-1.0);
PhiCaK_junc =  1.0*vffrt*(gamma_ki*ions_k_i*exp(1.0*vfrt)-gamma_ko*ko)/(exp(1.0*vfrt)-1.0);

%% sl/myo driving force
Io = 0.5*(ions_na_o + ko + clo + 4*ions_ca_o)/1000 ; % ionic strength outside. /1000 is for things being in micromolar
Ii = 0.5*(ions_na_i + ions_k_i + ions_cl_i + 4*ions_ca_i)/1000 ; % ionic strength inside. /1000 is for things being in micromolar % USING NONLOCAL CA
% The ionic strength is too high for basic DebHuc. We'll use Davies
dielConstant = 74; % water at 37°.
temp = 310; % body temp in kelvins.
constA = 1.82*10^6*(dielConstant*temp)^(-1.5);

gamma_cai = 10^(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_cao = 10^(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_nai = 10^(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_nao = 10^(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_ki = 10^(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_ko = 10^(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));

gammaCaoMyo = gamma_cao;
gammaCaiMyo = gamma_cai;

PhiCaL_i =  4.0*vffrt*(gamma_cai*ions_ca_i*exp(2.0*vfrt)-gamma_cao*ions_ca_o)/(exp(2.0*vfrt)-1.0); % USING NONLOCAL CA
PhiCaNa_i =  1.0*vffrt*(gamma_nai*ions_na_i*exp(1.0*vfrt)-gamma_nao*ions_na_o)/(exp(1.0*vfrt)-1.0);
PhiCaK_i =  1.0*vffrt*(gamma_ki*ions_k_i*exp(1.0*vfrt)-gamma_ko*ko)/(exp(1.0*vfrt)-1.0);

%% Calculating "pure" CDI
sigmoidTransition = 1 - 1/(1 + (1.86532 * ions_ca_junc/0.032));
tauTransition = 1.09670 + (1-sigmoidTransition)*141.42990;

rateRecovery = 0.02313;

d_ical_pureCDI_junc = -ical_pureCDI_junc*sigmoidTransition/tauTransition + (1-ical_pureCDI_junc) * rateRecovery;

% and for nondyadic
sigmoidTransition2 = 1 - 1/(1 + (1.86532 * ions_ca_sl/0.032));
tauTransition2 = 1.09670 + (1-sigmoidTransition2)*141.42990;

d_ical_pureCDI_sl = -ical_pureCDI_sl*sigmoidTransition2/tauTransition2 + (1-ical_pureCDI_sl) * rateRecovery;
%% The rest
PCa=1.5768e-04 * ICaL_PCaMultiplier;

if celltype==1
    PCa=PCa*1.025;
elseif celltype==2
    PCa=PCa*1.1;
end

PCap=1.1*PCa;
PCaNa=1.1737/1.8969*0.00125*PCa;
PCaK=1.1737/1.8969*3.574e-4*PCa;
PCaNap=1.1737/1.8969*0.00125*PCap;
PCaKp=1.1737/1.8969*3.574e-4*PCap;


%% PKA phosphorylation and other ones
dPss=1.0323*exp(-1.0553*exp(-0.0810*(v+12.62483)));   
dPss=min(dPss,1);
dd_P=(dPss-d_P)/td;
fss_P=1.0/(1.0+exp((v+19.58+6)/3.696)); 
dff_P=(fss_P-ff_P)/tff;
dfs_P=(fss_P-fs_P)/tfs;
fcass_P=fss_P;
dfcaf_P=(fcass_P-fcaf_P)/tfcaf;
dfcas_P=(fcass_P-fcas_P)/tfcas;
f_P=Aff*ff_P+Afs*fs_P;
fcap_P=Afcaf*fcaf_P+Afcas*fcas_P;

%% The rest (PKA)
PCa_P = PCa * 1.9;  
if celltype==1
    PCa_P=PCa_P*1.025;
elseif celltype==2
    PCa_P=PCa_P*1.1;
end

PCaNa_P=1.1737/1.8969*0.00125*PCa_P;
PCaK_P=1.1737/1.8969*3.574e-4*PCa_P;

% Both-P population takes dP, for f gate, takes ss from PKA and tau from CaMK (only fast component was modified)
fBPss = fss_P;
dfBPf = (fBPss-fBPf)/tffp;
fBP = Aff*fBPf+Afs*fs_P; % only fast component modified, slow component same as fPs
fcaBPss = fcass_P;
dfcaBPf = (fcaBPss-fcaBPf)/tfcafp;
fcaBP = Afcaf*fcaBPf+Afcas*fcas_P;

fICaL_P = fICaLP; % PKA-P fraction
fICaL_BP = fICaLp*fICaL_P;
fICaL_CaMKonly = fICaLp-fICaL_BP;
fICaL_PKAonly = fICaL_P-fICaL_BP;

ICaL_junc_NP=PCa*PhiCaL_junc*d*(f*(1.0-nca)+jca*fca*nca);
ICaL_junc_CaMK=PCap*PhiCaL_junc*d*(fp*(1.0-nca)+jca*fcap*nca);
ICaL_junc_PKA=PCa_P*PhiCaL_junc*d_P*(f_P*(1.0-nca)+jca*fcap_P*nca);
ICaL_junc_BP=PCa_P*PhiCaL_junc*d_P*(fBP*(1.0-nca)+jca*fcaBP*nca);
ICaL_i_NP=PCa*PhiCaL_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
ICaL_i_CaMK=PCap*PhiCaL_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
ICaL_i_PKA=PCa_P*PhiCaL_i*d_P*(f_P*(1.0-nca_i)+jca*fcap_P*nca_i);
ICaL_i_BP=PCa_P*PhiCaL_i*d_P*(fBP*(1.0-nca_i)+jca*fcaBP*nca_i);

ICaNa_junc_NP=PCaNa*PhiCaNa_junc*d*(f*(1.0-nca)+jca*fca*nca);
ICaNa_junc_CaMK=PCaNap*PhiCaNa_junc*d*(fp*(1.0-nca)+jca*fcap*nca);
ICaNa_junc_PKA=PCaNa_P*PhiCaNa_junc*d_P*(f_P*(1.0-nca)+jca*fcap_P*nca);
ICaNa_junc_BP=PCaNa_P*PhiCaNa_junc*d_P*(fBP*(1.0-nca)+jca*fcaBP*nca);
ICaNa_i_NP=PCaNa*PhiCaNa_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
ICaNa_i_CaMK=PCaNap*PhiCaNa_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
ICaNa_i_PKA=PCaNa_P*PhiCaNa_i*d_P*(f_P*(1.0-nca_i)+jca*fcap_P*nca_i);
ICaNa_i_BP=PCaNa_P*PhiCaNa_i*d_P*(fBP*(1.0-nca_i)+jca*fcaBP*nca_i);

ICaK_junc_NP=PCaK*PhiCaK_junc*d*(f*(1.0-nca)+jca*fca*nca);
ICaK_junc_CaMK=PCaKp*PhiCaK_junc*d*(fp*(1.0-nca)+jca*fcap*nca);
ICaK_junc_PKA=PCaK_P*PhiCaK_junc*d_P*(f_P*(1.0-nca)+jca*fcap_P*nca);
ICaK_junc_BP=PCaK_P*PhiCaK_junc*d_P*(fBP*(1.0-nca)+jca*fcaBP*nca);
ICaK_i_NP=PCaK*PhiCaK_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
ICaK_i_CaMK=PCaKp*PhiCaK_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);
ICaK_i_PKA=PCaK_P*PhiCaK_i*d_P*(f_P*(1.0-nca_i)+jca*fcap_P*nca_i);
ICaK_i_BP=PCaK_P*PhiCaK_i*d_P*(fBP*(1.0-nca_i)+jca*fcaBP*nca_i);


% 4 population combination and adding pure-CDI modulation
ICaL_junc = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaL_junc_NP + fICaL_CaMKonly*ICaL_junc_CaMK + fICaL_PKAonly*ICaL_junc_PKA + fICaL_BP*ICaL_junc_BP) * ical_pureCDI_junc;
ICaNa_junc = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaNa_junc_NP + fICaL_CaMKonly*ICaNa_junc_CaMK + fICaL_PKAonly*ICaNa_junc_PKA + fICaL_BP*ICaNa_junc_BP) * ical_pureCDI_junc;
ICaK_junc = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaK_junc_NP + fICaL_CaMKonly*ICaK_junc_CaMK + fICaL_PKAonly*ICaK_junc_PKA + fICaL_BP*ICaK_junc_BP) * ical_pureCDI_junc;
%
ICaL_sl = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaL_i_NP + fICaL_CaMKonly*ICaL_i_CaMK + fICaL_PKAonly*ICaL_i_PKA + fICaL_BP*ICaL_i_BP) * ical_pureCDI_sl;
ICaNa_sl = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaNa_i_NP + fICaL_CaMKonly*ICaNa_i_CaMK + fICaL_PKAonly*ICaNa_i_PKA + fICaL_BP*ICaNa_i_BP) * ical_pureCDI_sl;
ICaK_sl = ((1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaK_i_NP + fICaL_CaMKonly*ICaK_i_CaMK + fICaL_PKAonly*ICaK_i_PKA + fICaL_BP*ICaK_i_BP) * ical_pureCDI_sl;


% And we weight ICaL (in junc) and ICaL_i
ICaL_sl = ICaL_sl * (1-ICaL_fractionSS);
ICaNa_sl = ICaNa_sl * (1-ICaL_fractionSS);
ICaK_sl = ICaK_sl * (1-ICaL_fractionSS);
ICaL_junc = ICaL_junc * ICaL_fractionSS;
ICaNa_junc = ICaNa_junc * ICaL_fractionSS;
ICaK_junc = ICaK_junc * ICaL_fractionSS;
end


% A version which does not compute any PKA-dependent branches and is used only when PKA
% phosphorylation is zero. This is to save computing time - it makes a difference!
function [ICaL_junc,ICaNa_junc,ICaK_junc,ICaL_sl,ICaNa_sl,ICaK_sl,dd,dff,dfs,dfcaf,dfcas,...
    djca,dnca, dnca_i, dffp,dfcafp, dd_P,dff_P,dfs_P,dfcaf_P,dfcas_P,dfBPf,dfcaBPf, d_ical_pureCDI_junc, d_ical_pureCDI_sl, PhiCaL_junc, PhiCaL_i, gammaCaoMyo, gammaCaiMyo] = getICaL_noPKA(v, d,ff,fs,fcaf,fcas,jca,nca, nca_i,ffp,fcafp,...
    d_P,ff_P,fs_P,fcaf_P,fcas_P,fBPf,fcaBPf,...
    ical_pureCDI_junc, ical_pureCDI_sl, fICaLp, fICaLP, ions_ca_i, ions_ca_sl, ions_ca_junc, ions_ca_o, ions_na_i, ions_na_sl, ions_na_o, ions_k_i, ko, ions_cl_i, clo, celltype, ICaL_fractionSS, ICaL_PCaMultiplier, extraParams)

%physical constants
R=8314.0;
T=310.0;
F=96485.0;
vffrt=v*F*F/(R*T);
vfrt=v*F/(R*T);

%calculate ICaL, ICaNa, ICaK
dss=min(1.0763*exp(-1.0070*exp(-0.0829*(v+3.62483))), 1);
td= 1.5+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
dd=(dss-d)/td;
fss=1.0/(1.0+exp((v+19.58)/3.696));

tff=6.17111+1.0/(0.00126*exp(-(v+26.63596)/(9.69961))+0.00126*exp((v+26.63596)/(9.69961)));
tfs=2719.22489+1.0/(7.19411e-05*exp(-(v+5.74631)/(10.87690))+7.19411e-05*exp((v+5.74631)/(16.31535)));

Aff=0.52477;
Afs=1.0-Aff;
dff=(fss-ff)/tff;
dfs=(fss-fs)/tfs;
f=Aff*ff+Afs*fs;
fcass=fss;
tfcaf=13.50673+1.0/(0.15420*exp(-(v-1.31611)/(11.33960))+0.15420*exp((v-1.31611)/(11.33960)));
tfcas=177.95813+1.0/(4.73955e-04*exp((-v+ 0.79049)/(0.81777)) + 4.73955e-04*exp((v+2.40474)/(1.90812)));

Afcaf=0.3+0.6/(1.0+exp((v-9.24247)/(27.96201)));

Afcas=1.0-Afcaf;
dfcaf=(fcass-fcaf)/tfcaf;
dfcas=(fcass-fcas)/tfcas;
fca=Afcaf*fcaf+Afcas*fcas;

tjca = 66;
jcass = 1.0/(1.0+exp((v+17.66945)/(3.21501)));
djca=(jcass-jca)/tjca;
tffp=2.5*tff;
dffp=(fss-ffp)/tffp;
fp=Aff*ffp+Afs*fs;
tfcafp=2.5*tfcaf;
dfcafp=(fcass-fcafp)/tfcafp;
fcap=Afcaf*fcafp+Afcas*fcas;

%% SS nca
Kmn=0.00222;
k2n=957.85903;
km2n=jca*0.84191;
anca=1.0/(k2n/km2n+(1.0+Kmn/ions_ca_junc)^(3.80763));
dnca=anca*k2n-nca*km2n;

%% myoplasmic nca
anca_i = 1.0/(k2n/km2n+(1.0+Kmn/ions_ca_sl)^(3.80763));
dnca_i = anca_i*k2n-nca_i*km2n;

%% junctional driving force (nonlocal)
Io = 0.5*(ions_na_o + ko + clo + 4*ions_ca_o)/1000 ; % ionic strength outside. /1000 is for things being in micromolar
Ii = 0.5*(ions_na_sl + ions_k_i + ions_cl_i + 4*ions_ca_sl)/1000 ; % ionic strength inside. /1000 is for things being in micromolar % USING NONLOCAL CA
% The ionic strength is too high for basic DebHuc. We'll use Davies
dielConstant = 74; % water at 37°.
temp = 310; % body temp in kelvins.
constA = 1.82*10^6*(dielConstant*temp)^(-1.5);

gamma_cai = 10^(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_cao = 10^(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_nai = 10^(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_nao = 10^(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_ki = 10^(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_ko = 10^(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));


PhiCaL_junc =  4.0*vffrt*(gamma_cai*ions_ca_sl*exp(2.0*vfrt)-gamma_cao*ions_ca_o)/(exp(2.0*vfrt)-1.0); % USING NONLOCAL CA
PhiCaNa_junc =  1.0*vffrt*(gamma_nai*ions_na_sl*exp(1.0*vfrt)-gamma_nao*ions_na_o)/(exp(1.0*vfrt)-1.0);
PhiCaK_junc =  1.0*vffrt*(gamma_ki*ions_k_i*exp(1.0*vfrt)-gamma_ko*ko)/(exp(1.0*vfrt)-1.0);

%% sl/myo driving force
Io = 0.5*(ions_na_o + ko + clo + 4*ions_ca_o)/1000 ; % ionic strength outside. /1000 is for things being in micromolar
Ii = 0.5*(ions_na_i + ions_k_i + ions_cl_i + 4*ions_ca_i)/1000 ; % ionic strength inside. /1000 is for things being in micromolar % USING NONLOCAL CA
% The ionic strength is too high for basic DebHuc. We'll use Davies
dielConstant = 74; % water at 37°.
temp = 310; % body temp in kelvins.
constA = 1.82*10^6*(dielConstant*temp)^(-1.5);

gamma_cai = 10^(-constA * 4 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_cao = 10^(-constA * 4 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_nai = 10^(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_nao = 10^(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));
gamma_ki = 10^(-constA * 1 * (sqrt(Ii)/(1+sqrt(Ii))-0.3*Ii));
gamma_ko = 10^(-constA * 1 * (sqrt(Io)/(1+sqrt(Io))-0.3*Io));

gammaCaoMyo = gamma_cao;
gammaCaiMyo = gamma_cai;

PhiCaL_i =  4.0*vffrt*(gamma_cai*ions_ca_i*exp(2.0*vfrt)-gamma_cao*ions_ca_o)/(exp(2.0*vfrt)-1.0); % USING NONLOCAL CA
PhiCaNa_i =  1.0*vffrt*(gamma_nai*ions_na_i*exp(1.0*vfrt)-gamma_nao*ions_na_o)/(exp(1.0*vfrt)-1.0);
PhiCaK_i =  1.0*vffrt*(gamma_ki*ions_k_i*exp(1.0*vfrt)-gamma_ko*ko)/(exp(1.0*vfrt)-1.0);

%% Calculating "pure" CDI
sigmoidTransition = 1 - 1/(1 + (1.86532 * ions_ca_junc/0.032));
tauTransition = 1.09670 + (1-sigmoidTransition)*141.42990;

rateRecovery = 0.02313;

d_ical_pureCDI_junc = -ical_pureCDI_junc*sigmoidTransition/tauTransition + (1-ical_pureCDI_junc) * rateRecovery;

% and for nondyadic
sigmoidTransition2 = 1 - 1/(1 + (1.86532 * ions_ca_sl/0.032));
tauTransition2 = 1.09670 + (1-sigmoidTransition2)*141.42990;

d_ical_pureCDI_sl = -ical_pureCDI_sl*sigmoidTransition2/tauTransition2 + (1-ical_pureCDI_sl) * rateRecovery;
%% The rest
PCa=1.5768e-04 * ICaL_PCaMultiplier;
%
if celltype==1
    PCa=PCa*1.025;
elseif celltype==2
    PCa=PCa*1.1;
end

PCap=1.1*PCa;
PCaNa=1.1737/1.8969*0.00125*PCa;
PCaK=1.1737/1.8969*3.574e-4*PCa;
PCaNap=1.1737/1.8969*0.00125*PCap;
PCaKp=1.1737/1.8969*3.574e-4*PCap;

%
dd_P=0;
dff_P=0;
dfs_P=0;
dfcaf_P=0;
dfcas_P=0;
dfBPf = 0;
dfcaBPf = 0;
fICaL_CaMKonly = fICaLp;

ICaL_junc_NP=PCa*PhiCaL_junc*d*(f*(1.0-nca)+jca*fca*nca);
ICaL_junc_CaMK=PCap*PhiCaL_junc*d*(fp*(1.0-nca)+jca*fcap*nca);

ICaL_i_NP=PCa*PhiCaL_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
ICaL_i_CaMK=PCap*PhiCaL_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);

ICaNa_junc_NP=PCaNa*PhiCaNa_junc*d*(f*(1.0-nca)+jca*fca*nca);
ICaNa_junc_CaMK=PCaNap*PhiCaNa_junc*d*(fp*(1.0-nca)+jca*fcap*nca);
ICaNa_i_NP=PCaNa*PhiCaNa_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
ICaNa_i_CaMK=PCaNap*PhiCaNa_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);

ICaK_junc_NP=PCaK*PhiCaK_junc*d*(f*(1.0-nca)+jca*fca*nca);
ICaK_junc_CaMK=PCaKp*PhiCaK_junc*d*(fp*(1.0-nca)+jca*fcap*nca);

ICaK_i_NP=PCaK*PhiCaK_i*d*(f*(1.0-nca_i)+jca*fca*nca_i);
ICaK_i_CaMK=PCaKp*PhiCaK_i*d*(fp*(1.0-nca_i)+jca*fcap*nca_i);

% population combination and adding the pure
% CDI modulation
ICaL_junc = ((1-fICaL_CaMKonly)*ICaL_junc_NP + fICaL_CaMKonly*ICaL_junc_CaMK) * ical_pureCDI_junc;
ICaNa_junc = ((1-fICaL_CaMKonly)*ICaNa_junc_NP + fICaL_CaMKonly*ICaNa_junc_CaMK ) * ical_pureCDI_junc;
ICaK_junc = ((1-fICaL_CaMKonly)*ICaK_junc_NP + fICaL_CaMKonly*ICaK_junc_CaMK) * ical_pureCDI_junc;
%
ICaL_sl = ((1-fICaL_CaMKonly)*ICaL_i_NP + fICaL_CaMKonly*ICaL_i_CaMK) * ical_pureCDI_sl;
ICaNa_sl = ((1-fICaL_CaMKonly)*ICaNa_i_NP + fICaL_CaMKonly*ICaNa_i_CaMK) * ical_pureCDI_sl;
ICaK_sl = ((1-fICaL_CaMKonly)*ICaK_i_NP + fICaL_CaMKonly*ICaK_i_CaMK) * ical_pureCDI_sl;

% And we weight ICaL (in ss) and ICaL_i
ICaL_sl = ICaL_sl * (1-ICaL_fractionSS);
ICaNa_sl = ICaNa_sl * (1-ICaL_fractionSS);
ICaK_sl = ICaK_sl * (1-ICaL_fractionSS);
ICaL_junc = ICaL_junc * ICaL_fractionSS;
ICaNa_junc = ICaNa_junc * ICaL_fractionSS;
ICaK_junc = ICaK_junc * ICaL_fractionSS;

end


function [I_to, I_tof, I_tos, d_ito_xtos, d_ito_ytos, d_ito_xtof, d_ito_ytof, d_ito_xtos_p, d_ito_ytos_p, d_ito_xtof_p, d_ito_ytof_p] = getIto(membrane_v, ito_xtos, ito_ytos, ito_xtof, ito_ytof, ito_xtos_p, ito_ytos_p, ito_xtof_p, ito_ytof_p, ek, fItop, cellType, Itos_Multiplier, Itof_Multiplier)
% Fast and slow components of Ito
if cellType == 1 % EPI
    GtoSlow=0.02036 * Itos_Multiplier; 
    GtoFast=0.29856 * Itof_Multiplier; 
elseif cellType == 2 % MID
    GtoSlow= 0.04632 * Itos_Multiplier; 
    GtoFast=0.14928 * Itof_Multiplier;
elseif cellType == 0 % ENDO
    GtoSlow=0.07210 * Itos_Multiplier; 
    GtoFast=0.01276 * Itof_Multiplier; 
end

xtoss = 1/(1+exp(-(membrane_v-19.0)/13));
ytoss = 1/(1+exp((membrane_v+19.5)/5));
tauxtos = 9/(1+exp((membrane_v+3.0)/15))+0.5;
tauytos = 800/(1+exp((membrane_v+60.0)/10))+30;
d_ito_xtos = (xtoss - ito_xtos)/tauxtos;
d_ito_ytos = (ytoss - ito_ytos)/tauytos;

tauxtof = 8.5*exp(-((membrane_v+45)/50)^2)+0.5;
tauytof = 85*exp((-(membrane_v+40)^2/220))+7;
d_ito_xtof = (xtoss - ito_xtof)/tauxtof;
d_ito_ytof = (ytoss - ito_ytof)/tauytof;

% CaMKII effect on ss activation
xtoss_p = 1/(1+exp(-(membrane_v-29.0)/13));
d_ito_xtos_p = (xtoss_p - ito_xtos_p)/tauxtos;
d_ito_xtof_p = (xtoss_p - ito_xtof_p)/tauxtof;

% And CaMKII effect on ss inactivation
dti_develop=1.354+1.0e-4/(exp((membrane_v-167.4)/15.89)+exp(-(membrane_v-12.23)/0.2154));
dti_recover=1.0-0.5/(1.0+exp((membrane_v+70.0)/20.0));
tauytos_p = tauytos * dti_develop * dti_recover;
tauytof_p = tauytof * dti_develop * dti_recover;

d_ito_ytos_p = (ytoss - ito_ytos_p)/tauytos_p;
d_ito_ytof_p = (ytoss - ito_ytof_p)/tauytof_p;

I_tos = GtoSlow*(membrane_v-ek) * ((1-fItop)*ito_xtos*ito_ytos + fItop*ito_xtos_p*ito_ytos_p);    % [uA/uF]
I_tof = GtoFast*(membrane_v-ek)*((1-fItop) * ito_xtof*ito_ytof + fItop * ito_xtof_p*ito_ytof_p);
I_to = I_tos + I_tof;
end

function [I_K1] = getIK1_CRLP(v,  ko , EK, celltype, IK1_Multiplier)
aK1 = 4.094/(1+exp(0.1217*(v-EK-49.934)));
bK1 = (15.72*exp(0.0674*(v-EK-3.257))+exp(0.0618*(v-EK-594.31)))/(1+exp(-0.1629*(v-EK+14.207)));
K1ss = aK1/(aK1+bK1);

GK1=IK1_Multiplier  * 0.6992;
if celltype==1
    GK1=GK1*1.1;
elseif celltype==2
    GK1=GK1*1.3;
end
I_K1=GK1*sqrt(ko/5)*K1ss*(v-EK);
end


% Based on the Lu-Vandenberg model
function [IKr, dc0, dc1, dc2, do, di ] = getIKr(membrane_v,c0,c1, c2, o, i,...
    ko, EK, vfrt, celltype, IKr_Multiplier)

% transition rates
% from c0 to c1 in l-v model,
alpha = 0.1161 * exp(0.2990 * vfrt);
% from c1 to c0 in l-v/
beta =  0.2442 * exp(-1.604 * vfrt);

% from c1 to c2 in l-v/
alpha1 = 1.25 * 0.1235 ;
% from c2 to c1 in l-v/
beta1 =  0.1911;

% from c2 to o/           c1 to o
alpha2 =0.0578 * exp(0.9710 * vfrt); 
% from o to c2/
beta2 = 0.349e-3* exp(-1.062 * vfrt); 

% from o to i
alphai = 0.2533 * exp(0.5953 * vfrt); 
% from i to o
betai = 0.04568 * exp(-0.8209 * vfrt); 

% from c2 to i (from c1 in orig)
alphac2ToI = 0.52e-4 * exp(1.525 * vfrt); 
% from i to c2
betaItoC2 = (beta2 * betai * alphac2ToI)/(alpha2 * alphai); %
% for reason of backward compatibility of naming of an older version of a
% MM IKr, c3 in code is c0 in article diagram, c2 is c1, c1 is c2.

dc0 = c1 * beta - c0 * alpha; % delta for c0
dc1 = c0 * alpha + c2*beta1 - c1*(beta+alpha1); % c1
dc2 = c1 * alpha1 + o*beta2 + i*betaItoC2 - c2 * (beta1 + alpha2 + alphac2ToI); % subtraction is into c2, to o, to i. % c2
do = c2 * alpha2 + i*betai - o*(beta2+alphai);
di = c2*alphac2ToI + o*alphai - i*(betaItoC2 + betai);

GKr = 0.043 * sqrt(ko/5) * IKr_Multiplier; 
if celltype==1
    GKr=GKr*1.25; 
elseif celltype==2
    GKr=GKr*0.7; 
end

IKr = GKr * o  * (membrane_v-EK);
end

function [I_Ks, d_iks_xs_junc, d_iks_xs_sl ] = getIKs_Bartos(membrane_v, iks_xs_junc, iks_xs_sl, ions_ca_junc, ions_ca_sl, eks, Fjunc, Fsl, fIKs_PKA, cellType, IKs_Multiplier)
kPKA_Iks = fIKs_PKA;
gks_multiplier = 2.97002 * IKs_Multiplier;

if (cellType == 1)
    gks_multiplier = 1.4 * gks_multiplier;
elseif (cellType == 2)
    gks_multiplier = 0.5 * gks_multiplier;
end

% Modulation of I_Ks by beta-AR stimulation
gks_factor = 0.01;
P_g_0 = gks_factor*(0.2+0.2*kPKA_Iks);
P_g_max = gks_factor*(0.8+7*kPKA_Iks);
P_vh_0 = -1-10*kPKA_Iks;
P_vh_max = -12-9*kPKA_Iks;
P_tau_0 = 26+9*kPKA_Iks;
P_tau_max = 40+4*kPKA_Iks;

caks_junc = ions_ca_junc;
caks_sl = ions_ca_sl; 

gks_junc = P_g_0 + (P_g_max-P_g_0)/(1 + (150e-6/caks_junc)^1.3); % Regulated by PKA
gks_sl = P_g_0 + (P_g_max-P_g_0)/(1 + (150e-6/caks_sl)^1.3); % Regulated by PKA
VsXs_Ca_junc = P_vh_0 + (P_vh_max-P_vh_0)/(1 + (350e-6/caks_junc)^4); % Regulated by PKA
xsss_junc = 1/(1+exp(-(membrane_v-VsXs_Ca_junc)/25));
VsTs_Ca_junc = P_tau_0 + (P_tau_max-P_tau_0)/(1 + (150e-6/caks_junc)^3); % Regulated by PKA
tauxs_junc = 2*(50+(50+350*exp(-((membrane_v+30)^2)/4000))*1/(1+exp(-(membrane_v+VsTs_Ca_junc)/10)));
VsXs_Ca_sl = P_vh_0 + (P_vh_max-P_vh_0)/(1 + (350e-6/caks_sl)^4); % Regulated by PKA
xsss_sl = 1/(1+exp(-(membrane_v-VsXs_Ca_sl)/25));
VsTs_Ca_sl = P_tau_0 + (P_tau_max-P_tau_0)/(1 + (150e-6/caks_sl)^3); % Regulated by PKA
tauxs_sl = 2*(50+(50+350*exp(-((membrane_v+30)^2)/4000))*1/(1+exp(-(membrane_v+VsTs_Ca_sl)/10)));
d_iks_xs_junc = (xsss_junc-iks_xs_junc)/tauxs_junc;
d_iks_xs_sl = (xsss_sl-iks_xs_sl)/tauxs_sl;

I_ks_junc = Fjunc*gks_multiplier*gks_junc*iks_xs_junc^2*(membrane_v-eks);
I_ks_sl = Fsl*gks_multiplier*gks_sl*iks_xs_sl^2*(membrane_v-eks);

I_Ks = I_ks_junc+I_ks_sl;
end


function [ I_NCX_sl, I_NCX_junc] = getINaCa_ORd2011(v,F,R,T, na_ss, na_sl, na_i, na_o, ca_ss, ca_sl, ca_i, ca_o, celltype, INaCa_Multiplier, INaCa_fraction_junc)
% ORd-based NCX with modifications for better response to intracellular Na, and
% non-local driving force
zca = 2.0;
kna1=11.9712;
kna2=2.76;
kna3=88.767;
kasymm=19.4258;
wna= 3.2978e+04;
wca = 5.1756e+04;
wnaca= 2.7763e+03;
kcaon= 3.4164e+06;
kcaoff= 3.8532e+03;
qna=0.6718;
qca=0.0955;

hca=exp((qca*(v-8.3117)*F)/(R*T));
hna=exp((qna*(v-8.3117)*F)/(R*T));

h1=1+na_i/kna3*(1+hna);
h2=(na_i*hna)/(kna3*h1);
h3=1.0/h1;
h4=1.0+na_i/kna1*(1+na_i/kna2);
h5=na_i*na_i/(h4*kna1*kna2);
h6=1.0/h4;
h7=1.0+na_o/kna3*(1.0+1.0/hna);
h8=na_o/(kna3*hna*h7);
h9=1.0/h7;
h10=kasymm+1.0+na_o/kna1*(1.0+na_o/kna2);
h11=na_o*na_o/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*ca_o*kcaon;
k2=kcaoff;
k3p=h9*wca;
k3pp=h8*wnaca;
k3=k3p+k3pp;
k4p=h3*wca/hca;
k4pp=h2*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6*ca_i*kcaon;
k7=h5*h2*wna;
k8=h8*h11*wna;
x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
KmCaAct=150e-6;
allo=1.0/(1.0+(KmCaAct/ca_sl)^2);
zna=1.0;
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;
Gncx= 0.00179* INaCa_Multiplier;
if celltype==1
    Gncx=Gncx*1.1;
elseif celltype==2
    Gncx=Gncx*1.4;
end
I_NCX_sl=(1-INaCa_fraction_junc)*Gncx*allo*(zna*JncxNa+zca*JncxCa);

%calculate INaCa_junc
h1=1+na_sl/kna3*(1+hna);
h2=(na_sl*hna)/(kna3*h1);
h3=1.0/h1;
h4=1.0+na_sl/kna1*(1+na_sl/kna2);
h5=na_sl*na_sl/(h4*kna1*kna2);
h6=1.0/h4;
h7=1.0+na_o/kna3*(1.0+1.0/hna);
h8=na_o/(kna3*hna*h7);
h9=1.0/h7;
h10=kasymm+1.0+na_o/kna1*(1+na_o/kna2);
h11=na_o*na_o/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*ca_o*kcaon;
k2=kcaoff;
k3p=h9*wca;
k3pp=h8*wnaca;
k3=k3p+k3pp;
k4p=h3*wca/hca;
k4pp=h2*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6*ca_sl*kcaon;
k7=h5*h2*wna;
k8=h8*h11*wna;
x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
KmCaAct=150e-6;
allo=1.0/(1.0+(KmCaAct/ca_ss)^2);
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;
I_NCX_junc=INaCa_fraction_junc*Gncx*allo*(zna*JncxNa+zca*JncxCa);
end


function [J_SRCarel, J_SRleak, Jrel_ICaLdep, d_ryr_R, d_ryr_O, d_ryr_I, d_ryr_CaRI, d_ryr_R_p, d_ryr_O_p, d_ryr_I_p, d_ryr_CaRI_p, d_jrel_icaldep_act, d_jrel_icaldep_f1, d_jrel_icaldep_f2] = getJrel(ions_ca_junc, ions_ca_SR, I_Ca_junc, ryr_R, ryr_O, ryr_I, ryr_CaRI, ryr_R_p, ryr_O_p, ryr_I_p, jrel_icaldep_act, jrel_icaldep_f1, jrel_icaldep_f2, ryr_CaRI_p, camk_f_RyR, Jrel_Multiplier, extraParams)
% SR Ca release 
ks = 26.6 * Jrel_Multiplier;                 
koCa = 23.87221;               
kom = 0.16219;              
kiCa = 0.39871;              
kim = 0.04311;             
ec50SR = 0.75385;           
steepnessCaSR = 5.09473;
caExpFactor = 0.68655 ;
caTransFactor = 0.94428 ;
caExpFactor2 = 2.06273 ;
caTransFactor2 = 0.52967 ;


%% Directly ICaL-coupled RyRs
directRelMidpoint = 0.95271;
bt=12.47670;
a_rel=1.25 * bt ;
I_Ca_junc_positive = abs(I_Ca_junc);
I_Ca_junc_sigmoided = 1 - 1/(1+(I_Ca_junc_positive/0.45)^4.5);
Jrel_inf=a_rel*(I_Ca_junc_sigmoided)/(1.0+(directRelMidpoint/ions_ca_SR)^(7.72672));
tau_rel=bt/(1.0+0.0123/ions_ca_SR);

if tau_rel<0.001
    tau_rel=0.001;
end

d_jrel_icaldep_act = (Jrel_inf-jrel_icaldep_act)/tau_rel;

% and inactivation
tauInact = 64.11202;
Jrel_inact_inf = 1/(1 + (I_Ca_junc_sigmoided/1e-3));
d_jrel_icaldep_f1 = (Jrel_inact_inf - jrel_icaldep_f1)/tauInact;

% and slower inactivation
tauInact2 = 119.48978;
Jrel_inact_inf2 = 1/(1 + (I_Ca_junc_sigmoided/0.6e-3));
d_jrel_icaldep_f2 = (Jrel_inact_inf2 - jrel_icaldep_f2)/tauInact2;

Jrel_ICaLdep = (0.00174 * jrel_icaldep_act * jrel_icaldep_f1 * jrel_icaldep_f2 );

%% Main Ca-sensitive RyRs
MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/ions_ca_SR)^steepnessCaSR);
koSRCa = koCa/kCaSR;
kiSRCa = kiCa*kCaSR;

ecCaI = 0.001;
steepnessCaI = 5.93447; 
minCaI = 0.93249;
maxCaI = 30.13294;

baseRateCaI = 3.02320e-04;
sigmoidBaseCaI = minCaI + (maxCaI-minCaI)/(1+(ecCaI/ions_ca_junc)^steepnessCaI);
RI_to_CI = baseRateCaI * sigmoidBaseCaI;
CI_to_RI = 0.00248;

% not phosphorylated by CaMKII
RIcleft = 1 - ryr_R - ryr_O - ryr_I - ryr_CaRI;
d_ryr_R = (kim*RIcleft-kiSRCa*(caTransFactor*ions_ca_junc^caExpFactor)*ryr_R)-(caTransFactor2*koSRCa*ions_ca_junc^caExpFactor2*ryr_R-kom*ryr_O);   % R
d_ryr_O = (caTransFactor2*koSRCa*ions_ca_junc^caExpFactor2*ryr_R-kom*ryr_O)-(kiSRCa*(caTransFactor*ions_ca_junc^caExpFactor)*ryr_O-kim*ryr_I);% O
d_ryr_I = (kiSRCa*(caTransFactor*ions_ca_junc^caExpFactor)*ryr_O-kim*ryr_I)-(kom*ryr_I-caTransFactor2*koSRCa*ions_ca_junc^caExpFactor2*RIcleft);   % I
d_ryr_CaRI =  RI_to_CI*RIcleft - CI_to_RI * ryr_CaRI; % shift of 29 state numbers

J_SRCarel_np = ks*ryr_O*(ions_ca_SR-ions_ca_junc) + Jrel_ICaLdep;          % [mM/ms]

% And also a version of phosphorylated
caTransFactor2 = caTransFactor2 * 1.5; % *1.69
RIcleftP = 1 - ryr_R_p - ryr_O_p - ryr_I_p - ryr_CaRI_p;
d_ryr_R_p = (kim*RIcleftP-kiSRCa*(caTransFactor*ions_ca_junc^caExpFactor)*ryr_R_p)-(caTransFactor2*koSRCa*ions_ca_junc^caExpFactor2*ryr_R_p-kom*ryr_O_p);   % R
d_ryr_O_p = (caTransFactor2*koSRCa*ions_ca_junc^caExpFactor2*ryr_R_p-kom*ryr_O_p)-(kiSRCa*(caTransFactor*ions_ca_junc^caExpFactor)*ryr_O_p-kim*ryr_I_p);% O
d_ryr_I_p = (kiSRCa*(caTransFactor*ions_ca_junc^caExpFactor)*ryr_O_p-kim*ryr_I_p)-(kom*ryr_I_p-caTransFactor2*koSRCa*ions_ca_junc^caExpFactor2*RIcleftP);   % I
d_ryr_CaRI_p =  RI_to_CI*RIcleftP - CI_to_RI * ryr_CaRI_p;

J_SRCarel_p = ks*ryr_O_p*(ions_ca_SR-ions_ca_junc) + Jrel_ICaLdep;          % [mM/ms]

% Total release
J_SRCarel = J_SRCarel_p * camk_f_RyR + J_SRCarel_np * (1 - camk_f_RyR);

% Additional leak
nonlinearModifier = 0.2144*exp(1.83*ions_ca_SR);  % Leak should be nonlinear with load
CaMKIILeakMultiplier = 1 + 2*camk_f_RyR;   % And it should be promoted by CaMKII
J_SRleak = 1.59306e-06*(ions_ca_SR-ions_ca_junc) * nonlinearModifier * CaMKIILeakMultiplier;           % [mM/ms]
end

function J_serca = getJup(ions_ca_SR, ions_ca_i, Temp, camk_f_PLB, fPLB_PKA, casig_SERCA_act, cellType, Jup_Multiplier, extraParams)
% SR Ca uptake
Qpow = (Temp-310)/10;
Q10SRCaP = 2.6;          
Vmax_SRCaP = 0.00543 * Jup_Multiplier;  

Kmr = 2.31442;               
hillSRCaP =  1.02809;       
Kmf = 0.30672e-03;          

if (cellType == 1)
    Vmax_SRCaP = Vmax_SRCaP * 1.2; 
end

% PLB phosphorylation effect on affinity
phosphorylationTotal = camk_f_PLB + fPLB_PKA - camk_f_PLB*fPLB_PKA; % we assume the same effect, just making sure we don't count it twice.
Kmf_Phospho = Kmf * 0.5; %

% Direct Ca-based modulation
Km_SERCA_Ca = 0.4; 
Max_Vmax_SERCA_Ca = 1.11142;
Vmax_mult = 1 + Max_Vmax_SERCA_Ca / (1 + (Km_SERCA_Ca / casig_SERCA_act)^2);

J_serca_np = Q10SRCaP^Qpow*Vmax_SRCaP*Vmax_mult*((ions_ca_i/Kmf)^hillSRCaP-(ions_ca_SR/Kmr)^hillSRCaP)...
    /(1+(ions_ca_i/Kmf)^hillSRCaP+(ions_ca_SR/Kmr)^hillSRCaP);
J_serca_p = Q10SRCaP^Qpow*Vmax_SRCaP*Vmax_mult*((ions_ca_i/Kmf_Phospho)^hillSRCaP-(ions_ca_SR/Kmr)^hillSRCaP)...
    /(1+(ions_ca_i/Kmf_Phospho)^hillSRCaP+(ions_ca_SR/Kmr)^hillSRCaP);
J_serca = J_serca_np * (1 - phosphorylationTotal) + J_serca_p * phosphorylationTotal;
end

function [Ta, d_contraction_XS, d_contraction_XW, d_contraction_Ca_TRPN, d_contraction_TmBlocked, d_contraction_ZETAS, d_contraction_ZETAW] = getContractionLand(ions_ca_i, contraction_XS, contraction_XW, contraction_Ca_TRPN, contraction_TmBlocked, contraction_ZETAS, contraction_ZETAW, fTnI_PKA, fMyBPC_PKA, extraParams)
%==========================================================================
fracTnIpo = .0031;  % Derived quantity (TnI_PKAp(baseline)/TnItot)
fPKA_TnI = (1.45-0.45*(1-fTnI_PKA)/(1-fracTnIpo)); % multiplier for Ca unbinding from troponin.

PKAForceMultiplier = 1 + fMyBPC_PKA * 0.26;

% input coded here:
mode='intact';
lambda=1;
lambda_rate=0;

%EC parameters
perm50=0.35;
TRPN_n=1.65;
koff= 0.07854;
dr=0.25;
wfrac=0.5;
TOT_A=25;
ktm_unblock= 0.02626;
beta_1=-2.4;
beta_0=2.3;
gamma=0.0085;
gamma_wu=0.615;
phi=2.23;

if strcmp(mode,'skinned')
    nperm=2.2;
    ca50=2.5;
    Tref=40.5;
    nu=1;
    mu=1;
else
    nperm=2.036;
    ca50=fPKA_TnI * 0.76450;
    Tref= 80; %120; % Jakub updated
    nu=10.15996;
    mu=3.94046;
end

k_ws=0.004 * (1 + fMyBPC_PKA/2);
k_uw=0.026 ;

lambda_min=0.87;
lambda_max=1.2;

k_ws=k_ws*mu;
k_uw=k_uw*nu;

cdw=phi*k_uw*(1-dr)*(1-wfrac)/((1-dr)*wfrac);
cds=phi*k_ws*(1-dr)*wfrac/dr;
k_wu=k_uw*(1/wfrac-1)-k_ws;
k_su=k_ws*(1/dr-1)*wfrac;
A=(0.25*TOT_A)/((1-dr)*wfrac+dr)*(dr/0.25);

%XB model
lambda0=min(lambda_max,lambda);
Lfac=max(0,1+beta_0*(lambda0+min(lambda_min,lambda0)-(1+lambda_min)));

XU=(1-contraction_TmBlocked)-contraction_XW-contraction_XS; % unattached available xb = all - tm blocked - already prepowerstroke - already post-poststroke - no overlap
xb_ws=k_ws*contraction_XW ;
xb_uw=k_uw*XU* (1 + fMyBPC_PKA/2);
xb_wu=k_wu*contraction_XW;
xb_su=k_su*contraction_XS;

gamma_rate=gamma*max((contraction_ZETAS>0).*contraction_ZETAS,(contraction_ZETAS<-1).*(-contraction_ZETAS-1));
xb_su_gamma=gamma_rate*contraction_XS;
gamma_rate_w=gamma_wu*abs(contraction_ZETAW); % weak xbs don't like being strained
xb_wu_gamma=gamma_rate_w*contraction_XW;

d_contraction_XS=xb_ws-xb_su-xb_su_gamma;
d_contraction_XW=xb_uw-xb_wu-xb_ws-xb_wu_gamma;

ca50=ca50+beta_1*min(0.2,lambda-1);
d_contraction_Ca_TRPN=koff*(((ions_ca_i*1000)/ca50)^TRPN_n*(1-contraction_Ca_TRPN) - contraction_Ca_TRPN);

XSSS=dr*0.5;
XWSS=(1-dr)*wfrac*0.5;
ktm_block=ktm_unblock*(perm50^nperm)*0.5/(0.5-XSSS-XWSS);
d_contraction_TmBlocked=ktm_block*min(100,(contraction_Ca_TRPN^-(nperm/2)))*XU-ktm_unblock*(contraction_Ca_TRPN^(nperm/2))*contraction_TmBlocked;

%velocity dependence -- assumes distortion resets on W->S
d_contraction_ZETAS=A*lambda_rate-cds*contraction_ZETAS;% - gamma_rate * ZETAS;
d_contraction_ZETAW=A*lambda_rate-cdw*contraction_ZETAW;% - gamma_rate_w * ZETAW;


% Active Force
Ta=Lfac*PKAForceMultiplier * (Tref/dr)*((contraction_ZETAS+1)*contraction_XS+(contraction_ZETAW)*contraction_XW);
end




function [fICaLP,fIKsP,fPLBP,fTnIP,fINaP,fINaKP,fRyRP,fIKurP] = getEffectiveFraction(y, c)
%
% adapted from:  https://github.com/JQXGong/Ohara-beta-adrenergic
%
% Quantitative analysis of variability in an integrated model of
% human ventricular electrophysiology and B-adrenergic signaling
% by Jingqi Q.X. Gong, Monica E. Susilo, Anna Sher, Cynthia J. Musante, Eric A. Sobie
% DOI:https://doi.org/10.1016/j.yjmcc.2020.04.009
%


% Calculating effective fraction of phosphorylated substrates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ICaL
ICaLp = y(40);
fp_ICaL = (ICaLp + c(163)) / c(156); % Fraction of phosphorylated ICaL channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_ICaL = ifthenelse((fp_ICaL < 0.0) ,0.0001 , fp_ICaL);
fp_ICaL = ifthenelse((fp_ICaL > 1.0) ,0.9999 , fp_ICaL);

ical_f_hat_val = (fp_ICaL - c(166)) / (0.9273 - c(166));
ical_f_hat_val = ifthenelse((ical_f_hat_val < 0.0) , 0.0 , ical_f_hat_val); % Effective fraction of phosphorylated ICaL channels
ical_f_hat = ifthenelse((ical_f_hat_val > 1.0) , 1.0 , ical_f_hat_val); % Effective fraction of phosphorylated ICaL channels
%% IKs
IKsp = y(41);
fp_iks = (IKsp + c(145)) / c(144);% Fraction of phosphorylated IKs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_iks = ifthenelse((fp_iks < 0.0) ,0.0001 , fp_iks);
fp_iks = ifthenelse((fp_iks > 1.0) ,0.9999 , fp_iks);

iks_f_hat_val = (fp_iks - c(167)) / (0.7 - c(167));
iks_f_hat_val = ifthenelse((iks_f_hat_val < 0.0) , 0.0 , iks_f_hat_val);% Effective fraction of phosphorylated IKs channels
iks_f_hat = ifthenelse((iks_f_hat_val > 1.0) , 1.0 , iks_f_hat_val);% Effective fraction of phosphorylated IKs channels
%% Iup (PLB)
iup_f_plb = y(42);
% iup_f_pka_val = (iup_f_plb - 0.6591) / (0.9945 - 0.6591);
iup_f_pka_val = (iup_f_plb - 0.6662) / (0.9945 - 0.6662);
iup_f_pka_val = ifthenelse((iup_f_pka_val < 0.0) ,0.0 , iup_f_pka_val);
iup_f_pka = ifthenelse((iup_f_pka_val > 1.0) ,1.0 , iup_f_pka_val);
%% Tni
f_tni = y(43);
calcium_fhat_val = (f_tni - 0.67352) / (0.99918 - 0.67352);
calcium_fhat_val = ifthenelse((calcium_fhat_val < 0.0) , 0.0 , calcium_fhat_val);% Effective fraction of phosphorylated Troponin
calcium_fhat = ifthenelse((calcium_fhat_val > 1.0) , 1.0 , calcium_fhat_val);% Effective fraction of phosphorylated Troponin
%% INa
ina_f_ina = y(44);
ina_f_pka_val = (ina_f_ina - 0.23948) / (0.95014 - 0.23948);
% % ina_f_pka = ifthenelse((ina_f_pka_val < 0.0) , 0.0 , ina_f_pka_val);% Effective fraction of phosphorylated INa channels
ina_f_pka_val = ifthenelse((ina_f_pka_val < 0.0) , 0.0 , ina_f_pka_val);% Effective fraction of phosphorylated INa channels
ina_f_pka = ifthenelse((ina_f_pka_val > 1.0) , 1.0 , ina_f_pka_val);% Effective fraction of phosphorylated INa channels
%% INaK
f_inak = y(45);
inak_fhat_val = (f_inak - 0.12635) / (0.99801 - 0.12635);
inak_fhat_val = ifthenelse((inak_fhat_val < 0.0) , 0.0 , inak_fhat_val);% Effective fraction of phosphorylated INaK pumps
inak_fhat = ifthenelse((inak_fhat_val > 1.0) , 1.0 , inak_fhat_val);% Effective fraction of phosphorylated INaK pumps
%% RyR
RyRp = y(46);
fp_RyR = (RyRp + c(161)) / c(151);% Fraction of phosphorylated RyR channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_RyR = ifthenelse((fp_RyR < 0.0) ,0.0001 , fp_RyR);
fp_RyR = ifthenelse((fp_RyR > 1.0) ,0.9999 , fp_RyR);

irel_fhat_val = (fp_RyR - c(165)) / (0.9586 - c(165));
irel_fhat_val = ifthenelse((irel_fhat_val < 0.0) , 0.0 , irel_fhat_val);% Effective fraction of phosphorylated ryr channels
irel_fhat = ifthenelse((irel_fhat_val > 1.0) , 1.0 , irel_fhat_val);% Effective fraction of phosphorylated ryr channels
%% IKur
f_ikur = y(47);
ikur_fhat_val = (f_ikur -  5.89380e-02) / (0.39375 -  5.89380e-02);
ikur_fhat_val = ifthenelse((ikur_fhat_val < 0.0) , 0.0 , ikur_fhat_val);% Effective fraction of phosphorylated IKur channels
ikur_fhat = ifthenelse((ikur_fhat_val > 1.0) , 1.0 , ikur_fhat_val);% Effective fraction of phosphorylated IKur channels

% Applying the effective fraction of phosphorylated substrates into
% settings struct...  Note: the commented equations below are from the
% code downloaded from Rudy website... The ones above are from Myokit

fICaLP = ical_f_hat;%min(max(((y(98) + sigdata.ICaL_AKAP_PKA) / sigdata.ICaL_tot - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)) / (0.9273 - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)), 0), 1);
fIKsP = iks_f_hat;%min(max(((y(99) + sigdata.IKs_AKAP_PKA) / sigdata.IKs_tot - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)) / (0.7850 - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)), 0), 1);
fPLBP = iup_f_pka;%min(max((y(100) - 6.59100e-001) / (9.94500e-001 - 6.59100e-001), 0), 1);
fTnIP = calcium_fhat;%min(max((y(101) - 6.73519e-001) / (9.99180e-001 - 6.73519e-001), 0), 1);
fINaP = ina_f_pka;%min(max((y(102) - 2.39480e-001) / (9.50143e-001 - 2.39480e-001), 0), 1);
fINaKP = inak_fhat;%min(max((y(103) - 1.26345e-001) / (9.98014e-001 - 1.26345e-001), 0), 1);
fRyRP = irel_fhat;%min(max(((y(104) + sigdata.RyR_AKAP_PKA) / sigdata.RyR_tot - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)) / (0.9586 - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)), 0), 1);
fIKurP = ikur_fhat;%min(max((y(105) - 5.89380e-002) / (3.93747e-001 - 5.89380e-002), 0), 1);

end



function ydot = getPKASignalling(y, c)  % Remove t and pace because they are not used....

% adapted from:  https://github.com/JQXGong/Ohara-beta-adrenergic
%
% Quantitative analysis of variability in an integrated model of
% human ventricular electrophysiology and B-adrenergic signaling
% by Jingqi Q.X. Gong, Monica E. Susilo, Anna Sher, Cynthia J. Musante, Eric A. Sobie
% DOI:https://doi.org/10.1016/j.yjmcc.2020.04.009
%

% Create derivatives vector
ydot = zeros(size(y,1), size(y,2));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of the Signaling state variables
beta_cav_Gs_aGTP = y(1);   % 1: Gs_aGTP_CAV
beta_eca_Gs_aGTP = y(2);   % 2: Gs_aGTP_ECAV
beta_cyt_Gs_aGTP = y(3);   % 3: Gs_a_GTP_CYT
beta_cav_Gs_bg = y(4);     % 4: Gs_bg_CAV
beta_eca_Gs_bg = y(5);     % 5: Gs_bg_ECAV
beta_cyt_Gs_bg = y(6);     % 6: Gs_bg_CYT
beta_cav_Gs_aGDP = y(7);   % 7: Gs_aGDP_CAV
beta_eca_Gs_aGDP = y(8);   % 8: Gs_aGDP_ECAV
beta_cyt_Gs_aGDP = y(9);   % 9: Gs_aGDP_CYT

cAMP_cav = y(10);   % 10: cAMP_CAVVV
cAMP_eca = y(11);   % 11: cAMP_ECAV
cAMP_cyt = y(12);  % 12: cAMP_CYT

beta_cav_Rb1_pka_tot = y(13);  % 13: R_pkap_tot_CAV
beta_eca_Rb1_pka_tot = y(14);  % 14: R_pkap_tot_ECAV
beta_cyt_Rb1_pka_tot = y(15);  % 15: R_pkap_tot_CYT
beta_cav_Rb1_grk_tot = y(16);  % 16: R_grkp_tot_CAV
beta_eca_Rb1_grk_tot = y(17);  % 17: R_grkp_tot_ECAV
beta_cyt_Rb1_grk_tot = y(18);  % 18: R_grkp_tot_CYT

pka_cav_ARC = y(19);   % 19: RLC_CAV
pka_cav_A2RC = y(20);  % 20: L2RC_CAV
pka_cav_A2R = y(21);   % 21: L2R_CAV
pka_cav_C = y(22);     % 22: C_CAV
pka_cav_PKIC = y(23);  % 23: PKI_CAV
pka_eca_ARC = y(24);   % 24: RLC_ECAV
pka_eca_A2RC = y(25);  % 25: L2RC_ECAV
pka_eca_A2R = y(26);   % 26: L2R_ECAV
pka_eca_C = y(27);     % 27: C_ECAV
pka_eca_PKIC = y(28);  % 28: PKI_ECAV
pka_cyt_ARC = y(29);   % 29: RLC_CYT
pka_cyt_A2RC = y(30);  % 30: L2RC_CYT
pka_cyt_A2R = y(31);   % 31: L2R_CYT
pka_cyt_C = y(32);     % 32: C_CYT
pka_cyt_PKIC = y(33);  % 33: PKI_CYT

PDE3_P_cav = y(34);    %34   34: PDE3_P_CAV
PDE3_P_cyt = y(35);    %35   35: PDE3_P_CYT
PDE4_P_cav = y(36);    %36   36: PDE4_P_CAV
PDE4_P_eca = y(37);    %37   37: PDE4_P_ECAV
PDE4_P_cyt = y(38);    %38   38: PDE4_P_CYT

inhib1_p = y(39);  %39       39: Inhib1_P_CYT

ICaLp = y(40);     % 40: fLCC_P
IKsp = y(41);      % 41: fIKS_P
iup_f_plb = y(42); % 42: fPLB_P
f_tni = y(43);     % 43: fTnI_P
ina_f_ina = y(44); % 44: fINa_P
f_inak = y(45);    % 45: fINaK_P
RyRp = y(46);      % 46: fRyR_P
f_ikur = y(47);    % 47: fIKur_P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% update here to make the original PKA fractions bound between [0,1]
%%%% use 0.0001 and 0.9999
%%%% for AKAP related 3 channels implement the constraint in EffectiveFraction function
%%% ICaLp IKsp RyRp

%%% iup_f_plb
if iup_f_plb < 0.0
    iup_f_plb=0.0001;
elseif iup_f_plb > 1.0
    iup_f_plb= 0.9999;
end

%%% f_tni
if f_tni  < 0.0
    f_tni =0.0001;
elseif f_tni  > 1.0
    f_tni = 0.9999;
end

%%% ina_f_ina
if ina_f_ina  < 0.0
    ina_f_ina =0.0001;
elseif ina_f_ina > 1.0
    ina_f_ina = 0.9999;
end

%%% f_inak
if f_inak  < 0.0
    f_inak =0.0001;
elseif f_inak > 1.0
    f_inak = 0.9999;
end

%%% f_ikur
if f_ikur  < 0.0
    f_ikur =0.0001;
elseif f_ikur > 1.0
    f_ikur= 0.9999;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta_cav_Rb2_pka_tot = y(48);  % 48: Rb2_pkap_tot_CAV
beta_cav_Rb2_grk_tot = y(49);  % 49: Rb2_grkp_tot_CAV
beta_cav_Gi_aGTP = y(50);      % 50: Gi_aGTP_CAV
beta_cav_Gi_bg = y(51);        % 51: Gi_bg_CAV
beta_cav_Gi_aGDP = y(52);      % 52: Gi_aGDP_CAV
beta_eca_Rb2_pka_tot = y(53);  % 53: Rb2_pkap_tot_ECAV
beta_eca_Rb2_grk_tot = y(54);  % 54: Rb2_grkp_tot_ECAV
beta_eca_Gi_aGTP = y(55);      % 55: Gi_aGTP_ECAV
beta_eca_Gi_bg = y(56);        % 56: Gi_bg_ECAV
beta_eca_Gi_aGDP = y(57);      % 57: Gi_aGDP_ECAV

%%

%% CAVEOLAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total concentration of non-phosphorylated B1AR in the caveolar subspace
beta_cav_Rb1_np_tot = c(87) - beta_cav_Rb1_pka_tot - beta_cav_Rb1_grk_tot;
% Total concentration of non-phosphorylated B2AR in the caveolar subspace
beta_cav_Rb2_np_tot = c(89) - beta_cav_Rb2_pka_tot - beta_cav_Rb2_grk_tot;
% Concentration of Gi holoenzyme in the caveolar subspace
beta_cav_Gi_abg = c(63) * c(59) * c(5) - beta_cav_Gi_aGTP - beta_cav_Gi_aGDP;
% Concentration of Gs holoenzyme in the caveolar subspace
beta_cav_Gs_abg = c(61) * c(58) * c(5) - beta_cav_Gs_aGTP - beta_cav_Gs_aGDP;

beta_cav_Gs_f_d = beta_cav_Gs_abg * c(93) / c(92);
beta_cav_Gs_f_b = (c(94) + c(91)) / c(92) + beta_cav_Rb1_np_tot + beta_cav_Rb2_np_tot - beta_cav_Gs_abg;
beta_cav_Gs_f_c = (c(91) * (beta_cav_Rb1_np_tot - beta_cav_Gs_abg) + c(94) * (beta_cav_Rb2_np_tot - beta_cav_Gs_abg) + c(93)) / c(92);
beta_cav_Gs_f_rr = -beta_cav_Gs_f_d / 27.0 * beta_cav_Gs_f_b ^ 3.0 - beta_cav_Gs_f_b * beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_c / 108.0 + beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_d / 6.0 + beta_cav_Gs_f_c ^ 3.0 / 27.0 + beta_cav_Gs_f_d * beta_cav_Gs_f_d / 4.0;
beta_cav_Gs_f_yr = ifthenelse((beta_cav_Gs_f_rr > 0.0) , sqrt(beta_cav_Gs_f_rr) , 0.0) + beta_cav_Gs_f_d / 2.0 + beta_cav_Gs_f_b * beta_cav_Gs_f_c / 6.0 - beta_cav_Gs_f_b ^ 3.0 / 27.0;
beta_cav_Gs_f_yi = ifthenelse((beta_cav_Gs_f_rr < 0.0) , sqrt(-beta_cav_Gs_f_rr) , 0.0);
beta_cav_Gs_f_mag = (beta_cav_Gs_f_yr * beta_cav_Gs_f_yr + beta_cav_Gs_f_yi * beta_cav_Gs_f_yi) ^ (1.0 / 6.0);
beta_cav_Gs_f_arg = atan(beta_cav_Gs_f_yi / beta_cav_Gs_f_yr) / 3.0;
beta_cav_Gs_f_x = (beta_cav_Gs_f_c / 3.0 - beta_cav_Gs_f_b * beta_cav_Gs_f_b / 9.0) / (beta_cav_Gs_f_mag * beta_cav_Gs_f_mag);
beta_cav_Gs_f_r = beta_cav_Gs_f_mag * cos(beta_cav_Gs_f_arg) * (1.0 - beta_cav_Gs_f_x) - beta_cav_Gs_f_b / 3.0;
beta_cav_Gs_f_i = beta_cav_Gs_f_mag * sin(beta_cav_Gs_f_arg) * (1.0 + beta_cav_Gs_f_x);

% Concentration of free Gs in the caveolar subspace
beta_cav_Gs_f = sqrt(beta_cav_Gs_f_r * beta_cav_Gs_f_r + beta_cav_Gs_f_i * beta_cav_Gs_f_i);
% Concentration of free non-phosphorylated beta1AR in caveolar subspace
beta_cav_Rb1_f = beta_cav_Rb1_np_tot / (1.0 + c(1) / c(65) + beta_cav_Gs_f * (c(67) + c(1)) / (c(66) * c(67)));
% Concentration of non-phosphorylated Ligand / Receptor1 complexes in caveolar subspace
beta_cav_LRb1 = c(1) * beta_cav_Rb1_f / c(65);
% Concentration of non-phosphorylated Ligand / Receptor1 / G-protein complexes in caveolar subspace
beta_cav_LRb1Gs = c(1) * beta_cav_Rb1_f * beta_cav_Gs_f / (c(66) * c(67));
% Concentration of free non-phosphorylated beta2AR in caveolar subspace
beta_cav_Rb2_f = beta_cav_Rb2_np_tot / (1.0 + c(1) / c(72) + beta_cav_Gs_f * (c(69) + c(1)) / (c(71) * c(69)));
% Concentration of non-phosphorylated Ligand / Receptor2 complexes in caveolar subspace
beta_cav_LRb2 = c(1) * beta_cav_Rb2_f / c(72);

% Concentration of non-phosphorylated Receptor2 / G-protein complexes in caveolar subspace
beta_cav_Rb2Gs = beta_cav_Rb2_f * beta_cav_Gs_f / c(71);
% Concentration of non-phosphorylated Receptor1 / G-protein complexes in caveolar subspace
beta_cav_Rb1Gs = beta_cav_Rb1_f * beta_cav_Gs_f / c(66);
% Concentration of non-phosphorylated Ligand / Receptor2 / G-protein complexes in caveolar subspace
beta_cav_LRb2Gs = c(1) * beta_cav_Rb2_f * beta_cav_Gs_f / (c(71) * c(69));

% Concentration of total PKA-phosphorylated beta1 receptors
ydot(13) = 0.001 * (c(84) * pka_cav_C * beta_cav_Rb1_np_tot - c(85) * beta_cav_Rb1_pka_tot);
% Concentration of total GRK-phosphorylated beta1 receptors
ydot(16) = 0.001 * (c(83) * c(86) * (beta_cav_LRb1 + beta_cav_LRb1Gs) - c(82) * beta_cav_Rb1_grk_tot);
% Concentration of total PKA-phosphorylated beta2 receptors
ydot(48) = 0.001 * (c(84) * pka_cav_C * beta_cav_Rb2_np_tot - c(85) * beta_cav_Rb2_pka_tot);
% Concentration of total GRK-phosphorylated beta2 receptors
ydot(49) = 0.001 * (c(83) * c(86) * (beta_cav_LRb2 + beta_cav_LRb2Gs) - c(82) * beta_cav_Rb2_grk_tot);

%% EXTRACAVEOLAR %%%%%%%%%%%%%%%%%%
%
% beta_eca
%
% Concentration of Gs holoenzyme in the extracaveolar space
beta_eca_Gs_abg = c(60) * c(58) * c(6) - beta_eca_Gs_aGTP - beta_eca_Gs_aGDP;
% Total concentration of non-phosphorylated B2AR in the extracaveolar space
beta_eca_Rb2_np_tot = c(97) - beta_eca_Rb2_pka_tot - beta_eca_Rb2_grk_tot;
% Total concentration of non-phosphorylated B1AR in the extracaveolar space
beta_eca_Rb1_np_tot = c(98) - beta_eca_Rb1_pka_tot - beta_eca_Rb1_grk_tot;
beta_eca_Gs_f_d = beta_eca_Gs_abg * c(101) / c(102);
beta_eca_Gs_f_c = (c(103) * (beta_eca_Rb1_np_tot - beta_eca_Gs_abg) + c(100) * (beta_eca_Rb2_np_tot - beta_eca_Gs_abg) + c(101)) / c(102);
beta_eca_Gs_f_b = (c(100) + c(103)) / c(102) + beta_eca_Rb1_np_tot + beta_eca_Rb2_np_tot - beta_eca_Gs_abg;
beta_eca_Gs_f_rr = -beta_eca_Gs_f_d / 27.0 * beta_eca_Gs_f_b ^ 3.0 - beta_eca_Gs_f_b * beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_c / 108.0 + beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_d / 6.0 + beta_eca_Gs_f_c ^ 3.0 / 27.0 + beta_eca_Gs_f_d * beta_eca_Gs_f_d / 4.0;
beta_eca_Gs_f_yi = ifthenelse((beta_eca_Gs_f_rr < 0.0) , sqrt(-beta_eca_Gs_f_rr) , 0.0);
beta_eca_Gs_f_yr = ifthenelse((beta_eca_Gs_f_rr > 0.0) , sqrt(beta_eca_Gs_f_rr) , 0.0) + beta_eca_Gs_f_d / 2.0 + beta_eca_Gs_f_b * beta_eca_Gs_f_c / 6.0 - beta_eca_Gs_f_b ^ 3.0 / 27.0;
beta_eca_Gs_f_mag = (beta_eca_Gs_f_yr * beta_eca_Gs_f_yr + beta_eca_Gs_f_yi * beta_eca_Gs_f_yi) ^ (1.0 / 6.0);
beta_eca_Gs_f_arg = atan(beta_eca_Gs_f_yi / beta_eca_Gs_f_yr) / 3.0;
beta_eca_Gs_f_x = (beta_eca_Gs_f_c / 3.0 - beta_eca_Gs_f_b * beta_eca_Gs_f_b / 9.0) / (beta_eca_Gs_f_mag * beta_eca_Gs_f_mag);
beta_eca_Gs_f_i = beta_eca_Gs_f_mag * sin(beta_eca_Gs_f_arg) * (1.0 + beta_eca_Gs_f_x);
beta_eca_Gs_f_r = beta_eca_Gs_f_mag * cos(beta_eca_Gs_f_arg) * (1.0 - beta_eca_Gs_f_x) - beta_eca_Gs_f_b / 3.0;
% Concentration of free Gs in the caveolar subspace
beta_eca_Gs_f = sqrt(beta_eca_Gs_f_r * beta_eca_Gs_f_r + beta_eca_Gs_f_i * beta_eca_Gs_f_i);
% Concentration of free non-phosphorylated beta1AR in the extracaveolar space
beta_eca_Rb1_f = beta_eca_Rb1_np_tot / (1.0 + c(1) / c(65) + beta_eca_Gs_f * (c(67) + c(1)) / (c(66) * c(67)));
% Concentration of free non-phosphorylated beta2AR in the extracaveolar space
beta_eca_Rb2_f = beta_eca_Rb2_np_tot / (1.0 + c(1) / c(72) + beta_eca_Gs_f * (c(69) + c(1)) / (c(71) * c(69)));
% Concentration of non-phosphorylated Ligand / Receptor1 complexes in the extracaveolar space
beta_eca_LRb1 = c(1) * beta_eca_Rb1_f / c(65);
% Concentration of non-phosphorylated Ligand / Receptor2 complexes in the extracaveolar space
beta_eca_LRb2 = c(1) * beta_eca_Rb2_f / c(72);
% Concentration of non-phosphorylated Ligand / Receptor2 / G-protein complexes in the extracaveolar space
beta_eca_LRb2Gs = c(1) * beta_eca_Rb2_f * beta_eca_Gs_f / (c(71) * c(69));
% Concentration of non-phosphorylated Ligand / Receptor1 / G-protein complexes in the extracaveolar space
beta_eca_LRb1Gs = c(1) * beta_eca_Rb1_f * beta_eca_Gs_f / (c(66) * c(67));
% Concentration of non-phosphorylated Receptor2 / G-protein complexes in the extracaveolar space
beta_eca_Rb2Gs = beta_eca_Rb2_f * beta_eca_Gs_f / c(71);
% Concentration of non-phosphorylated Receptor1 / G-protein complexes in the extracaveolar space
beta_eca_Rb1Gs = beta_eca_Rb1_f * beta_eca_Gs_f / c(66);
beta_eca_RGs_tot = beta_eca_Rb1Gs + c(96) * beta_eca_Rb2Gs;
beta_eca_LRGs_tot = beta_eca_LRb1Gs + c(96) * beta_eca_LRb2Gs;
% Concentration of Gi holoenzyme in the extracaveolar space
beta_eca_Gi_abg = c(64) * c(59) * c(6) - beta_eca_Gi_aGTP - beta_eca_Gi_aGDP;


beta_eca_Rb2_pka_f_c = -beta_eca_Rb2_pka_tot * c(70) * c(73);
beta_eca_Rb2_pka_f_b = beta_eca_Gi_abg * (c(1) + c(73)) - beta_eca_Rb2_pka_tot * (c(73) + c(1)) + c(70) * c(73) * (1.0 + c(1) / c(68));
% Concentration of free PKA-phosphorylated beta2AR in the extracaveolar space
beta_eca_Rb2_pka_f = (-beta_eca_Rb2_pka_f_b + sqrt(beta_eca_Rb2_pka_f_b * beta_eca_Rb2_pka_f_b - 4.0 * c(99) * beta_eca_Rb2_pka_f_c)) / (2.0 * c(99));

% Concentration of total PKA-phosphorylated beta1 receptors in the extracaveolar space
ydot(14) = 0.001 * (c(84) * pka_eca_C * beta_eca_Rb1_np_tot - c(85) * beta_eca_Rb1_pka_tot);
% Concentration of total GRK-phosphorylated beta1 receptors in the extracaveolar space
ydot(17) = 0.001 * (c(83) * c(95) * (beta_eca_LRb1 + beta_eca_LRb1Gs) - c(82) * beta_eca_Rb1_grk_tot);
% Concentration of total PKA-phosphorylated beta2 receptors in the extracaveolar space
ydot(53) = 0.001 * (c(84) * pka_eca_C * beta_eca_Rb2_np_tot - c(85) * beta_eca_Rb2_pka_tot);
% Concentration of total GRK-phosphorylated beta2 receptors in the extracaveolar space
ydot(54) = 0.001 * (c(83) * c(95) * (beta_eca_LRb2 + beta_eca_LRb2Gs) - c(82) * beta_eca_Rb2_grk_tot);


%% CYTOPLASM %%%%%%%%%%%
% Concentration of Gs holoenzyme in the cytoplasm
beta_cyt_Gs_abg = c(62) * c(58) * c(7) - beta_cyt_Gs_aGTP - beta_cyt_Gs_aGDP;
% Total concentration of non-phosphorylated beta-1 AR in the cytoplasm
beta_cyt_Rb1_np_tot = c(104) - beta_cyt_Rb1_pka_tot - beta_cyt_Rb1_grk_tot;

beta_cyt_Rb1_np_f_b = beta_cyt_Gs_abg * (c(67) + c(1)) - beta_cyt_Rb1_np_tot * (c(67) + c(1)) + c(66) * c(67) * (1.0 + c(1) / c(65));
beta_cyt_Rb1_np_f_c = -beta_cyt_Rb1_np_tot * c(67) * c(66);
% Concentration of free non-phosphorylated beta-1 AR in the cytoplasm
Rb1_np_f = (-beta_cyt_Rb1_np_f_b + sqrt(beta_cyt_Rb1_np_f_b * beta_cyt_Rb1_np_f_b - 4.0 * c(106) * beta_cyt_Rb1_np_f_c)) / (2.0 * c(106));
% Concentration of free (non-complexed) Gi in the cytoplasm
beta_cyt_Gs_f = beta_cyt_Gs_abg / (1.0 + Rb1_np_f / c(66) * (1.0 + c(1) / c(67)));
% Concentration of non-phosphorylated ligand / beta-1 AR complexes in the cytoplasm
LRb1_np = c(1) * Rb1_np_f / c(65);
% Concentration of non-phosphorylated ligand / beta-1 AR / Gs complexes in the cytoplasm
LRb1Gs_np = c(1) * Rb1_np_f * beta_cyt_Gs_f / (c(66) * c(67));
% Concentration of non-phosphorylated beta-1 AR / Gs complexes in the cytoplasm
Rb1Gs_np = beta_cyt_Gs_f * Rb1_np_f / c(66);

% Concentration of total PKA-phosphorylated receptors in the cytoplasm
ydot(15) = 0.001 * (c(84) * pka_cyt_C * beta_cyt_Rb1_np_tot - c(85) * beta_cyt_Rb1_pka_tot);
% Concentration of total GRK-phosphorylated receptors in the cytoplasm
ydot(18) = 0.001 * (c(83) * c(105) * (LRb1_np + LRb1Gs_np) - c(82) * beta_cyt_Rb1_grk_tot);


%% function Mod_GprotAct() in the other version %%%

beta_cav_RGs_tot = beta_cav_Rb1Gs + c(88) * beta_cav_Rb2Gs;
beta_cav_LRGs_tot = beta_cav_LRb1Gs + c(88) * beta_cav_LRb2Gs;
beta_cav_Rb2_pka_f_c = -beta_cav_Rb2_pka_tot * c(70) * c(73);
beta_cav_Rb2_pka_f_b = beta_cav_Gi_abg * (c(1) + c(73)) - beta_cav_Rb2_pka_tot * (c(73) + c(1)) + c(70) * c(73) * (1.0 + c(1) / c(68));
% Concentration of free PKA-phosphorylated beta2AR in the caveolar subspace
beta_cav_Rb2_pka_f = (-beta_cav_Rb2_pka_f_b + sqrt(beta_cav_Rb2_pka_f_b * beta_cav_Rb2_pka_f_b - 4.0 * c(90) * beta_cav_Rb2_pka_f_c)) / (2.0 * c(90));
% Concentration of free (non-complexed) Gi in the caveolar subspace
beta_cav_Gi_f = beta_cav_Gi_abg / (1.0 + beta_cav_Rb2_pka_f / c(70) * (1.0 + c(1) / c(73)));
% Concentration of PKA phosphorylated b2AR/Gi complexes in caveolar subspace
beta_cav_Rb2Gi = beta_cav_Rb2_pka_f * beta_cav_Gi_f / c(70);
% Concentration of PKA phosphorylated Ligand/b2AR/Gi complexes in caveolar subspace
beta_cav_LRb2Gi = beta_cav_Rb2Gi * c(1) / c(73);
% Concentration of free (non-complexed) Gi in the extracaveolar space
beta_eca_Gi_f = beta_eca_Gi_abg / (1.0 + beta_eca_Rb2_pka_f / c(70) * (1.0 + c(1) / c(73)));
% Concentration of PKA phosphorylated b2AR/Gi complexes in extracaveolar space
beta_eca_Rb2Gi = beta_eca_Rb2_pka_f * beta_eca_Gi_f / c(70);
% Concentration of PKA phosphorylated Ligand/b2AR/Gi complexes in extracaveolar space
beta_eca_LRb2Gi = c(1) / c(73) * beta_eca_Rb2Gi;

% Concentration of active Gs alpha subunit in caveolar subspace (tri-phosphate)
ydot(1) = 0.001 * (c(75) * beta_cav_RGs_tot + c(74) * beta_cav_LRGs_tot - c(78) * beta_cav_Gs_aGTP);
% Concentration of active Gi alpha subunit in caveolar subspace
ydot(50) = 0.001 * (c(77) * beta_cav_Rb2Gi + c(76) * beta_cav_LRb2Gi - c(79) * beta_cav_Gi_aGTP);
% Concentration of active Gs alpha subunit in the extracaveolar space (tri-phosphate)
ydot(2) = 0.001 * (c(75) * beta_eca_RGs_tot + c(74) * beta_eca_LRGs_tot - c(78) * beta_eca_Gs_aGTP);
% Concentration of active Gi alpha subunit in the extracaveolar space
ydot(55) = 0.001 * (c(77) * beta_eca_Rb2Gi + c(76) * beta_eca_LRb2Gi - c(79) * beta_eca_Gi_aGTP);
% Concentration of active Gs alpha subunit in cytoplasm
ydot(3) = 0.001 * (c(75) * Rb1Gs_np + c(74) * LRb1Gs_np - c(78) * beta_cyt_Gs_aGTP);

% Concentration of active Gs beta-gamma subunit in caveolar subspace
ydot(4) = 0.001 * (c(75) * beta_cav_RGs_tot + c(74) * beta_cav_LRGs_tot - c(80) * beta_cav_Gs_bg * beta_cav_Gs_aGDP);
% Concentration of active Gi beta-gamma subunit in caveolar subspace
ydot(51) = 0.001 * (c(77) * beta_cav_Rb2Gi + c(76) * beta_cav_LRb2Gi - c(81) * beta_cav_Gi_bg * beta_cav_Gi_aGDP);
% Concentration of active Gs beta-gamma subunit in the extracaveolar space
ydot(5) = 0.001 * (c(75) * beta_eca_RGs_tot + c(74) * beta_eca_LRGs_tot - c(80) * beta_eca_Gs_bg * beta_eca_Gs_aGDP);
% Concentration of active Gi beta-gamma subunit in the extracaveolar space
ydot(56) = 0.001 * (c(77) * beta_eca_Rb2Gi + c(76) * beta_eca_LRb2Gi - c(81) * beta_eca_Gi_bg * beta_eca_Gi_aGDP);
% Concentration of active Gs beta-gamma subunit in cytoplasm
ydot(6) = 0.001 * (c(75) * Rb1Gs_np + c(74) * LRb1Gs_np - c(80) * beta_cyt_Gs_bg * beta_cyt_Gs_aGDP);

% Concentration of inactive Gs alpha subunit in caveolar subspace (di-phosphate)
ydot(7) = 0.001 * (c(78) * beta_cav_Gs_aGTP - c(80) * beta_cav_Gs_bg * beta_cav_Gs_aGDP);
% Concentration of inactive Gi alpha subunit in caveolar subspace
ydot(52) = 0.001 * (c(79) * beta_cav_Gi_aGTP - c(81) * beta_cav_Gi_bg * beta_cav_Gi_aGDP);
% Concentration of inactive Gs alpha subunit in the extracaveolar space (di-phosphate)
ydot(8) = 0.001 * (c(78) * beta_eca_Gs_aGTP - c(80) * beta_eca_Gs_bg * beta_eca_Gs_aGDP);
% Concentration of inactive Gi alpha subunit in the extracaveolar space
ydot(57) = 0.001 * (c(79) * beta_eca_Gi_aGTP - c(81) * beta_eca_Gi_bg * beta_eca_Gi_aGDP);
% Concentration of inactive Gs alpha subunit in cytoplasm
ydot(9) = 0.001 * (c(78) * beta_cyt_Gs_aGTP - c(80) * beta_cyt_Gs_bg * beta_cyt_Gs_aGDP);

%% function Mod_AC()

% only calculating constants, but does not compute derivative of state
% variables

%% function Mod_PKA() in the other version %%%%%%%%%%
% Concentration of free PKA RC subunits in the caveolar compartment
pka_cav_RCf = c(12) - pka_cav_ARC - pka_cav_A2RC - pka_cav_A2R;

% Caveolar concentration of PKA RC dimer with 1 cAMP molecule bound
ydot(19) = 0.001 * (c(19) * pka_cav_RCf * cAMP_cav - c(22) * pka_cav_ARC - c(20) * pka_cav_ARC * cAMP_cav + c(23) * pka_cav_A2RC);
% Caveolar concentration of PKA RC dimer with 2 cAMP molecules bound
ydot(20) = 0.001 * (c(20) * pka_cav_ARC * cAMP_cav - (c(23) + c(21)) * pka_cav_A2RC + c(24) * pka_cav_A2R * pka_cav_C);
% Caveolar concentration of PKA R subunit with 2 cAMP molecules bound
ydot(21) = 0.001 * (c(21) * pka_cav_A2RC - c(24) * pka_cav_A2R * pka_cav_C);
% Caveolar concentration of free PKA catalytic subunit
ydot(22) = 0.001 * (c(21) * pka_cav_A2RC - c(24) * pka_cav_A2R * pka_cav_C + c(18) * pka_cav_PKIC - c(17) * (c(14) - pka_cav_PKIC) * pka_cav_C);
% Caveolar concentration of free PKI inactivated PKA C subunit
ydot(23) = 0.001 * (c(17) * (c(14) - pka_cav_PKIC) * pka_cav_C - c(18) * pka_cav_PKIC);

% Concentration of free PKA RC subunits in the Extracaveolar compartment
pka_eca_RCf = c(11) - pka_eca_ARC - pka_eca_A2RC - pka_eca_A2R;
% Extracaveolar rate of change in free cAMP through binding by PKA
pka_eca_dcAMP = -c(19) * pka_eca_RCf * cAMP_eca + c(25) * pka_eca_ARC - c(20) * pka_eca_ARC * cAMP_eca + c(26) * pka_eca_A2RC;

% Extracaveolar concentration of PKA RC dimer with 1 cAMP molecule bound
ydot(24) = 0.001 * (c(19) * pka_eca_RCf * cAMP_eca - c(25) * pka_eca_ARC - c(20) * pka_eca_ARC * cAMP_eca + c(26) * pka_eca_A2RC);
% Extracaveolar concentration of PKA RC dimer with 2 cAMP molecules bound
ydot(25) = 0.001 * (c(20) * pka_eca_ARC * cAMP_eca - (c(26) + c(21)) * pka_eca_A2RC + c(27) * pka_eca_A2R * pka_eca_C);
% Extracaveolar concentration of PKA R subunit with 2 cAMP molecules bound
ydot(26) = 0.001 * (c(21) * pka_eca_A2RC - c(27) * pka_eca_A2R * pka_eca_C);
% Extracaveolar concentration of free PKA catalytic subunit
ydot(27) = 0.001 * (c(21) * pka_eca_A2RC - c(27) * pka_eca_A2R * pka_eca_C + c(18) * pka_eca_PKIC - c(17) * (c(15) - pka_eca_PKIC) * pka_eca_C);
% Extracaveolar concentration of free PKI inactivated PKA C subunit
ydot(28) = 0.001 * (c(17) * (c(15) - pka_eca_PKIC) * pka_eca_C - c(18) * pka_eca_PKIC);

% Concentration of free PKA RC subunits in the Cytosolic compartment
pka_cyt_RCf = c(13) - pka_cyt_ARC - pka_cyt_A2RC - pka_cyt_A2R;

% Cytosolic concentration of PKA RC dimer with 1 cAMP molecule bound
ydot(29) = 0.001 * (c(19) * pka_cyt_RCf * cAMP_cyt - c(28) * pka_cyt_ARC - c(20) * pka_cyt_ARC * cAMP_cyt + c(29) * pka_cyt_A2RC);
% Cytosolic concentration of PKA RC dimer with 2 cAMP molecules bound
ydot(30) = 0.001 * (c(20) * pka_cyt_ARC * cAMP_cyt - (c(29) + c(21)) * pka_cyt_A2RC + c(30) * pka_cyt_A2R * pka_cyt_C);
% Cytosolic concentration of PKA R subunit with 2 cAMP molecules bound
ydot(31) = 0.001 * (c(21) * pka_cyt_A2RC - c(30) * pka_cyt_A2R * pka_cyt_C);
% Cytosolic concentration of free PKA catalytic subunit
ydot(32) = 0.001 * (c(21) * pka_cyt_A2RC - c(30) * pka_cyt_A2R * pka_cyt_C + c(18) * pka_cyt_PKIC - c(17) * (c(16) - pka_cyt_PKIC) * pka_cyt_C);
% Cytosolic concentration of free PKI inactivated PKA C subunit
ydot(33) = 0.001 * (c(17) * (c(16) - pka_cyt_PKIC) * pka_cyt_C - c(18) * pka_cyt_PKIC);

%% function Mod_cAMP() in the other version


% Caveolar rate of change in free cAMP through binding by PKA
pka_cav_dcAMP = -c(19) * pka_cav_RCf * cAMP_cav + c(22) * pka_cav_ARC - c(20) * pka_cav_ARC * cAMP_cav + c(23) * pka_cav_A2RC;

% Cytosolic rate of change in free cAMP through binding by PKA
pka_cyt_dcAMP = -c(19) * pka_cyt_RCf * cAMP_cyt + c(28) * pka_cyt_ARC - c(20) * pka_cyt_ARC * cAMP_cyt + c(29) * pka_cyt_A2RC;

%PDE
% Rate of cAMP degradation by PDE2 in cytosolic subspace
dcAMP_PDE2_cyt = c(52) * c(41) / (1.0 + c(44) / cAMP_cyt);
% Rate of cAMP degradation by PDE2 in extracaveolar subspace
dcAMP_PDE2_eca = c(51) * c(41) / (1.0 + c(44) / cAMP_eca);
% Rate of cAMP degradation by PDE2 in caveolar subspace
dcAMP_PDE2_cav = c(50) * c(41) / (1.0 + c(44) / cAMP_cav);
% Rate of cAMP degradation by PDE3 in caveolar subspace
dcAMP_PDE3_cav = (c(53) + (c(47) - 1.0) * PDE3_P_cav) * c(42) / (1.0 + c(45) / cAMP_cav);
% Rate of cAMP degradation by PDE4 in cytosolic subspace
dcAMP_PDE4_cyt = (c(57) + (c(47) - 1.0) * PDE4_P_cyt) * c(43) / (1.0 + c(46) / cAMP_cyt);
% Rate of cAMP degradation by PDE4 in extracaveolar subspace
dcAMP_PDE4_eca = (c(56) + (c(47) - 1.0) * PDE4_P_eca) * c(43) / (1.0 + c(46) / cAMP_eca);
% Rate of cAMP degradation by PDE4 in caveolar subspace
dcAMP_PDE4_cav = (c(55) + (c(47) - 1.0) * PDE4_P_cav) * c(43) / (1.0 + c(46) / cAMP_cav);
% Rate of cAMP degradation by PDE3 in cytosolic subspace
dcAMP_PDE3_cyt = (c(54) + (c(47) - 1.0) * PDE3_P_cyt) * c(42) / (1.0 + c(45) / cAMP_cyt);

camp_cAMP_cyt_pde = dcAMP_PDE2_cyt + dcAMP_PDE3_cyt + dcAMP_PDE4_cyt;
camp_cAMP_cyt_j1 = c(9) * (cAMP_cav - cAMP_cyt) / c(4);
camp_cAMP_cyt_j2 = c(10) * (cAMP_eca - cAMP_cyt) / c(4);

camp_cAMP_eca_pde = dcAMP_PDE2_eca + dcAMP_PDE4_eca;
camp_cAMP_eca_j2 = c(10) * (cAMP_eca - cAMP_cyt) / c(3);
camp_cAMP_eca_j1 = c(8) * (cAMP_cav - cAMP_eca) / c(3);

camp_cAMP_cav_pde = dcAMP_PDE2_cav + dcAMP_PDE3_cav + dcAMP_PDE4_cav;
camp_cAMP_cav_j2 = c(9) * (cAMP_cav - cAMP_cyt) / c(2);
camp_cAMP_cav_j1 = c(8) * (cAMP_cav - cAMP_eca) / c(2);
% ac
%
ac_kAC47_cyt_gsa = beta_cyt_Gs_aGTP ^ c(107);
kAC47_cyt = c(116) * (c(114) + ac_kAC47_cyt_gsa / (c(110) + ac_kAC47_cyt_gsa));
ac_kAC56_cav_gsa = beta_cav_Gs_aGTP ^ c(108);
gsi = beta_cav_Gs_aGTP ^ c(109);
kAC56_cav = c(117) * (c(115) + ac_kAC56_cav_gsa / (c(111) + ac_kAC56_cav_gsa)) * (1.0 - (1.0 - c(118) * gsi / (c(113) + gsi)) * beta_cav_Gi_bg / (c(112) + beta_cav_Gi_bg));
ac_kAC47_eca_gsa = beta_eca_Gs_aGTP ^ c(107);
kAC47_eca = c(116) * (c(114) + ac_kAC47_eca_gsa / (c(110) + ac_kAC47_eca_gsa));
ac_kAC56_cyt_gsa = beta_cyt_Gs_aGTP ^ c(108);
kAC56_cyt = c(117) * (c(115) + ac_kAC56_cyt_gsa / (c(111) + ac_kAC56_cyt_gsa));

% Rate of cAMP production by AC type 4/7 in cytoplasm
dcAMP_AC47_cyt = kAC47_cyt * c(119) * c(121);
% Rate of cAMP production by AC type 5/6 in cytoplasm
dcAMP_AC56_cyt = kAC56_cyt * c(123) * c(121);
% Rate of cAMP production by AC type 5/6 in caveolar subspace
dcAMP_AC56_cav = kAC56_cav * c(120) * c(121);
% Rate of cAMP production by AC type 4/7 in extracaveolar subspace
dcAMP_AC47_eca = kAC47_eca * c(122) * c(121);


% Caveolar concentration of cAMP
ydot(10) = 0.001 * (pka_cav_dcAMP + dcAMP_AC56_cav - camp_cAMP_cav_pde - camp_cAMP_cav_j1 - camp_cAMP_cav_j2);
% Extracaveolar concentration of cAMP
ydot(11) = 0.001 * (pka_eca_dcAMP + dcAMP_AC47_eca - camp_cAMP_eca_pde + camp_cAMP_eca_j1 - camp_cAMP_eca_j2);
% Cytosolic concentration of cAMP
ydot(12) = 0.001 * (pka_cyt_dcAMP + dcAMP_AC47_cyt + dcAMP_AC56_cyt - camp_cAMP_cyt_pde + camp_cAMP_cyt_j1 + camp_cAMP_cyt_j2);


%% Mod_PDE_Phosphorylation() function

% Concentration of phosphorylated PDE3 in the caveolar subspace
ydot(34) = 0.001 * (c(48) * pka_cav_C * (c(53) - PDE3_P_cav) - c(49) * PDE3_P_cav);
% Concentration of phosphorylated PDE3 in the cytosolic subspace
ydot(35) = 0.001 * (c(48) * pka_cyt_C * (c(54) - PDE3_P_cyt) - c(49) * PDE3_P_cyt);
% Concentration of phosphorylated PDE4 in the caveolar subspace
ydot(36) = 0.001 * (c(48) * pka_cav_C * (c(55) - PDE4_P_cav) - c(49) * PDE4_P_cav);
% Concentration of phosphorylated PDE4 in the extracaveolar subspace
ydot(37) = 0.001 * (c(48) * pka_eca_C * (c(56) - PDE4_P_eca) - c(49) * PDE4_P_eca);
% Concentration of phosphorylated PDE4 in the cytosolic subspace
ydot(38) = 0.001 * (c(48) * pka_cyt_C * (c(57) - PDE4_P_cyt) - c(49) * PDE4_P_cyt);

%% Mod_PP1_Inhibition()

pp1_PP1f_cyt_sum = c(38) - c(37) + inhib1_p;
% Concentration of uninhibited PP1 in the cytosolic compartment
PP1f_cyt = 0.5 * (sqrt(pp1_PP1f_cyt_sum ^ 2.0 + 4.0 * c(38) * c(37)) - pp1_PP1f_cyt_sum);
di = c(39) - inhib1_p;
% Concentration of phosphorylated PP1 inhibitor 1 (cytoplasmic)
ydot(39) = 0.001 * (c(31) * pka_cyt_C * di / (c(33) + di) - c(32) * c(40) * inhib1_p / (c(34) + inhib1_p));

%% Mod_Channel_Phosphorylation()

% Substrates without AKAP
% Fraction of phosphorylated PLB
ydot(42) = 0.001 * (c(132) * pka_cyt_C * (1.0 - iup_f_plb) / (c(134) + 1.0 - iup_f_plb) - c(133) * PP1f_cyt * iup_f_plb / (c(135) + iup_f_plb));
% Fraction of phosphorylated Troponin
ydot(43) = 0.001 * (c(147) * pka_cyt_C * (1.0 - f_tni) / (c(149) + 1.0 - f_tni) - c(148) * c(40) * f_tni / (c(150) + f_tni));
% Fraction of phosphorylated INa channels
ydot(44) = 0.001 * (c(130) * pka_cav_C * (1.0 - ina_f_ina) / (c(128) + 1.0 - ina_f_ina) - c(131) * c(36) * ina_f_ina / (c(129) + ina_f_ina));
% Fraction of phosphorylated INaK
ydot(45) = 0.001 * (c(124) * pka_cav_C * (1.0 - f_inak) / (c(126) + 1.0 - f_inak) - c(125) * c(36) * f_inak / (c(127) + f_inak));
% Fraction of phosphorylated IKur channels
ydot(47) = 0.001 * (c(136) * pka_eca_C * (1.0 - f_ikur) / (c(138) + 1.0 - f_ikur) - c(137) * c(35) * f_ikur / (c(139) + f_ikur));


%Substrates with AKAP
iks_sig_IKsp_dif = c(146) - IKsp;
% Concentration of phosphorylated IKs channels
ydot(41) = 0.001 * (c(140) * pka_eca_C * iks_sig_IKsp_dif / (c(142) + iks_sig_IKsp_dif) - c(141) * c(35) * IKsp / (c(143) + IKsp));
akap_sig_RyRp_dif = c(162) - RyRp;
ydot(46) = 0.001 * (c(152) * pka_cav_C * akap_sig_RyRp_dif / (c(154) + akap_sig_RyRp_dif) - c(153) * c(36) * RyRp / (c(155) + RyRp));
akap_sig_ICaLp_dif = c(164) - ICaLp;
% Concentration of phosphorylated L-type Calcium channels
ydot(40) = 0.001 * (c(157) * pka_cav_C * akap_sig_ICaLp_dif / (c(159) + akap_sig_ICaLp_dif) - c(158) * c(36) * ICaLp / (c(160) + ICaLp));

end

% Emulate the ternary operator x = (cond) ? iftrue : iffalse
function y = ifthenelse(cond, iftrue, iffalse)
if (cond)
    y = iftrue;
else
    y = iffalse;
end
end
% UNPUBLISHED PRE-REVIEW VERSION, CONFIDENTIAL.