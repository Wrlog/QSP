$PROB
Wajima 2009 coagulation network

// Humoral coagulation network with warfarin, vitamin K and heparin.
// Wajima T, Isbister GK, Duffull SB. A comprehensive model for the humoral
// coagulation network in humans. Clin Pharmacol Ther 2009;86:290-298.
//
// Converted by tools/build_coagulation.py from BioModels BIOMD0000000340
// (CC0). Species F, D, P and II are renamed Fibrin, Ddimer, Plasmin and
// FII (II is reserved in mrgsolve).
// The SBML dosing events are replaced by dose records: warfarin (mg) into
// A_warf; heparin through the heparin_infusion parameter (nM/h).
// AUC_Fibrin (added) integrates fibrin for the in-vitro clotting tests.
//
// Units: time in hours, factors in nM, warfarin in mg (A_warf) and mg/L (C_warf).

$PARAM @annotated
INVITRO                      : 0            : 1 = start from a supplied plasma state (clotting tests)
I_max                        : 1            : I_max
IC50                         : 0.34         : IC50
II0                          : 1394.4       : II(0)
VII0                         : 10           : VII(0)
IX0                          : 89.6         : IX(0)
X0                           : 174.3        : X(0)
PC0                          : 60           : PC(0)
PS0                          : 300          : PS(0)
VKH20                        : 0.1          : VKH2(0)
d_II                         : 0.01         : d_II
d_VII                        : 0.12         : d_VII
d_IX                         : 0.029        : d_IX
d_X                          : 0.018        : d_X
d_PC                         : 0.05         : d_PC
d_PS                         : 0.0165       : d_PS
vitaminK_Vc                  : 24           : vitaminK_Vc
d_VK2                        : 0.0228       : d_VK2
d_VKH2                       : 0.228        : d_VKH2
d_VKO                        : 0.228        : d_VKO
VK0                          : 1            : VK(0)
VKO0                         : 0.1          : VKO(0)
vitaminK_k21_Vc              : 0.000508333333333333 : vitaminK_k21/Vc
vitaminK_k12                 : 0.0587       : vitaminK_k12
heparin_ke                   : 0.693        : heparin_ke
warfarin_ka                  : 1            : warfarin_ka
warfarin_Vd                  : 10           : warfarin_Vd
warfarin_CL                  : 0.2          : warfarin_CL
warfarin_ke                  : 0.02         : warfarin_ke
d_XII                        : 0.012        : d_XII
d_VIII                       : 0.058        : d_VIII
d_XI                         : 0.1          : d_XI
d_V                          : 0.043        : d_V
d_Fg                         : 0.032        : d_Fg
d_XIII                       : 0.0036       : d_XIII
d_Pg                         : 0.05         : d_Pg
d_Tmod                       : 0.05         : d_Tmod
d_TFPI                       : 20           : d_TFPI
d_Pk                         : 0.05         : d_Pk
XII0                         : 375          : XII(0)
VIII0                        : 0.7          : VIII(0)
XI0                          : 30.6         : XI(0)
V0                           : 26.7         : V(0)
Fg0                          : 8945.5       : Fg(0)
XIII0                        : 70.3         : XIII(0)
Pg0                          : 2154.3       : Pg(0)
Tmod0                        : 50           : Tmod(0)
TFPI0                        : 2.5          : TFPI(0)
Pk0                          : 450          : Pk(0)
R1                           : 0.140845070422535 : R1
R2                           : 1            : R2
c44                          : 0.119718309859155 : c44
c45                          : 0.85         : c45
c46                          : 0.85         : c46
d_VK                         : 0.2052       : d_VK
warfarin_daily_dose          : 4            : warfarin_daily_dose
heparin_infusion             : 0            : heparin_infusion
heparin_bolus                : 0            : heparin_bolus
heparin_infusion_duration_hr : 24           : heparin_infusion_duration [hr]
r1__v                        : 50000        : v, reaction r1 
r1__k                        : 1E-6         : k, reaction r1 
r2__v                        : 50           : v, reaction r2 
r2__k                        : 1            : k, reaction r2 
r3__v                        : 7            : v, reaction r3 
r3__k                        : 10           : k, reaction r3 
r4__v                        : 7            : v, reaction r4 
r4__k                        : 1            : k, reaction r4 
r5__v                        : 10           : v, reaction r5 
r5__k                        : 10           : k, reaction r5 
r6__v                        : 0.1          : v, reaction r6 
r6__k                        : 10           : k, reaction r6 
r7__v                        : 0.02         : v, reaction r7 
r7__k                        : 10           : k, reaction r7 
r8__v                        : 2            : v, reaction r8 
r8__k                        : 0.1          : k, reaction r8 
r9__v                        : 1E-9         : v, reaction r9 
r9__k                        : 10           : k, reaction r9 
r10_v                        : 50000        : v, reaction r10
r10_k                        : 10           : k, reaction r10
r11_v                        : 50           : v, reaction r11
r11_k                        : 1            : k, reaction r11
r12_v                        : 100          : v, reaction r12
r12_k                        : 10           : k, reaction r12
r13_v                        : 9            : v, reaction r13
r13_k                        : 500          : k, reaction r13
r14_v                        : 20000        : v, reaction r14
r14_k                        : 0.5          : k, reaction r14
r15_v                        : 500          : v, reaction r15
r15_k                        : 500          : k, reaction r15
r16_v                        : 7            : v, reaction r16
r16_k                        : 10           : k, reaction r16
r17_v                        : 7            : v, reaction r17
r17_k                        : 10           : k, reaction r17
r18_v                        : 7            : v, reaction r18
r18_k                        : 100          : k, reaction r18
r19_v                        : 1            : v, reaction r19
r19_k                        : 1            : k, reaction r19
r20_v                        : 7            : v, reaction r20
r20_k                        : 1            : k, reaction r20
r21_v                        : 7            : v, reaction r21
r21_k                        : 5000         : k, reaction r21
r22_v                        : 5            : v, reaction r22
r22_k                        : 10000        : k, reaction r22
r23_v                        : 2            : v, reaction r23
r23_k                        : 1            : k, reaction r23
r24_v                        : 7            : v, reaction r24
r24_k                        : 1            : k, reaction r24
r25_v                        : 2            : v, reaction r25
r25_k                        : 1            : k, reaction r25
r26_c                        : 0.01         : c, reaction r26
r27_c                        : 0.5          : c, reaction r27
r28_c                        : 0.5          : c, reaction r28
r29_c                        : 0.5          : c, reaction r29
r30_c                        : 0.1          : c, reaction r30
r31_c                        : 0.5          : c, reaction r31
r32_c                        : 0.5          : c, reaction r32
r33_v                        : 70           : v, reaction r33
r33_k                        : 1            : k, reaction r33
r34_v                        : 900          : v, reaction r34
r34_k                        : 200          : k, reaction r34
r35_v                        : 70           : v, reaction r35
r35_k                        : 1            : k, reaction r35
r36_v                        : 1000         : v, reaction r36
r36_k                        : 1            : k, reaction r36
r37_c                        : 0.5          : c, reaction r37
r38_v                        : 1            : v, reaction r38
r38_k                        : 10           : k, reaction r38
r39_v                        : 1            : v, reaction r39
r39_k                        : 10           : k, reaction r39
r40_v                        : 0.2          : v, reaction r40
r40_k                        : 10           : k, reaction r40
r41_v                        : 7            : v, reaction r41
r41_k                        : 1            : k, reaction r41
r42_v                        : 70           : v, reaction r42
r42_k                        : 1            : k, reaction r42
r43_v                        : 7            : v, reaction r43
r43_k                        : 1            : k, reaction r43
dF_k1                        : 0.05         : k1, reaction dF
dXF_k1                       : 0.05         : k1, reaction dXF
dIIa_k1                      : 67.4         : k1, reaction dIIa
dTF_k1                       : 0.05         : k1, reaction dTF
dVa_k1                       : 20           : k1, reaction dVa
dVIIa_k1                     : 20           : k1, reaction dVIIa
dVIIIa_k1                    : 20           : k1, reaction dVIIIa
dXa_k1                       : 20           : k1, reaction dXa
dIXa_k1                      : 20           : k1, reaction dIXa
dXIIa_k1                     : 20           : k1, reaction dXIIa
dXIIIa_k1                    : 0.69         : k1, reaction dXIIIa
dK_k1                        : 20           : k1, reaction dK
dP_k1                        : 20           : k1, reaction dP
dAPC_k1                      : 20.4         : k1, reaction dAPC
dFDP_k1                      : 3.5          : k1, reaction dFDP
dD_k1                        : 0.1          : k1, reaction dD
dVIIa_TF_k1                  : 20           : k1, reaction dVIIa_TF
dVII_TF_k1                   : 0.7          : k1, reaction dVII_TF
dAPC_PS_k1                   : 20           : k1, reaction dAPC_PS
dVa_Xa_k1                    : 20           : k1, reaction dVa_Xa
dIXa_VIIIa_k1                : 20           : k1, reaction dIXa_VIIIa
dIIa_Tmod_k1                 : 20           : k1, reaction dIIa_Tmod
dXa_TFPI_k1                  : 20           : k1, reaction dXa_TFPI
dVIIa_TF_Xa_TFPI_k1          : 20           : k1, reaction dVIIa_TF_Xa_TFPI
dTAT_k1                      : 0.2          : k1, reaction dTAT
dCA_k1                       : 0.05         : k1, reaction dCA
dXIa_k1                      : 20           : k1, reaction dXIa
dXI_k1                       : 0.1          : k1, reaction dXI

$CMT @annotated
IIa              : IIa (nM)
VIII             : VIII (nM)
VIIIa            : VIIIa (nM)
APC_PS           : APC_PS (nM)
IX               : IX (nM)
IXa              : IXa (nM)
XIa              : XIa (nM)
XI               : XI (nM)
XIIa             : XIIa (nM)
VII              : VII (nM)
VIIa             : VIIa (nM)
X                : X (nM)
Xa               : Xa (nM)
IXa_VIIIa        : IXa_VIIIa (nM)
V                : V (nM)
Va               : Va (nM)
FII              : II (nM)
Fibrin           : F (nM)
Fg               : Fg (nM)
Plasmin          : P (nM)
XF               : XF (nM)
XIII             : XIII (nM)
Pg               : Pg (nM)
APC              : APC (nM)
IIa_Tmod         : IIa_Tmod (nM)
PC               : PC (nM)
Tmod             : Tmod (nM)
TF               : TF (nM)
VIIa_TF          : VIIa_TF (nM)
VII_TF           : VII_TF (nM)
Xa_TFPI          : Xa_TFPI (nM)
TFPI             : TFPI (nM)
PS               : PS (nM)
VKH2             : VKH2 (nM)
Va_Xa            : Va_Xa (nM)
CA               : CA (nM)
XII              : XII (nM)
K                : K (nM)
ATIII_Heparin    : ATIII_Heparin (nM)
Xa_ATIII_Heparin : Xa_ATIII_Heparin (nM)
VK               : VK (nM)
C_warf           : Warfarin plasma concentration (mg/L)
VKO              : VKO (nM)
Pk               : Pk (nM)
FDP              : FDP (nM)
Ddimer           : D (nM)
TAT              : TAT (nM)
VIIa_TF_Xa_TFPI  : VIIa_TF_Xa_TFPI (nM)
XIIIa            : XIIIa (nM)
IIa_ATIII_Heparin : IIa_ATIII_Heparin (nM)
IXa_ATIII_Heparin : IXa_ATIII_Heparin (nM)
A_warf           : Warfarin in the gut (mg)
VK_p             : VK_p (nM)
AUC_Fibrin       : Integral of fibrin, for the clotting-time tests (nM*h)

$MAIN
if (INVITRO < 0.5) {
  IIa_0 = 0;
  VIII_0 = VIII0;
  VIIIa_0 = 0;
  APC_PS_0 = 0;
  IX_0 = IX0;
  IXa_0 = 0;
  XIa_0 = 0;
  XI_0 = XI0;
  XIIa_0 = 0;
  VII_0 = VII0;
  VIIa_0 = 0;
  X_0 = X0;
  Xa_0 = 0;
  IXa_VIIIa_0 = 0;
  V_0 = V0;
  Va_0 = 0;
  FII_0 = II0;
  Fibrin_0 = 0;
  Fg_0 = Fg0;
  Plasmin_0 = 0;
  XF_0 = 0;
  XIII_0 = XIII0;
  Pg_0 = Pg0;
  APC_0 = 0;
  IIa_Tmod_0 = 0;
  PC_0 = PC0;
  Tmod_0 = Tmod0;
  TF_0 = 0;
  VIIa_TF_0 = 0;
  VII_TF_0 = 0;
  Xa_TFPI_0 = 0;
  TFPI_0 = TFPI0;
  PS_0 = PS0;
  VKH2_0 = VKH20;
  Va_Xa_0 = 0;
  CA_0 = 0;
  XII_0 = XII0;
  K_0 = 0;
  ATIII_Heparin_0 = heparin_bolus;
  Xa_ATIII_Heparin_0 = 0;
  VK_0 = VK0;
  C_warf_0 = 0;
  VKO_0 = VKO0;
  Pk_0 = Pk0;
  FDP_0 = 0;
  Ddimer_0 = 0;
  TAT_0 = 0;
  VIIa_TF_Xa_TFPI_0 = 0;
  XIIIa_0 = 0;
  IIa_ATIII_Heparin_0 = 0;
  IXa_ATIII_Heparin_0 = 0;
  A_warf_0 = 0;
  VK_p_0 = ((VK0 * vitaminK_k12) / vitaminK_k21_Vc);
}

$ODE
double DP = (FDP + Ddimer);
double v_r1 = (1.0 * (((r1__v * VIII * IIa) / (r1__k + IIa))));
double v_r2 = (1.0 * (((r2__v * VIIIa * APC_PS) / (r2__k + APC_PS))));
double v_r3 = (1.0 * (((r3__v * IX * XIa) / (r3__k + XIa))));
double v_r4 = (1.0 * (((r4__v * XI * XIIa) / (r4__k + XIIa))));
double v_r5 = (1.0 * (((r5__v * XI * IIa) / (r5__k + IIa))));
double v_r6 = (1.0 * (((r6__v * VII * IIa) / (r6__k + IIa))));
double v_r7 = (1.0 * (((r7__v * X * IXa) / (r7__k + IXa))));
double v_r8 = (1.0 * (((r8__v * X * IXa_VIIIa) / (r8__k + IXa_VIIIa))));
double v_r9 = (1.0 * (((r9__v * X * VIIa) / (r9__k + VIIa))));
double v_r10 = (1.0 * (((r10_v * V * IIa) / (r10_k + IIa))));
double v_r11 = (1.0 * (((r11_v * Va * APC_PS) / (r11_k + APC_PS))));
double v_r12 = (1.0 * (((r12_v * FII * Va_Xa) / (r12_k + Va_Xa))));
double v_r13 = (1.0 * (((r13_v * FII * Xa) / (r13_k + Xa))));
double v_r14 = (1.0 * (((r14_v * Fg * IIa) / (r14_k + IIa))));
double v_r15 = (1.0 * (((r15_v * Fg * Plasmin) / (r15_k + Plasmin))));
double v_r16 = (1.0 * (((r16_v * Fibrin * XIIIa) / (r16_k + XIIIa))));
double v_r17 = (1.0 * (((r17_v * Fibrin * Plasmin) / (r17_k + Plasmin))));
double v_r18 = (1.0 * (((r18_v * XF * Plasmin) / (r18_k + Plasmin))));
double v_r19 = (1.0 * (((r19_v * XF * APC_PS) / (r19_k + APC_PS))));
double v_r20 = (1.0 * (((r20_v * XIII * IIa) / (r20_k + IIa))));
double v_r21 = (1.0 * (((r21_v * Pg * IIa) / (r21_k + IIa))));
double v_r22 = (1.0 * (((r22_v * Pg * Fibrin) / (r22_k + Fibrin))));
double v_r23 = (1.0 * (((r23_v * Pg * APC_PS) / (r23_k + APC_PS))));
double v_r24 = (1.0 * (((r24_v * PC * IIa_Tmod) / (r24_k + IIa_Tmod))));
double v_r25 = (1.0 * (((r25_v * Va_Xa * APC_PS) / (r25_k + APC_PS))));
double v_r26 = (1.0 * (((IXa * VIIIa) / r26_c)));
double v_r27 = (1.0 * (((Va * Xa) / r27_c)));
double v_r28 = (1.0 * (((IIa * Tmod) / r28_c)));
double v_r29 = (1.0 * (((VIIa * TF) / r29_c)));
double v_r30 = (1.0 * (((VII * TF) / r30_c)));
double v_r31 = (1.0 * (((VIIa_TF * Xa_TFPI) / r31_c)));
double v_r32 = (1.0 * (((Xa * TFPI) / r32_c)));
double v_r33 = (1.0 * (((r33_v * VII_TF * Xa) / (r33_k + Xa))));
double v_r34 = (1.0 * (((r34_v * X * VIIa_TF) / (r34_k + VIIa_TF))));
double v_r35 = (1.0 * (((r35_v * IX * VIIa_TF) / (r35_k + VIIa_TF))));
double v_r36 = (1.0 * (((r36_v * VII_TF * TF) / (r36_k + TF))));
double v_r37 = (1.0 * (((APC * PS) / r37_c)));
double v_r38 = (1.0 * (((r38_v * VII * Xa) / (r38_k + Xa))));
double v_r39 = (1.0 * (((r39_v * VII * VIIa_TF) / (r39_k + VIIa_TF))));
double v_r40 = (1.0 * (((r40_v * VII * IXa) / (r40_k + IXa))));
double v_r41 = (1.0 * (((r41_v * XII * CA) / (r41_k + CA))));
double v_r42 = (1.0 * (((r42_v * XII * K) / (r42_k + K))));
double v_r43 = (1.0 * (((r43_v * Pk * XIIa) / (r43_k + XIIa))));
double v_r44 = (1.0 * (((IIa * ATIII_Heparin) / c44)));
double v_r45 = (1.0 * (((Xa * ATIII_Heparin) / c45)));
double v_r46 = (1.0 * (((IXa * ATIII_Heparin) / c46)));
double v_r47 = (1.0 * ((d_VK2 * VK * (1.0 - ((I_max * C_warf) / (IC50 + C_warf))))));
double v_r48 = (1.0 * ((d_VKO * VKO * (1.0 - ((I_max * C_warf) / (IC50 + C_warf))))));
double v_pII_VKH2 = (1.0 * (((d_II * II0 * VKH2) / VKH20)));
double v_pVII_VKH2 = (1.0 * (((d_VII * VII0 * VKH2) / VKH20)));
double v_pIX_VKH2 = (1.0 * (((d_IX * IX0 * VKH2) / VKH20)));
double v_pX_VKH2 = (1.0 * (((d_X * X0 * VKH2) / VKH20)));
double v_pPC_VKH2 = (1.0 * (((d_PC * PC0 * VKH2) / VKH20)));
double v_pPS_VKH2 = (1.0 * (((d_PS * PS0 * VKH2) / VKH20)));
double v_dFg = (1.0 * d_Fg * Fg);
double v_dF = (1.0 * dF_k1 * Fibrin);
double v_dXF = (1.0 * dXF_k1 * XF);
double v_dII = (1.0 * d_II * FII);
double v_dIIa = (1.0 * dIIa_k1 * IIa);
double v_dTF = (1.0 * dTF_k1 * TF);
double v_dV = (1.0 * d_V * V);
double v_dVa = (1.0 * dVa_k1 * Va);
double v_dVII = (1.0 * d_VII * VII);
double v_dVIIa = (1.0 * dVIIa_k1 * VIIa);
double v_dVIII = (1.0 * d_VIII * VIII);
double v_dVIIIa = (1.0 * dVIIIa_k1 * VIIIa);
double v_dX = (1.0 * d_X * X);
double v_dXa = (1.0 * dXa_k1 * Xa);
double v_dIX = (1.0 * d_IX * IX);
double v_dIXa = (1.0 * dIXa_k1 * IXa);
double v_dXII = (1.0 * d_XII * XII);
double v_dXIIa = (1.0 * dXIIa_k1 * XIIa);
double v_dXIII = (1.0 * d_XIII * XIII);
double v_dXIIIa = (1.0 * dXIIIa_k1 * XIIIa);
double v_dPk = (1.0 * d_Pk * Pk);
double v_dK = (1.0 * dK_k1 * K);
double v_dPg = (1.0 * d_Pg * Pg);
double v_dP = (1.0 * dP_k1 * Plasmin);
double v_dPC = (1.0 * d_PC * PC);
double v_dAPC = (1.0 * dAPC_k1 * APC);
double v_dPS = (1.0 * d_PS * PS);
double v_dFDP = (1.0 * dFDP_k1 * FDP);
double v_dD = (1.0 * dD_k1 * Ddimer);
double v_dTFPI = (1.0 * d_TFPI * TFPI);
double v_dVIIa_TF = (1.0 * dVIIa_TF_k1 * VIIa_TF);
double v_dVII_TF = (1.0 * dVII_TF_k1 * VII_TF);
double v_dAPC_PS = (1.0 * dAPC_PS_k1 * APC_PS);
double v_dVa_Xa = (1.0 * dVa_Xa_k1 * Va_Xa);
double v_dIXa_VIIIa = (1.0 * dIXa_VIIIa_k1 * IXa_VIIIa);
double v_dTmod = (1.0 * d_Tmod * Tmod);
double v_dIIa_Tmod = (1.0 * dIIa_Tmod_k1 * IIa_Tmod);
double v_dXa_TFPI = (1.0 * dXa_TFPI_k1 * Xa_TFPI);
double v_dVIIa_TF_Xa_TFPI = (1.0 * dVIIa_TF_Xa_TFPI_k1 * VIIa_TF_Xa_TFPI);
double v_dTAT = (1.0 * dTAT_k1 * TAT);
double v_dCA = (1.0 * dCA_k1 * CA);
double v_dXIa = (1.0 * dXIa_k1 * XIa);
double v_dVKH2 = (1.0 * d_VKH2 * VKH2);
double v_VK_transport = (1.0 * ((vitaminK_k12 * VK) - (vitaminK_k21_Vc * VK_p)));
double v_eHeparin = (1.0 * heparin_ke * ATIII_Heparin);
double v_eHeparinXa = (1.0 * heparin_ke * Xa_ATIII_Heparin);
double v_eHeparinIXa = (1.0 * heparin_ke * IXa_ATIII_Heparin);
double v_eHeparinIIa = (1.0 * heparin_ke * IIa_ATIII_Heparin);
double v_dXI = (1.0 * dXI_k1 * XI);
double v_pXII = (1.0 * ((XII0 * d_XII)));
double v_pVIII = (1.0 * ((VIII0 * d_VIII)));
double v_pXI = (1.0 * ((XI0 * d_XI)));
double v_pV = (1.0 * ((V0 * d_V)));
double v_pFg = (1.0 * ((Fg0 * d_Fg)));
double v_pXIII = (1.0 * ((XIII0 * d_XIII)));
double v_pPg = (1.0 * ((Pg0 * d_Pg)));
double v_pTmod = (1.0 * ((Tmod0 * d_Tmod)));
double v_pTFPI = (1.0 * ((TFPI0 * d_TFPI)));
double v_pPk = (1.0 * ((Pk0 * d_Pk)));
double v_pVK = (1.0 * ((VK0 * d_VK)));
double v_dVK = (1.0 * d_VK * VK);
double v_pHeparin = (1.0 * (heparin_infusion));

dxdt_IIa = v_r12 + v_r13 - v_r28 - v_r44 - v_dIIa;
dxdt_VIII = -v_r1 - v_dVIII + v_pVIII;
dxdt_VIIIa = v_r1 - v_r2 - v_r26 - v_dVIIIa;
dxdt_APC_PS = v_r37 - v_dAPC_PS;
dxdt_IX = -v_r3 - v_r35 + v_pIX_VKH2 - v_dIX;
dxdt_IXa = v_r3 - v_r26 + v_r35 - v_r46 - v_dIXa;
dxdt_XIa = v_r4 + v_r5 - v_dXIa;
dxdt_XI = -v_r4 - v_r5 - v_dXI + v_pXI;
dxdt_XIIa = v_r41 + v_r42 - v_dXIIa;
dxdt_VII = -v_r6 - v_r30 - v_r38 - v_r39 - v_r40 + v_pVII_VKH2 - v_dVII;
dxdt_VIIa = v_r6 - v_r29 + v_r38 + v_r39 + v_r40 - v_dVIIa;
dxdt_X = -v_r7 - v_r8 - v_r9 - v_r34 + v_pX_VKH2 - v_dX;
dxdt_Xa = v_r7 + v_r8 + v_r9 - v_r27 - v_r32 + v_r34 - v_r45 - v_dXa;
dxdt_IXa_VIIIa = v_r26 - v_dIXa_VIIIa;
dxdt_V = -v_r10 - v_dV + v_pV;
dxdt_Va = v_r10 - v_r11 - v_r27 - v_dVa;
dxdt_FII = -v_r12 - v_r13 + v_pII_VKH2 - v_dII;
dxdt_Fibrin = v_r14 - v_r16 - v_r17 - v_dF;
dxdt_Fg = -v_r14 - v_r15 - v_dFg + v_pFg;
dxdt_Plasmin = v_r21 + v_r22 + v_r23 - v_dP;
dxdt_XF = v_r16 - v_r18 - v_r19 - v_dXF;
dxdt_XIII = -v_r20 - v_dXIII + v_pXIII;
dxdt_Pg = -v_r21 - v_r22 - v_r23 - v_dPg + v_pPg;
dxdt_APC = v_r24 - v_r37 - v_dAPC;
dxdt_IIa_Tmod = v_r28 - v_dIIa_Tmod;
dxdt_PC = -v_r24 + v_pPC_VKH2 - v_dPC;
dxdt_Tmod = -v_r28 - v_dTmod + v_pTmod;
dxdt_TF = -v_r29 - v_r30 - v_dTF;
dxdt_VIIa_TF = v_r29 - v_r31 + v_r33 + v_r36 - v_dVIIa_TF;
dxdt_VII_TF = v_r30 - v_r33 - v_r36 - v_dVII_TF;
dxdt_Xa_TFPI = -v_r31 + v_r32 - v_dXa_TFPI;
dxdt_TFPI = -v_r32 - v_dTFPI + v_pTFPI;
dxdt_PS = -v_r37 + v_pPS_VKH2 - v_dPS;
dxdt_VKH2 = v_r47 - v_dVKH2;
dxdt_Va_Xa = -v_r25 + v_r27 - v_dVa_Xa;
dxdt_CA = -v_dCA;
dxdt_XII = -v_r41 - v_r42 - v_dXII + v_pXII;
dxdt_K = v_r43 - v_dK;
dxdt_ATIII_Heparin = -v_r44 - v_r45 - v_r46 - v_eHeparin + v_pHeparin;
dxdt_Xa_ATIII_Heparin = v_r45 - v_eHeparinXa;
dxdt_VK = -v_r47 + v_r48 - v_VK_transport + v_pVK - v_dVK;
dxdt_VKO = -v_r48 + v_dVKH2;
dxdt_Pk = -v_r43 - v_dPk + v_pPk;
dxdt_FDP = v_r15 + v_r17 + v_dFg + v_dF - v_dFDP;
dxdt_Ddimer = v_r18 + v_r19 + v_dXF - v_dD;
dxdt_TAT = v_dIIa - v_dTAT;
dxdt_VIIa_TF_Xa_TFPI = v_r31 - v_dVIIa_TF_Xa_TFPI;
dxdt_XIIIa = v_r20 - v_dXIIIa;
dxdt_IIa_ATIII_Heparin = v_r44 - v_eHeparinIIa;
dxdt_IXa_ATIII_Heparin = v_r46 - v_eHeparinIXa;
dxdt_VK_p = v_VK_transport;
dxdt_C_warf = (((warfarin_ka * A_warf) / warfarin_Vd) - (warfarin_ke * C_warf));
dxdt_A_warf = ((-warfarin_ka) * A_warf);
dxdt_AUC_Fibrin = Fibrin;
