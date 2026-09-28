$PROB
4GI model: glucose, insulin, GLP-1, glucagon and GIP, with liraglutide

// Bosch R, Petrone M, Arends R, Vicini P, Sijbrands EJG, Hoefman S,
// Snelder N. A novel integrated QSP model of in vivo human glucose
// regulation to support the development of a glucagon/GLP-1 dual agonist.
// CPT Pharmacometrics Syst Pharmacol 2022;11:302-317 (open access).
//
// Transcribed from the final NONMEM control stream in the article's
// Supplementary Material 2 ($PK and $DES), with the final estimates as
// parameter defaults. Study-specific branches (meal bioavailabilities,
// infusion durations, per-subject baselines) are not carried over: meals
// enter as glucose amounts already multiplied by their bioavailability,
// and baselines are parameters.
//
// Liraglutide PK: Watson et al. 2010 (one compartment, first-order SC
// absorption), free fraction 0.5%, GLP-1R EC50 6 pM (as in the paper).
//
// PAT = 1 selects the healthy-volunteer parameter set, PAT = 0 type 2
// diabetes, as in the control stream.
//
// Units: time in hours; glucose mmol and mM; hormones pmol and pM;
// liraglutide pmol (1 mg = 1e9 / 3751.2 pmol).

$PARAM @annotated
PAT       : 0         : 1 = healthy, 0 = type 2 diabetes
BW        : 90        : Body weight (kg), for liraglutide PK
BSLglc    : 7.8       : Baseline glucose (mM)
BSLins    : 52        : Baseline insulin (pM)
BSLglp    : 4.61      : Baseline GLP-1 (pM)
BSLglg    : 28.2      : Baseline glucagon (pM)
BSLgip    : 10        : Baseline GIP (pM)

// glucose
CLglc_T2  : 1.72      : Insulin-independent glucose clearance, T2DM (L/h)
CLglci_T2 : 0.0256    : Insulin-dependent glucose clearance, T2DM ((L/h)/pM)
CLglc_HV  : 5.36      : Insulin-independent glucose clearance, healthy (L/h)
CLglci_HV : 0.072     : Insulin-dependent glucose clearance, healthy ((L/h)/pM)
Qglc      : 26.5      : Intercompartmental clearance, glucose (L/h)
VCglc     : 9.33      : Central volume, glucose (L)
VPglc     : 8.56      : Peripheral volume, glucose (L)
KAglc     : 0.853283  : Glucose absorption rate (1/h)
Keglc     : 0.281284  : Glucose buffer to gut transit rate (1/h)
Kelglc    : 1.93412   : Gut transit rate (1/h)

// insulin
CLins     : 73.2      : Insulin clearance (L/h)
VCins     : 6.09      : Insulin central volume (L)
lKE0ins   : -0.158729 : log insulin effect-compartment rate (1/h)

// GLP-1
VCglp     : 16.0441   : GLP-1 central volume (L)
lVM       : 7.96952   : log Vmax, GLP-1 elimination (pmol/h)
lKM       : 4.90603   : log Km, GLP-1 elimination (pM)

// glucagon
CLglg     : 453.174   : Glucagon clearance (L/h)
VCglg     : 64.6416   : Glucagon central volume (L)

// GIP
CLgip     : 86.8      : GIP clearance (L/h)
VCgip     : 9.21      : GIP central volume (L)
Qgip      : 49.4      : GIP intercompartmental clearance (L/h)
VPgip     : 22.8      : GIP peripheral volume (L)

// food and hormone interactions
FDGLP     : 0.0101601 : Food-driven GLP-1 production (1/mmol)
FDGLP_2   : 0         : Glucose-driven GLP-1 production (1/mmol)
POW_1     : 1         : Power, glucose-driven GLP-1 production
FDGIP     : 0.0342994 : Food-driven GIP production (1/mmol)
FDGLG     : 0.00328611 : Food-driven glucagon production (1/mmol)
FDINS     : 0         : Food-driven first-phase insulin
GLCINS_S  : 2.4636    : Glucose-dependent insulin secretion (power)
POW_2H    : 0.924585  : Glucose feedback on glucagon (power)
lEMAX_1   : 2.36968   : Emax, GLP-1 on insulin secretion (on exp scale: EMAX_1 = exp(lEMAX_1))
lEC50_1   : 3.28677   : log EC50, GLP-1 on insulin secretion (pM)
HILL_1    : 1.79297   : Hill, GLP-1 on insulin secretion
EMAX_2    : 1         : Emax, GLP-1 slowing glucose absorption
lEC50_2   : 4.96625   : log EC50, GLP-1 on absorption (pM)
HILL_2    : 1         : Hill, GLP-1 on absorption
EMAX_3    : 1         : Emax, GLP-1 inhibition of glucagon
lEC50_3   : 4.60247   : log EC50, GLP-1 on glucagon (pM)
HILL_3    : 1         : Hill, GLP-1 on glucagon
EMAX_4    : 6.72498   : Emax, glucagon on glucose production
lEC50_4   : 4.59      : log EC50, glucagon on glucose production (pM)
HILL_4    : 1         : Hill, glucagon on glucose production
GIPINS    : 1         : GIP effect on insulin secretion
POW_3_HV  : 0.285572  : Power, GIP on insulin secretion (healthy; 0 in T2DM)
POW_4     : 0.109068  : Power, GIP on glucagon

// liraglutide (Watson 2010) and potency
KAdrug    : 0.154     : Liraglutide absorption rate (1/h)
CLdrug_kg : 0.013     : Liraglutide clearance (L/h/kg)
VCdrug_kg : 0.16      : Liraglutide volume (L/kg)
fu        : 0.005     : Liraglutide free fraction
EC50lira  : 6         : Liraglutide GLP-1R EC50 (pM)
ECGLP     : 1.919     : Endogenous GLP-1 in vitro EC50 (pM)

$CMT @annotated
Dglc    : Glucose dose in the stomach (mmol)
CTglc   : Central glucose (mmol)
CTins   : Central insulin (pmol)
CTglp   : Central GLP-1 (pmol)
CTglg   : Central glucagon (pmol)
CTgip   : Central GIP (pmol)
Pglc    : Peripheral glucose (mmol)
AUCG    : Glucose AUC (mM*h)
Bglc    : Glucose buffer (mmol)
Eins    : Insulin effect compartment (pM)
Pgip    : Peripheral GIP (pmol)
GLCgut1 : Gut transit 1 (mmol)
GLCgut2 : Gut transit 2 (mmol)
GLCgut3 : Gut transit 3 (mmol)
Ddrug   : Liraglutide SC depot (pmol)
Cdrugc  : Liraglutide central (pmol)

$MAIN
double KE0ins = exp(lKE0ins);
double VM = exp(lVM);
double KM = exp(lKM);
double EMAX_1 = exp(lEMAX_1);
double EC50_1 = exp(lEC50_1);
double EC50_2 = exp(lEC50_2);
double EC50_3 = exp(lEC50_3);
double EC50_4 = exp(lEC50_4);
double CLdrug = CLdrug_kg * BW;
double VCdrug = VCdrug_kg * BW;
double ECGLP1 = EC50lira / ECGLP * EC50_1;
double ECGLP2 = EC50lira / ECGLP * EC50_2;
double ECGLP3 = EC50lira / ECGLP * EC50_3;

// healthy volunteers (PAT = 1) vs type 2 diabetes (PAT = 0)
double CLglc = PAT > 0.5 ? CLglc_HV : CLglc_T2;
double CLglci = PAT > 0.5 ? CLglci_HV : CLglci_T2;
double POW_3 = PAT > 0.5 ? POW_3_HV : 0.0;

CTglc_0 = BSLglc * VCglc;
CTins_0 = BSLins * VCins;
CTglp_0 = BSLglp * VCglp;
CTglg_0 = BSLglg * VCglg;
CTgip_0 = BSLgip * VCgip;
Pglc_0 = BSLglc * VPglc;
Eins_0 = BSLins;
Pgip_0 = BSLgip * VPgip;

$ODE
double Cglc = CTglc / VCglc;
double Cins = CTins / VCins;
double Cglp = CTglp / VCglp;
double Cglg = CTglg / VCglg;
double Cgip = CTgip / VCgip;
double Cdrug = Cdrugc / VCdrug;
double Cdrugf = Cdrug * fu;

// glucose feedback on glucagon; off below baseline in T2DM
double POW_2 = POW_2H;
if (PAT < 0.5 && Cglc > 0 && Cglc < BSLglc) POW_2 = 0;
double glcEFFglg = pow(BSLglc / Cglc, POW_2);

// GLP-1 (and liraglutide) on glucose-dependent insulin secretion
double GLPINS_S0 = EMAX_1 * (pow(BSLglp / EC50_1, HILL_1) / (1 + pow(BSLglp / EC50_1, HILL_1)));
double NUM1 = pow(Cglp / EC50_1, HILL_1) + pow(Cdrugf / ECGLP1, HILL_1);
double GLPINS_S = EMAX_1 * NUM1 / (1 + NUM1);
double glpEFFins = (1 + GLPINS_S) / (1 + GLPINS_S0);

// GLP-1 (and liraglutide) slowing glucose absorption
double NUM2 = pow(Cglp / EC50_2, HILL_2) + pow(Cdrugf / ECGLP2, HILL_2);
double GLPGLU_AI = EMAX_2 * NUM2 / (1 + NUM2);

// GLP-1 (and liraglutide) inhibiting glucagon secretion
double GLPGLG_I0 = EMAX_3 * (pow(BSLglp / EC50_3, HILL_3) / (1 + pow(BSLglp / EC50_3, HILL_3)));
double NUM3 = pow(Cglp / EC50_3, HILL_3) + pow(Cdrugf / ECGLP3, HILL_3);
double GLPGLG_I = EMAX_3 * NUM3 / (1 + NUM3);
double glpEFFglg = (1 - GLPGLG_I) / (1 - GLPGLG_I0);

// glucagon driving glucose production
double GLGGLC_S0 = EMAX_4 * (pow(BSLglg / EC50_4, HILL_4) / (1 + pow(BSLglg / EC50_4, HILL_4)));
double GLGGLC_S = EMAX_4 * (pow(Cglg / EC50_4, HILL_4) / (1 + pow(Cglg / EC50_4, HILL_4)));
double glgEFFglc = (1 + GLGGLC_S) / (1 + GLGGLC_S0);

// GIP on insulin and glucagon secretion
double GIPINS_S = GIPINS * pow(fmax(Cgip, 0.0), POW_3);
double GIPINS_S0 = GIPINS * pow(BSLgip, POW_3);
double gipEFFglg = CTgip > 0 ? pow(Cgip / BSLgip, POW_4) : 1.0;

double STglc = GLPINS_S + GIPINS_S;
double STglc0 = GLPINS_S0 + GIPINS_S0;

// food effects, driven by glucose in the buffer and the gut
double FDGLP_S = Bglc > 0 ? FDGLP * Bglc : 0.0;
double FDGLP_S2 = (GLCgut3 > 0 && PAT > 0.5) ? FDGLP_2 * pow(GLCgut3, POW_1) : 0.0;
double FDGIP_S = Bglc > 0 ? FDGIP * Bglc : 0.0;
double FDGLG_S = Bglc > 0 ? FDGLG * Bglc : 0.0;
double FDINS_S = (Bglc > 0 && PAT > 0.5) ? FDINS * Bglc : 0.0;

// baseline production rates
double KINglc = BSLglc * (CLglc + CLglci * BSLins);
double KINins = BSLins * CLins / (1 + STglc0 * pow(BSLglc, GLCINS_S));
double KINglp = VM * (BSLglp * VCglp) / (KM + BSLglp);
double KINglg = BSLglg * CLglg;
double KINgip = BSLgip * CLgip;

double KAglc2 = KAglc * (1 - GLPGLU_AI);
double K27 = Qglc / VCglc;
double K72 = Qglc / VPglc;
double K612 = Qgip / VCgip;
double K126 = Qgip / VPgip;

dxdt_Dglc = -KAglc2 * Dglc;
dxdt_CTglc = KAglc * Bglc + KINglc * glgEFFglc - K27 * CTglc + K72 * Pglc -
             (CLglc / VCglc) * CTglc - (CLglci * Eins / VCglc) * CTglc;
dxdt_CTins = KINins * (1 + FDINS_S) * (1 + STglc * pow(fmax(Cglc, 0.0), GLCINS_S)) - (CLins / VCins) * CTins;
dxdt_CTglp = KINglp * (1 + FDGLP_S) * (1 + FDGLP_S2) - VM * CTglp / (KM + Cglp);
dxdt_CTglg = KINglg * (1 + FDGLG_S) * glcEFFglg * glpEFFglg * gipEFFglg - (CLglg / VCglg) * CTglg;
dxdt_CTgip = KINgip * (1 + FDGIP_S) - (CLgip / VCgip) * CTgip - K612 * CTgip + K126 * Pgip;
dxdt_Pglc = K27 * CTglc - K72 * Pglc;
dxdt_AUCG = CTglc / VCglc;
dxdt_Bglc = KAglc2 * Dglc - KAglc * Bglc - Keglc * Bglc;
dxdt_Eins = KE0ins * (CTins / VCins - Eins);
dxdt_Pgip = K612 * CTgip - K126 * Pgip;
dxdt_GLCgut1 = Keglc * Bglc - Kelglc * GLCgut1;
dxdt_GLCgut2 = Kelglc * (GLCgut1 - GLCgut2);
dxdt_GLCgut3 = Kelglc * (GLCgut2 - GLCgut3);
dxdt_Ddrug = -KAdrug * Ddrug;
dxdt_Cdrugc = KAdrug * Ddrug - (CLdrug / VCdrug) * Cdrugc;

$TABLE
double GLC = CTglc / VCglc;
double INS = CTins / VCins;
double GLP1 = CTglp / VCglp;
double GLG = CTglg / VCglg;
double GIP = CTgip / VCgip;
double LIRA = Cdrugc / VCdrug;

$CAPTURE GLC INS GLP1 GLG GIP LIRA
