$PROB
// Monoclonal antibody with target-mediated drug disposition,
// quasi-steady-state approximation (Gibiansky et al. 2008), plus a
// biomarker whose production is driven by free target.
//
//   - two-compartment antibody PK, IV (dose into CENT) or SC (dose into DEPOT)
//   - membrane target synthesised at KDEG * R0, degraded at KDEG when free
//   - drug-target complex internalised at KINT
//   - free drug from total drug and total target by the QSS quadratic
//   - biomarker: dB/dt = KOUT * B0 * (1 - IMAX + IMAX * (Rfree / R0)^GAMMA) - KOUT * B
//
// Individual parameters (allometry, variability) are computed in R
// (R/population.R) and supplied per subject, so this file and the base-R
// engine in R/tmdd_engine.R are fed identical inputs.
//
// Units: time in days, drug in nmol / nM, target in nM. SC bioavailability
// F is applied to the absorption flux, so doses go in as nmol.

$PARAM @annotated
CL    : 0.20   : Linear clearance (L/day)
V1    : 3.0    : Central volume (L)
Q     : 0.50   : Intercompartmental clearance (L/day)
V2    : 2.5    : Peripheral volume (L)
KA    : 0.25   : SC absorption rate constant (1/day)
F     : 0.70   : SC bioavailability
MW    : 150000 : Molecular weight (g/mol)
KSS   : 1.0    : Quasi-steady-state constant (nM)
KINT  : 1.0    : Complex internalisation rate (1/day)
KDEG  : 0.2    : Free target degradation rate (1/day)
R0    : 2.0    : Baseline target (nM)
B0    : 100    : Baseline biomarker
KOUT  : 0.3    : Biomarker turnover rate (1/day)
GAMMA : 1.0    : Biomarker sensitivity to free target
IMAX  : 0.8    : Maximum fractional suppression of biomarker production

$CMT @annotated
DEPOT : Subcutaneous depot (nmol)
CENT  : Total drug in central compartment (nmol)
PERI  : Free drug in peripheral compartment (nmol)
RTOT  : Total target (nM)
BIO   : Biomarker

$MAIN
RTOT_0 = R0;
BIO_0 = B0;

$ODE
double Ctot = CENT / V1;
double b = Ctot - RTOT - KSS;
double disc = sqrt(b * b + 4.0 * KSS * Ctot);
double Cf = (b >= 0) ? 0.5 * (b + disc) : 2.0 * KSS * Ctot / (disc - b);
double Cplx = Ctot - Cf;
double Rfree = RTOT - Cplx;
if (Rfree < 0) Rfree = 0;

dxdt_DEPOT = -KA * DEPOT;
dxdt_CENT  = F * KA * DEPOT - CL * Cf - Q * Cf + Q * PERI / V2 - KINT * Cplx * V1;
dxdt_PERI  = Q * Cf - Q * PERI / V2;
dxdt_RTOT  = KDEG * R0 - KDEG * RTOT - (KINT - KDEG) * Cplx;
dxdt_BIO   = KOUT * B0 * (1.0 - IMAX + IMAX * pow(Rfree / R0, GAMMA)) - KOUT * BIO;

$TABLE
double CTOT = CENT / V1;
double bb = CTOT - RTOT - KSS;
double dd = sqrt(bb * bb + 4.0 * KSS * CTOT);
double CFREE = (bb >= 0) ? 0.5 * (bb + dd) : 2.0 * KSS * CTOT / (dd - bb);
double RO = (RTOT > 0) ? (CTOT - CFREE) / RTOT : 0;
double CONC = CFREE * MW / 1e6;

$CAPTURE @annotated
CFREE : Free drug (nM)
CONC  : Free drug (ug/mL)
RO    : Receptor occupancy (fraction)
