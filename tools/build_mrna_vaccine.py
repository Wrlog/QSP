"""Build models/mrna_vaccine_dasti2025.cpp from the COSBI MATLAB model.

Usage: python tools/build_mrna_vaccine.py

Source: the tissue layer of the multiscale QSP model of mRNA vaccines,
github.com/cosbi-research/QSPmRNAVaccines (tissue_layer/model_equations.m,
Parameters.m, simulation.m and the fitted parameter vectors Liang_fit.mat and
BNT162b2_fit.mat / BNT162b2_fit_over60.mat / mRNA-1273_fit.mat), published as
Dasti et al., CPT Pharmacometrics Syst Pharmacol 2025. COSBI-SSLA licence:
non-commercial use; this derivative carries the same licence.

The MATLAB code is translated line by line. Two changes, both numerical:
- The free antigen concentration in the lymph node, found with fzero in the
  original, is found with a fixed number of Newton steps from zero. The
  equation is increasing and concave, so the iterates rise monotonically to
  the root without overshoot; the test suite checks them against a
  bracketing solver.
- MATLAB restarts time at 0 after the second dose, so its "t < 0.1" guard on
  the T-cell proliferation function applies after each dose; here that is
  written with the time since the last dose (doses at 0, TD2 and TD3).
"""
from math import exp, log, pi, sqrt
from pathlib import Path

J = 17
OUT = Path(__file__).resolve().parent.parent / 'models' / 'mrna_vaccine_dasti2025.cpp'

LIANG = ['kdeg', 'K_mRNA', 'SF_Gamma', 'OmegaNP', 'EtaNP', 'GammaNP', 'DeltaNP', 'MuNP', 'XiNP',
         'OmegaMN', 'EtaMN', 'GammaMN', 'DeltaMN', 'MuMN', 'XiMN', 'OmegamDC', 'EtamDC', 'GammamDC',
         'DeltamDC', 'MumDC', 'XimDC', 'OmegapDC', 'EtapDC', 'GammapDC', 'DeltapDC', 'MupDC', 'XipDC',
         'BetamDC', 'BetapDC']
LIANG_FIT = [1.1869218703e-01, 2.0724385730e-01, 2.9587638911e+03, 4.1028128869e-04, 8.6040977998e-04,
             2.3216956755e-06, 1.8223457287e-01, 2.2638280770e-02, 1.1127157934e-01, 6.0836876827e-03,
             3.8097713758e-03, 4.0593260315e-06, 4.3040377103e+00, 5.2198127317e-01, 1.1946876459e-04,
             2.4491812268e-01, 5.2203808310e-03, 5.0855147569e-06, 4.9219018002e+00, 1.4684932978e-01,
             2.3517215461e-01, 1.5486738954e-01, 4.9994470187e-04, 4.4309199931e-06, 2.9159871988e+01,
             2.7072645268e+01, 3.9174034080e-03, 1.1402524274e+00, 1.3501733075e+00]
AB = ['CC_N', 'SigmaMB', 'RhoAB_N', 'RhoAB_M', 'LambdaSP', 'BetaA', 'expAg', 'SigmaNT', 'SigmaMT', 'RhoAT',
      'BetaLP', 'BetaSP', 'mRNA_max', 'ksat', 'delayB']
BNT162B2_FIT = [6.6178599775e+01, 8.6211997917e+00, 6.0612000850e+00, 6.3442991786e+00, 4.0482099838e+01,
                9.1799901351e-02, 3.3713258652e+05, 2.9465880005e+02, 9.3584379878e+02, 4.9808998282e+00,
                1.4300002628e-02, 1.3830001542e-01, 3.1416002938e+00, 5.5832999644e+01, 7.9505908077e-07]
# Liang et al. data, first row: initial injection-site cells (Liang_data.mat)
MN0_IS, NP0_IS, MDC0_IS, PDC0_IS = 300.0, 0.0, 1864.373, 0.0

FIXED = dict(V_LN=5e-4, V_BL=5.0, NAV=6.022e23, kdmrna=log(2) / (10 / 24), k_slope=1.0,
             BetaNP=24 * log(2) / 7, NP0_BL=3.19e9 * 5, BetaMN=24 * log(2) / 24, MN0_BL=(461 + 291) / 2 * 1e6 * 5,
             BetamIDC=0.0924, mIDC0_BL=(9.97 + 14.80) / 2 * 1e6 * 5, BetapIDC=0.0924, pIDC0_BL=7e6 * 5,
             NP0_IS=NP0_IS, MN0_IS=MN0_IS, mIDC0_IS=MDC0_IS, pIDC0_IS=PDC0_IS,
             p_L=0.2, p_M=0.65, p_H=1.0,
             k_trM_mDC=5.984386009847674, k_trH_mDC=4.023456291503178, k_atrH_mDC=2.001144541718650,
             k_atrM_mDC=1.182525483399790, k_atrL_mDC=0.393936426260313, k_trM_pDC=5.985073797735384,
             k_trH_pDC=4.025034796217001, k_atrH_pDC=2.000405122440113, k_atrM_pDC=1.182561965521271,
             k_atrL_pDC=0.393949613591618, NmaxMHC_mDC=1.226433347931808e+04, NmaxMHC_pDC=1.226433347931808e+04,
             NT0=1.445e3, KNT=400.0, KMT=40.0, BetaNT=0.0029, BetaAT=0.18, BetaMT=2.7397e-4, BetaFT=0.18, f1=0.5,
             BRN=75000.0, K_R=1.0, SigmaNB=2.48, BetaAB=0.2518, BetaMB=7.8278e-5, BetaNB=0.029, g1=0.5, g2=0.4,
             AlphaA=8.64e8, MW_Ab=150000.0)


def nb0():
    x = [-3 + 6 / 16 * i for i in range(J)]
    d = [exp(-v * v / 2) / sqrt(2 * pi) for v in x]
    s = sum(d)
    return [v / s * 1.3e9 * 0.000004 for v in d]


KA = [1e-6 * 2 ** (j - (J + 1) / 2) for j in range(1, J + 1)]


def main():
    dc_is = ['DC_IS', 'DCL_IS', 'DCAgLon_IS', 'DCAgMon_IS', 'DCAgH_IS', 'DCAgMoff_IS', 'DCAgLoff_IS', 'DCoff_IS']
    dc_ln = ['DCL_LN', 'DCAgLon_LN', 'DCAgMon_LN', 'DCAgH_LN', 'DCAgMoff_LN', 'DCAgLoff_LN', 'DC_LN']
    cmt = (['mRNA', 'NP_IS', 'NPL_IS', 'NPAg_IS', 'MN_IS', 'MNL_IS', 'MNAg_IS'] +
           ['m' + s for s in dc_is] + ['p' + s for s in dc_is] +
           ['NPL_LN', 'NPAg_LN', 'MNL_LN', 'MNAg_LN'] + ['m' + s for s in dc_ln] + ['p' + s for s in dc_ln] +
           ['NT', 'AT_N', 'MT', 'AT_M', 'FT'])
    bnames = ['NB', 'AB_N', 'MB', 'AB_M', 'SP', 'LP']
    for b in bnames:
        cmt += ['%s%d' % (b, j) for j in range(1, J + 1)]
    cmt += ['NP_BL', 'MN_BL', 'mDC_BL', 'pDC_BL']
    for b in ['SP_BL', 'LP_BL', 'Ab_BL', 'GCB']:
        cmt += ['%s%d' % (b, j) for j in range(1, J + 1)]
    assert len(cmt) == 50 + 10 * J
    rng = range(1, J + 1)
    L = []
    w = L.append
    w('$PROB\nmRNA vaccine: innate uptake, antigen presentation, T and B cells, antibodies\n')
    w('''// Multiscale QSP model for mRNA vaccines, COSBI, Fondazione The Microsoft
// Research - University of Trento Centre for Computational and Systems Biology.
// Dasti A, et al. A multiscale quantitative systems pharmacology model for the
// development and optimization of mRNA vaccines. CPT Pharmacometrics Syst
// Pharmacol 2025. Code: github.com/cosbi-research/QSPmRNAVaccines (tissue layer).
//
// Copyright (c) 2024, Fondazione The Microsoft Research - University of Trento
// Centre for Computational and Systems Biology (COSBI). All rights reserved.
// Distributed under the COSBI Licence Terms (COSBI-SSLA), for non-commercial
// purposes only; see apps/mrna-vaccine/LICENSE. THE SOFTWARE COMES "AS IS",
// WITH NO WARRANTIES, AND COSBI AND ITS CONTRIBUTORS ARE NOT LIABLE FOR ANY
// DAMAGES RELATED TO IT (see the licence for the full terms).
//
// MODIFIED, 2026-09-28: the MATLAB tissue layer (model_equations.m,
// Parameters.m) translated to an mrgsolve model file by
// tools/build_mrna_vaccine.py. Free lymph-node antigen found by Newton
// iteration instead of fzero; the t < 0.1 proliferation guard written with the
// time since the last dose. GENERATED - do not edit by hand.
//
// Injection site (IS): LNPs are taken up by neutrophils, monocytes, myeloid and
// plasmacytoid dendritic cells, which express the antigen and migrate to the
// lymph node (LN). Dendritic cells mature through low/medium/high antigen
// presentation. Presented antigen activates helper T cells; free antigen binds
// B-cell receptors of 17 affinity classes (Ka = 1e-6 * 2^(j-9) per pM), which
// with T-cell help become germinal-centre, memory, and short- and long-lived
// plasma cells that secrete antibody into blood.
//
// Units: time day; mRNA pmol (dose into mRNA: 1e12 * grams / 1377479.8);
// cells; antibody pmol per class in blood. IGG output in ng/mL.
''')
    w('$PARAM')
    for k, v in zip(LIANG, LIANG_FIT):
        w('%s = %r' % (k, v))
    for k, v in zip(AB, BNT162B2_FIT):
        w('%s = %r' % (k, v))
    for k, v in FIXED.items():
        w('%s = %r' % (k, float(v)))
    w('TD2 = 1e9\nTD3 = 1e9')
    w('\n$CMT ' + ' '.join(cmt))
    w('\n$MAIN')
    w('double AlphaNP_IS = NP0_IS * (BetaNP + OmegaNP);')
    w('double AlphaNP_BL = -OmegaNP * NP0_IS + BetaNP * NP0_BL;')
    w('double AlphaMN_IS = MN0_IS * (BetaMN + OmegaMN);')
    w('double AlphaMN_BL = -OmegaMN * MN0_IS + BetaMN * MN0_BL;')
    w('double AlphamIDC_IS = mIDC0_IS * (BetamIDC + OmegamDC);')
    w('double AlphamIDC_BL = -OmegamDC * mIDC0_IS + BetamIDC * mIDC0_BL;')
    w('double AlphapIDC_IS = pIDC0_IS * (BetapIDC + OmegapDC);')
    w('double AlphapIDC_BL = -OmegapDC * pIDC0_IS + BetapIDC * pIDC0_BL;')
    w('double maxMHC_mDC = NmaxMHC_mDC / NAV * 1e12;')
    w('double maxMHC_pDC = NmaxMHC_pDC / NAV * 1e12;')
    w('double CC_M = CC_N * 10.0;')
    w('double LambdaLP = LambdaSP;')
    w('NP_IS_0 = NP0_IS;\nMN_IS_0 = MN0_IS;\nmDC_IS_0 = mIDC0_IS;\npDC_IS_0 = pIDC0_IS;\nNT_0 = NT0;')
    for j, v in zip(rng, nb0()):
        w('NB%d_0 = %r;' % (j, v))
    w('NP_BL_0 = NP0_BL;\nMN_BL_0 = MN0_BL;\nmDC_BL_0 = mIDC0_BL;\npDC_BL_0 = pIDC0_BL;')

    w('\n$ODE')
    for x in 'mp':
        w('double %sDCAgL = %sDCAgLon_LN + %sDCAgLoff_LN;' % (x, x, x))
        w('double %sDCAgM = %sDCAgMon_LN + %sDCAgMoff_LN;' % (x, x, x))
        w('double %sDC_mat0 = %sDCAgL * p_L + %sDCAgM * p_M + %sDCAgH_LN * p_H;' % (x, x, x, x))
        w('double %sDC_tot = %sDCAgL + %sDCAgM + %sDCAgH_LN;' % (x, x, x, x))
        w('double %sDC_mat = %sDC_mat0 < p_L ? 0.0 : %sDC_mat0;' % (x, x, x))
        w('double Nmedio_%sDC = %sDC_mat < p_L ? 0.0 : %sDC_mat * NmaxMHC_%sDC / %sDC_tot;' % (x, x, x, x, x))
        w('double sum_%sDC_T = %sDC_tot / (%sDC_tot + NT + AT_N + AT_M + MT);' % (x, x, x))
        w('double D_%sDC_N = sum_%sDC_T * Nmedio_%sDC / (Nmedio_%sDC + KNT);' % (x, x, x, x))
        w('double D_%sDC_M = sum_%sDC_T * Nmedio_%sDC / (Nmedio_%sDC + KMT);' % (x, x, x, x))
        w('double E_%sDC_N = sum_%sDC_T * (Nmedio_%sDC - KNT) / (Nmedio_%sDC + KNT);' % (x, x, x, x))
        w('double E_%sDC_M = sum_%sDC_T * (Nmedio_%sDC - KMT) / (Nmedio_%sDC + KMT);' % (x, x, x, x))
    w('double DCsum = mDC_tot + pDC_tot;')
    w('double weight_mDC = DCsum > 0 ? mDC_tot / DCsum : 0.0;')
    w('double weight_pDC = DCsum > 0 ? pDC_tot / DCsum : 0.0;')
    w('double D_N = D_mDC_N * weight_mDC + D_pDC_N * weight_pDC;')
    w('double D_M = D_mDC_M * weight_mDC + D_pDC_M * weight_pDC;')
    w('double E_N0 = E_mDC_N * weight_mDC + E_pDC_N * weight_pDC;')
    w('double E_M0 = E_mDC_M * weight_mDC + E_pDC_M * weight_pDC;')
    w('double tlast2 = SOLVERTIME >= TD2 ? TD2 : 0.0;')
    w('double tlast = SOLVERTIME >= TD3 ? TD3 : tlast2;')
    w('double E_N = SOLVERTIME - tlast < 0.1 ? fmax(E_N0, 0.0) : E_N0;')
    w('double E_M = SOLVERTIME - tlast < 0.1 ? fmax(E_M0, 0.0) : E_M0;')
    # B-cell receptor occupancy
    w('double Bmult = BRN / NAV * 1e12;')
    for j in rng:
        w('double BCR%d = (NB%d + AB_N%d + GCB%d + AB_M%d + MB%d) * Bmult;' % (j, j, j, j, j, j))
    w('double Ag = expAg * (maxMHC_mDC * ((fabs(p_L - p_M) / 2 + p_L) * mDCAgL + (fabs(p_M - p_H) / 2 + p_M) * mDCAgM '
      '+ (fabs(p_H - 1.0) / 2 + p_H) * mDCAgH_LN) + maxMHC_pDC * ((fabs(p_L - p_M) / 2 + p_L) * pDCAgL '
      '+ (fabs(p_M - p_H) / 2 + p_M) * pDCAgM + (fabs(p_H - 1.0) / 2 + p_H) * pDCAgH_LN));')
    w('double Agc = Ag / V_LN;')
    for j, k in zip(rng, KA):
        w('double cB%d = %r * BCR%d / V_LN;' % (j, k, j))
    csum = ' + '.join('cB%d' % j for j in rng)
    w('double Agf0 = Agc / (1.0 + %s);' % csum)
    NEWTON = 8
    for i in range(NEWTON):
        x = 'Agf%d' % i
        g = ' + '.join('cB%d / (1.0 + %r * %s)' % (j, k, x) for j, k in zip(rng, KA))
        dg = ' + '.join('cB%d / pow(1.0 + %r * %s, 2)' % (j, k, x) for j, k in zip(rng, KA))
        w('double Agf%d = %s - (%s * (1.0 + %s) - Agc) / (1.0 + %s);' % (i + 1, x, x, g, dg))
    w('double Agfree = Agf%d;' % NEWTON)
    for j, k in zip(rng, KA):
        w('double ro%d = %r * Agfree / (1.0 + %r * Agfree);' % (j, k, k))
        w('double R%d = ro%d * BRN;' % (j, j))
        w('double F%d = R%d / (K_R + R%d);' % (j, j, j))
        w('double G%d = (1.0 - ro%d) * F%d;' % (j, j, j))
        w('double H%d = (R%d - K_R) / (R%d + K_R);' % (j, j, j))
    bsum = ' + '.join('NB{0} + AB_N{0} + GCB{0} + AB_M{0} + MB{0}'.format(j) for j in rng)
    w('double Bsum = %s;' % bsum)
    w('double P_N = CC_N * FT / (CC_N * FT + Bsum);')
    w('double P_M = CC_M * FT / (CC_M * FT + Bsum);')
    w('double upt = mRNA / (K_mRNA + mRNA);')
    # ODEs
    w('dxdt_mRNA = -GammapDC * mRNA * pDC_IS - GammamDC * mRNA * mDC_IS - GammaMN * mRNA * MN_IS - GammaNP * mRNA * NP_IS'
      ' - kdmrna * mRNA - kdeg * mRNA - ksat * mRNA * 0.5 * (1.0 + tanh((mRNA - mRNA_max) / k_slope));')
    for c, a, bb in (('NP', 'AlphaNP_IS', 'BetaNP'), ('MN', 'AlphaMN_IS', 'BetaMN')):
        w('dxdt_{0}_IS = {1} + Eta{0} * upt * {0}_BL - Gamma{0} * SF_Gamma * mRNA * {0}_IS - {2} * {0}_IS - Omega{0} * {0}_IS;'.format(c, a, bb))
        w('dxdt_{0}L_IS = Gamma{0} * SF_Gamma * mRNA * {0}_IS - Delta{0} * {0}L_IS - {1} * {0}L_IS - Mu{0} * {0}L_IS;'.format(c, bb))
        w('dxdt_{0}Ag_IS = Delta{0} * {0}L_IS - {1} * {0}Ag_IS - Xi{0} * {0}Ag_IS;'.format(c, bb))
    for x in 'mp':
        X = x + 'DC'
        w('dxdt_{0}_IS = Alpha{1}IDC_IS + Eta{0} * upt * {0}_BL - Gamma{0} * SF_Gamma * mRNA * {0}_IS - Beta{1}IDC * {0}_IS - Omega{0} * {0}_IS;'.format(X, x))
        w('dxdt_{0}L_IS = Gamma{0} * SF_Gamma * mRNA * {0}_IS - Delta{0} * {0}L_IS - Beta{0} * {0}L_IS - Mu{0} * {0}L_IS;'.format(X))
        w('dxdt_{0}AgLon_IS = Delta{0} * {0}L_IS - k_trM_{0} * {0}AgLon_IS - Beta{0} * {0}AgLon_IS - Xi{0} * {0}AgLon_IS;'.format(X))
        w('dxdt_{0}AgMon_IS = k_trM_{0} * {0}AgLon_IS - k_trH_{0} * {0}AgMon_IS - Beta{0} * {0}AgMon_IS - Xi{0} * {0}AgMon_IS;'.format(X))
        w('dxdt_{0}AgH_IS = k_trH_{0} * {0}AgMon_IS - k_atrH_{0} * {0}AgH_IS - Beta{0} * {0}AgH_IS - Xi{0} * {0}AgH_IS;'.format(X))
        w('dxdt_{0}AgMoff_IS = k_atrH_{0} * {0}AgH_IS - k_atrM_{0} * {0}AgMoff_IS - Beta{0} * {0}AgMoff_IS - Xi{0} * {0}AgMoff_IS;'.format(X))
        w('dxdt_{0}AgLoff_IS = k_atrM_{0} * {0}AgMoff_IS - k_atrL_{0} * {0}AgLoff_IS - Beta{0} * {0}AgLoff_IS - Xi{0} * {0}AgLoff_IS;'.format(X))
        w('dxdt_{0}off_IS = k_atrL_{0} * {0}AgLoff_IS - Beta{0} * {0}off_IS - Xi{0} * {0}off_IS;'.format(X))
    for c, bb in (('NP', 'BetaNP'), ('MN', 'BetaMN')):
        w('dxdt_{0}L_LN = Mu{0} * {0}L_IS - Delta{0} * {0}L_LN - {1} * {0}L_LN;'.format(c, bb))
        w('dxdt_{0}Ag_LN = Xi{0} * {0}Ag_IS + Delta{0} * {0}L_LN - {1} * {0}Ag_LN;'.format(c, bb))
    for x in 'mp':
        X = x + 'DC'
        w('dxdt_{0}L_LN = Mu{0} * {0}L_IS - Delta{0} * {0}L_LN - Beta{0} * {0}L_LN;'.format(X))
        w('dxdt_{0}AgLon_LN = Xi{0} * {0}AgLon_IS + Delta{0} * {0}L_LN - k_trM_{0} * {0}AgLon_LN - Beta{0} * {0}AgLon_LN;'.format(X))
        w('dxdt_{0}AgMon_LN = Xi{0} * {0}AgMon_IS + k_trM_{0} * {0}AgLon_LN - k_trH_{0} * {0}AgMon_LN - Beta{0} * {0}AgMon_LN;'.format(X))
        w('dxdt_{0}AgH_LN = Xi{0} * {0}AgH_IS + k_trH_{0} * {0}AgMon_LN - k_atrH_{0} * {0}AgH_LN - Beta{0} * {0}AgH_LN;'.format(X))
        w('dxdt_{0}AgMoff_LN = Xi{0} * {0}AgMoff_IS + k_atrH_{0} * {0}AgH_LN - k_atrM_{0} * {0}AgMoff_LN - Beta{0} * {0}AgMoff_LN;'.format(X))
        w('dxdt_{0}AgLoff_LN = Xi{0} * {0}AgLoff_IS + k_atrM_{0} * {0}AgMoff_LN - k_atrL_{0} * {0}AgLoff_LN - Beta{0} * {0}AgLoff_LN;'.format(X))
        w('dxdt_{0}_LN = Xi{0} * {0}off_IS + k_atrL_{0} * {0}AgLoff_LN - Beta{0} * {0}_LN;'.format(X))
    w('dxdt_NT = BetaNT * (NT0 - NT) - SigmaNT * D_N * NT;')
    w('dxdt_AT_N = SigmaNT * D_N * NT + RhoAT * E_N * AT_N - BetaAT * AT_N;')
    w('dxdt_MT = RhoAT * (1.0 - E_N) * f1 * AT_N + RhoAT * (1.0 - E_M) * f1 * AT_M - SigmaMT * D_M * MT - BetaMT * MT;')
    w('dxdt_AT_M = SigmaMT * D_M * MT + RhoAT * E_M * AT_M - BetaAT * AT_M;')
    w('dxdt_FT = RhoAT * (1.0 - E_N) * (1.0 - f1) * AT_N + RhoAT * (1.0 - E_M) * (1.0 - f1) * AT_M - BetaFT * FT;')
    for j in rng:
        w('dxdt_NB{0} = BetaNB * ({1!r} - NB{0}) - SigmaNB * F{0} * P_N * NB{0};'.format(j, nb0()[j - 1]))
        w('dxdt_AB_N{0} = SigmaNB * G{0} * P_N * NB{0} - delayB * AB_N{0} - BetaAB * AB_N{0};'.format(j))
        w('dxdt_GCB{0} = delayB * AB_N{0} + RhoAB_N * H{0} * P_N * GCB{0} - BetaAB * GCB{0};'.format(j))
        w('dxdt_MB{0} = RhoAB_N * (1.0 - H{0}) * P_N * g1 * GCB{0} + RhoAB_M * (1.0 - H{0}) * P_M * g1 * AB_M{0} - SigmaMB * F{0} * P_M * MB{0} - BetaMB * MB{0};'.format(j))
        w('dxdt_AB_M{0} = SigmaMB * G{0} * P_M * MB{0} + RhoAB_M * H{0} * P_M * AB_M{0} - BetaAB * AB_M{0};'.format(j))
        w('dxdt_SP{0} = RhoAB_N * (1.0 - H{0}) * P_N * g2 * GCB{0} + RhoAB_M * (1.0 - H{0}) * P_M * g2 * AB_M{0} - (BetaSP + LambdaSP) * SP{0};'.format(j))
        w('dxdt_LP{0} = RhoAB_N * (1.0 - H{0}) * P_N * (1.0 - g1 - g2) * GCB{0} + RhoAB_M * (1.0 - H{0}) * P_M * (1.0 - g1 - g2) * AB_M{0} - (BetaLP + LambdaLP) * LP{0};'.format(j))
    w('dxdt_NP_BL = AlphaNP_BL - EtaNP * upt * NP_BL - BetaNP * NP_BL + OmegaNP * NP_IS;')
    w('dxdt_MN_BL = AlphaMN_BL - EtaMN * upt * MN_BL - BetaMN * MN_BL + OmegaMN * MN_IS;')
    w('dxdt_mDC_BL = AlphamIDC_BL - EtamDC * upt * mDC_BL - BetamIDC * mDC_BL + OmegamDC * mDC_IS;')
    w('dxdt_pDC_BL = AlphapIDC_BL - EtapDC * upt * pDC_BL - BetapIDC * pDC_BL + OmegapDC * pDC_IS;')
    for j in rng:
        w('dxdt_SP_BL{0} = -BetaSP * SP_BL{0} + LambdaSP * SP{0};'.format(j))
        w('dxdt_LP_BL{0} = -BetaLP * LP_BL{0} + LambdaLP * LP{0};'.format(j))
        w('dxdt_Ab_BL{0} = -BetaA * Ab_BL{0} + AlphaA * (SP_BL{0} / NAV * 1e12) + AlphaA * (LP_BL{0} / NAV * 1e12);'.format(j))

    w('\n$TABLE')
    abt = ' + '.join('Ab_BL%d' % j for j in rng)
    w('double ABTOT = %s;' % abt)
    w('double IGG = ABTOT * 1e-12 * MW_Ab * 1e9 / (V_BL * 1e3);')
    w('double AFFINITY = ABTOT > 0 ? (%s) / ABTOT : 0.0;' % ' + '.join('%d * Ab_BL%d' % (j - 9, j) for j in rng))
    w('double GCBTOT = %s;' % ' + '.join('GCB%d' % j for j in rng))
    w('double MBTOT = %s;' % ' + '.join('MB%d' % j for j in rng))
    w('double SPTOT = %s;' % ' + '.join('SP_BL%d + SP%d' % (j, j) for j in rng))
    w('double LPTOT = %s;' % ' + '.join('LP_BL%d + LP%d' % (j, j) for j in rng))
    w('double APC_IS = NPAg_IS + MNAg_IS + mDCAgLon_IS + mDCAgMon_IS + mDCAgH_IS + mDCAgMoff_IS + mDCAgLoff_IS'
      ' + pDCAgLon_IS + pDCAgMon_IS + pDCAgH_IS + pDCAgMoff_IS + pDCAgLoff_IS;')
    w('double APC_LN = mDCAgLon_LN + mDCAgMon_LN + mDCAgH_LN + mDCAgMoff_LN + mDCAgLoff_LN'
      ' + pDCAgLon_LN + pDCAgMon_LN + pDCAgH_LN + pDCAgMoff_LN + pDCAgLoff_LN;')
    w('double THELP = FT;')
    w('double TACT = AT_N + AT_M;')
    w('double TMEM = MT;')
    w('\n$CAPTURE IGG AFFINITY GCBTOT MBTOT SPTOT LPTOT APC_IS APC_LN THELP TACT TMEM')
    with open(OUT, 'w', newline='\n') as f:
        f.write('\n'.join(L) + '\n')
    print(len(cmt), 'states written to', OUT.name)


if __name__ == '__main__':
    main()
