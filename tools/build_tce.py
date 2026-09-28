"""Build models/tce_mosunetuzumab_hosseini2020.cpp from the SimBiology model.

Usage: python tools/build_tce.py
Input: tools/tce_hosseini2020_source.json (from tools/extract_simbiology.py,
run on TDBr26_6_paper.sbproj in the Supplementary Software of Hosseini et
al., npj Syst Biol Appl 2020;6:28).

Parameters: the project's base values with the variants the authors'
MainRun_Par_paper.m activates, in order, plus the human physiology and human
PK variants. The result is checked against the paper's Supplementary Table 2
(human column) before anything is written.

Scope: the CD20xCD3 bispecific (mosunetuzumab) given intravenously. The
rituximab, blinatumomab, subcutaneous-dosing and BAFF sub-models are dropped
(their drugs are never dosed here; BAFF feeds nothing back).
"""
import ast
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = json.load(open(ROOT / 'tools' / 'tce_hosseini2020_source.json'))
OUT = ROOT / 'models' / 'tce_mosunetuzumab_hosseini2020.cpp'

# MainRun_Par_paper.m activates 5, 9, 11, 13, 14 and 24 for every study;
# 20 is the human physiology and 27 the human linear PK.
VARIANTS = [5, 9, 11, 13, 14, 20, 24, 27]
TUMOR_VARIANT = 28     # 'Tumor compartment - DLBCL', offered as a switch

# Supplementary Table 2, human column (param_table.xlsx)
PUBLISHED = dict(VmT=0.9, fTact=0.25, KmBT_act=0.716, ndrugactT=0.8, S=1.4, KdrugactT=130.083, VmB=0.95,
                 KmTB_kill=0.75, fKmTB_kill=10, nkill=1.027, KdrugB=1.302, Cl_tdb=5.4, Cld_tdb=24.17,
                 Vc_tdb=36.78, Vp_tdb=173.47, kBapop=0.02, kBkill=275.189, fBkill=0.1, kBprolif=0.05,
                 kTprolif=0.7, kTaexit=0.123, kTact=9.828, fTadeact=0.005, fTap=3.773, kTaapop=0.061,
                 fTaprolif=2.106, fTrapop=0.2, kTrexit=0.045, fdrug=6.3839, fa0=0.002, fAICD=1.497,
                 fTa0deact=0.001, fTa0apop=2.019, finj=1, Trpbref_perml=5e5, Bpbref_perml=5e5, Vpb=5000,
                 Vtissue=210, KTrp=200, KBp=333, Kp=0.14, Vtissue2=400, KTrp2=190, KBp2=190, Kp2=0.07,
                 Vtissue3=500, KTrp3=60, KBp3=80, Kp3=0.07, kBmat_kBapop_ratio=0.25, B19no20_B1920_ratio=0.25,
                 IL6_tiss_contribution=0.003, thalfIL6=20)

DROP_SPECIES = {'BAFF', 'RTXc_ugperkg', 'RTXp_ugperkg', 'Blinc_ug', 'TDBsc_ugperkg', 'TDBc_ugperml_AUC'}
ZERO = ['RTXc_ugperml', 'RTXt_ugperml', 'RTXt2_ugperml', 'RTXt3_ugperml', 'RTXtumor_ugperml',
        'Blinc_ngperml', 'Blint_ngperml', 'Blint2_ngperml', 'Blint3_ngperml', 'Blintumor_ngperml']


def parameters():
    p = dict(SRC['parameters'])
    for v in VARIANTS:
        p.update(SRC['variants'][v - 1]['values'])
    for k, v in PUBLISHED.items():
        assert abs(p[k] - v) <= 1e-9 * max(1, abs(v)), (k, p[k], v)
    return p


def side(s):
    out = {}
    for term in s.split('+'):
        term = term.strip()
        if not term or term == 'null':
            continue
        m = re.match(r'(\d+)\s+(\w+)$', term)
        k, name = (int(m.group(1)), m.group(2)) if m else (1, term)
        out[name] = out.get(name, 0) + k
    return out


def to_c(expr):
    """SimBiology expression -> C (pow for ^, fmax/fmin, double literals)."""
    tree = ast.parse(expr.replace('^', '**'), mode='eval').body

    def emit(n):
        if isinstance(n, ast.BinOp):
            l, r = emit(n.left), emit(n.right)
            if isinstance(n.op, ast.Pow):
                return 'pow(%s, %s)' % (l, r)
            op = {ast.Add: '+', ast.Sub: '-', ast.Mult: '*', ast.Div: '/'}[type(n.op)]
            return '(%s %s %s)' % (l, op, r)
        if isinstance(n, ast.UnaryOp):
            return '(-%s)' % emit(n.operand) if isinstance(n.op, ast.USub) else emit(n.operand)
        if isinstance(n, ast.Constant):
            return repr(float(n.value))
        if isinstance(n, ast.Name):
            return n.id
        if isinstance(n, ast.Call):
            f = {'max': 'fmax', 'min': 'fmin'}.get(n.func.id, n.func.id)
            return '%s(%s)' % (f, ', '.join(emit(a) for a in n.args))
        if isinstance(n, ast.Compare):
            op = {ast.Gt: '>', ast.Lt: '<', ast.GtE: '>=', ast.LtE: '<='}[type(n.ops[0])]
            return '(%s %s %s)' % (emit(n.left), op, emit(n.comparators[0]))
        raise ValueError(ast.dump(n))
    return emit(tree)


def names(expr):
    return set(re.findall(r'(?<![0-9.])[A-Za-z_][A-Za-z0-9_]*', expr))


def main():
    p = parameters()
    # reactions -> net stoichiometry
    rx = []
    for r in SRC['reactions']:
        if r['rate'].startswith('mw') or not r['rate'].strip():
            continue
        lhs, rhs = re.split(r'<->|->', r['reaction'])
        net = side(rhs)
        for k, v in side(lhs).items():
            net[k] = net.get(k, 0) - v
        net = {k: v for k, v in net.items() if v}
        if set(net) & DROP_SPECIES:
            continue
        rx.append((net, r['rate']))
    states = sorted({s for net, _ in rx for s in net}, key=lambda s: [x for x, _ in rx].index(next(n for n, _ in rx if s in n)))

    rules = {}
    for r in SRC['rules']:
        lhs, rhs = [x.strip() for x in r.split('=', 1)]
        lhs = lhs.split('.')[-1]
        rules[lhs] = rhs
    rules['TDBc_ugperml'] = '(TDBc_ugperkg/Vc_tdb>1e-5)*TDBc_ugperkg/Vc_tdb'   # PK_v26 with PKflag = 1
    for z in ZERO:
        rules[z] = '0'
    init = {k: v for k, v in rules.items() if k in states}
    assign = {k: v for k, v in rules.items() if k not in states and k not in DROP_SPECIES
              and k not in ('Baffconsumption',)}

    # order the repeated assignments by dependency, keeping only what is needed
    rates_text = ' '.join(r for _, r in rx)
    needed, frontier = set(), names(rates_text) | {'Bpb_perml', 'totTpb_perml', 'Tafraction_pb', 'IL6combo',
                                                    'totTtiss', 'totTtiss2', 'totTtiss3', 'Btiss_perml',
                                                    'Btiss2_perml', 'B1920tiss3_perml', 'TDBt_ugperml', 'Bpb_norm'}
    while frontier:
        n = frontier.pop()
        if n in assign and n not in needed:
            needed.add(n)
            frontier |= names(assign[n])
    order, done = [], set()
    while len(order) < len(needed):
        for n in sorted(needed - done):
            if not (names(assign[n]) & (needed - done - {n})):
                order.append(n); done.add(n)
    used = (names(rates_text) | set().union(*[names(assign[n]) for n in order]) |
            set().union(*[names(v) for v in init.values()]))
    pars = sorted((used - set(states) - set(assign) - {'max', 'min', 'log', 'exp'}) & set(p), key=str.lower)
    missing = (used - set(p) - set(states) - set(assign) - {'max', 'min', 'log', 'exp'})
    assert not missing, missing

    L = []
    L.append('$PROB\nCD20xCD3 T-cell engager (mosunetuzumab): T-cell activation, B-cell killing and IL-6\n')
    L.append('''// Hosseini I, Gadkar K, Stefanich E, Li CC, Sun LL, Chu YW, Ramanujan S.
// Mitigating the risk of cytokine release syndrome in a Phase I trial of
// CD20/CD3 bispecific antibody mosunetuzumab in NHL: impact of translational
// system modeling. npj Syst Biol Appl 2020;6:28 (open access).
//
// GENERATED by tools/build_tce.py from the authors' SimBiology project
// (Supplementary Software, TDBr26_6_paper.sbproj). Do not edit by hand.
//
// A bispecific antibody bridges CD3 on T cells and CD20 on B cells. T cells
// are activated in proportion to drug concentration and to the local B:T
// ratio, activated T cells kill B cells, proliferate, and die by
// activation-induced cell death; activated cells return to a post-activated
// state that can be reactivated. Each dose also causes transient T-cell
// margination out of blood (injection and drug effects). IL-6 is produced by
// activated T cells engaging B cells, so the first dose - when B cells are
// plentiful - drives the cytokine peak; step-up dosing depletes B cells at a
// low dose first. Compartments: peripheral blood (pb), spleen (tiss), lymph
// nodes (tiss2), bone marrow (tiss3, with CD19+CD20- precursors), and an
// optional tumour (tumor_on).
//
// Parameters: base project values with the variants activated in the
// authors' MainRun_Par_paper.m (%s), human physiology and human linear PK;
// checked against Supplementary Table 2 (human column) by the builder.
// kIL6prod was calibrated in cynomolgus monkeys only, so IL-6 is best read
// relative to its own peak. Rituximab, blinatumomab, subcutaneous dosing and
// BAFF are left out.
//
// Units: time day; cells (per compartment); drug amount ug/kg (dose TDBc_ugperkg,
// IV bolus), concentrations ug/mL. Give injection_effect and drug_effect an
// amount of 1 with every dose.
''' % ', '.join(map(str, VARIANTS)))
    L.append('$PARAM')
    tumor = SRC['variants'][TUMOR_VARIANT - 1]['values']
    for k in pars:
        v = p[k]
        L.append('%s = %r' % (k, float(v)))
    L.append('\n$CMT ' + ' '.join(states))
    L.append('\n$MAIN')
    for s in states:
        if s in init:
            L.append('%s_0 = %s;' % (s, to_c(init[s])))
    L.append('\n$ODE')
    for n in order:
        L.append('double %s = %s;' % (n, to_c(assign[n])))
    for s in states:
        terms = []
        for net, rate in rx:
            if s in net:
                c = net[s]
                r = to_c(rate)
                terms.append(('+ ' if c > 0 else '- ') + ('%d * ' % abs(c) if abs(c) != 1 else '') + r)
        L.append('dxdt_%s = %s;' % (s, ' '.join(terms).lstrip('+ ') if terms else '0'))
    outputs = ['TDBc_ugperml', 'Bpb_perml', 'totTpb_perml', 'Tafraction_pb', 'IL6combo', 'Btiss_perml',
               'Btiss2_perml', 'B1920tiss3_perml', 'Bpb_norm', 'Btumor_perml']
    L.append('\n$TABLE')
    for n in order:
        L.append('%s = %s;' % (n, to_c(assign[n])))
    L.append('\n$CAPTURE ' + ' '.join(o for o in outputs if o in order))
    with open(OUT, 'w', newline='\n') as f:
        f.write('\n'.join(L) + '\n')
    print(len(states), 'states,', len(rx), 'reactions,', len(order), 'assignments,', len(pars), 'parameters')
    print('tumour variant:', tumor)


if __name__ == '__main__':
    main()
