"""Rebuild models/coagulation_wajima2009.cpp from BioModels BIOMD0000000340.

Download the SBML first:
  curl -L "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000340?filename=BIOMD0000000340_url.xml" -o BIOMD0000000340.xml
then: python tools/build_coagulation.py BIOMD0000000340.xml
"""
import sys
sys.path.insert(0, "tools")
from sbml_to_mrgsolve import Converter

c = Converter(sys.argv[1], rename={"F": "Fibrin", "D": "Ddimer", "P": "Plasmin", "II": "FII"})
# heparin_infusion is switched by an SBML event; here it is an ordinary
# parameter, a constant input.
for p in c.m.findall("s:listOfParameters/s:parameter", c.ns):
    if p.get("id") == "heparin_infusion":
        p.set("constant", "true")
notes = [
    "Humoral coagulation network with warfarin, vitamin K and heparin.",
    "Wajima T, Isbister GK, Duffull SB. A comprehensive model for the humoral",
    "coagulation network in humans. Clin Pharmacol Ther 2009;86:290-298.",
    "",
    "Converted by tools/build_coagulation.py from BioModels BIOMD0000000340",
    "(CC0). Species F, D, P and II are renamed Fibrin, Ddimer, Plasmin and",
    "FII (II is reserved in mrgsolve).",
    "The SBML dosing events are replaced by dose records: warfarin (mg) into",
    "A_warf; heparin through the heparin_infusion parameter (nM/h).",
    "AUC_Fibrin (added) integrates fibrin for the in-vitro clotting tests.",
    "",
    "Units: time in hours, factors in nM, warfarin in mg (A_warf) and mg/L (C_warf).",
]
txt = c.convert("Wajima 2009 coagulation network", notes)
txt = txt.replace("C_warf           : C_warf (nM)", "C_warf           : Warfarin plasma concentration (mg/L)")
txt = txt.replace("A_warf           : A_warf (nM)", "A_warf           : Warfarin in the gut (mg)")
txt = txt.replace("\n$MAIN", "AUC_Fibrin       : Integral of fibrin, for the clotting-time tests (nM*h)\n\n$MAIN")
txt = txt.rstrip("\n") + "\ndxdt_AUC_Fibrin = Fibrin;\n"
# Initial values are set in $MAIN, which in mrgsolve overrides init(). The
# in-vitro clotting tests start from a given plasma state, so guard them.
head, rest = txt.split("\n$MAIN\n", 1)
main, ode = rest.split("\n$ODE\n", 1)
inits = [l for l in main.splitlines() if l.strip()]
main = "if (INVITRO < 0.5) {\n" + "\n".join("  " + l for l in inits) + "\n}\n"
txt = head + "\n$MAIN\n" + main + "\n$ODE\n" + ode
txt = txt.replace("$PARAM @annotated\n", "$PARAM @annotated\nINVITRO                      : 0            : 1 = start from a supplied plasma state (clotting tests)\n", 1)
open("models/coagulation_wajima2009.cpp", "w").write(txt)
print("wrote models/coagulation_wajima2009.cpp")
