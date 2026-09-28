"""Convert an SBML (level 2) reaction network into an mrgsolve model file.

Used to build models/coagulation_wajima2009.cpp from BioModels
BIOMD0000000340 (Wajima et al., Clin Pharmacol Ther 2009;86:290). Kept in
the repository so the conversion can be repeated and checked.

What it does:
  - species become compartments; initial values come from the species or
    from initial assignments (written as CMT_0 = ... in $MAIN)
  - global parameters and each reaction's local parameters (prefixed with
    the reaction id) become $PARAM entries
  - function definitions are expanded inline
  - rate rules become ODEs; assignment rules become $ODE variables
  - events are dropped: dosing is done with mrgsolve dose records instead

Usage: python sbml_to_mrgsolve.py model.xml out.cpp [--rename F=Fibrin,...]
"""

import re
import sys
import xml.etree.ElementTree as ET

M = "{http://www.w3.org/1998/Math/MathML}"


def ns_for(root):
    return {"s": root.tag.split("}")[0][1:], "m": M[1:-1]}


class Converter:
    def __init__(self, path, rename=None, drop_params=()):
        root = ET.parse(path).getroot()
        self.ns = ns_for(root)
        self.m = root.find("s:model", self.ns)
        self.rename = rename or {}
        self.drop_params = set(drop_params)
        self.functions = {}
        for f in self.m.findall("s:listOfFunctionDefinitions/s:functionDefinition", self.ns):
            lam = f.find("m:math", self.ns)[0]
            args = [b.find(M + "ci").text.strip() for b in lam.findall(M + "bvar")]
            body = [c for c in lam if c.tag != M + "bvar"][0]
            self.functions[f.get("id")] = (args, body)
        self.compartments = {c.get("id"): c.get("size") for c in
                             self.m.findall("s:listOfCompartments/s:compartment", self.ns)}

    def name(self, x):
        x = x.strip()
        return self.rename.get(x, x)

    def expr(self, e, local=None, subst=None):
        local = local or {}
        subst = subst or {}
        tag = e.tag.replace(M, "")
        if tag == "math":
            return self.expr(list(e)[0], local, subst)
        if tag == "ci":
            v = e.text.strip()
            if v in subst:
                return subst[v]
            if v in local:
                return local[v]
            if v in self.compartments:
                return "1.0" if float(self.compartments[v]) == 1 else self.compartments[v]
            return self.name(v)
        if tag == "cn":
            txt = e.text.strip() if e.text else ""
            if e.get("type") == "e-notation":
                parts = [t.strip() for t in e.itertext() if t.strip()]
                txt = "%se%s" % (parts[0], parts[1])
            if re.fullmatch(r"-?\d+", txt):
                txt += ".0"
            return txt
        if tag == "csymbol":
            return "SOLVERTIME"
        if tag == "apply":
            ch = list(e)
            op = ch[0].tag.replace(M, "")
            args = [self.expr(c, local, subst) for c in ch[1:]]
            if op == "ci":
                fname = ch[0].text.strip()
                fargs, body = self.functions[fname]
                return "(" + self.expr(body, {}, dict(zip(fargs, args))) + ")"
            if op == "plus":
                return "(" + " + ".join(args) + ")"
            if op == "times":
                return "(" + " * ".join(args) + ")"
            if op == "minus":
                return "(-" + args[0] + ")" if len(args) == 1 else "(" + " - ".join(args) + ")"
            if op == "divide":
                return "(" + args[0] + " / " + args[1] + ")"
            if op == "power":
                return "pow(" + args[0] + ", " + args[1] + ")"
            if op in ("exp", "log", "sqrt"):
                return op + "(" + args[0] + ")"
            raise ValueError("unsupported MathML operator: " + op)
        raise ValueError("unsupported MathML element: " + tag)

    def convert(self, title, header_notes):
        ns = self.ns
        species = self.m.findall("s:listOfSpecies/s:species", ns)
        params = [p for p in self.m.findall("s:listOfParameters/s:parameter", ns)
                  if p.get("constant") != "false" and p.get("id") not in self.drop_params]
        rules = self.m.find("s:listOfRules", ns)
        rate_rules = [r for r in rules if r.tag.endswith("rateRule")] if rules is not None else []
        assign_rules = [r for r in rules if r.tag.endswith("assignmentRule")] if rules is not None else []
        assigned = {r.get("variable") for r in assign_rules}
        init_assign = {a.get("symbol"): a for a in
                       self.m.findall("s:listOfInitialAssignments/s:initialAssignment", ns)}

        cmts = [s for s in species if s.get("id") not in assigned]
        extra_state = [r.get("variable") for r in rate_rules
                       if r.get("variable") not in {s.get("id") for s in species}]

        out = ["$PROB", title, ""]
        out += ["// " + l if l else "//" for l in header_notes]
        out += ["", "$PARAM @annotated"]
        for p in params:
            out.append("%-28s : %-12s : %s" % (self.name(p.get("id")), p.get("value"),
                                                 p.get("name") or p.get("id")))
        # Rate-rule states that the SBML declares as parameters (not species).
        for v in extra_state:
            out = [l for l in out if not l.startswith("%-28s :" % self.name(v))]

        reactions = self.m.findall("s:listOfReactions/s:reaction", ns)
        local_lines = []
        for r in reactions:
            kl = r.find("s:kineticLaw", ns)
            for p in kl.findall("s:listOfParameters/s:parameter", ns):
                local_lines.append("%-28s : %-12s : %s, reaction %s" % (
                    r.get("id") + "_" + p.get("id"), p.get("value"), p.get("id"), r.get("name") or r.get("id")))
        out += local_lines

        out += ["", "$CMT @annotated"]
        for s in cmts:
            out.append("%-16s : %s (nM)" % (self.name(s.get("id")), s.get("name") or s.get("id")))
        for v in extra_state:
            out.append("%-16s : %s" % (self.name(v), v))

        out += ["", "$MAIN"]
        for s in cmts:
            sid = s.get("id")
            if sid in init_assign:
                val = self.expr(init_assign[sid].find("m:math", ns))
            else:
                val = s.get("initialConcentration") or s.get("initialAmount") or "0"
            out.append("%s_0 = %s;" % (self.name(sid), val))

        out += ["", "$ODE"]
        for r in assign_rules:
            out.append("double %s = %s;" % (self.name(r.get("variable")),
                                           self.expr(r.find("m:math", ns))))
        flux = {s.get("id"): [] for s in cmts}
        for r in reactions:
            kl = r.find("s:kineticLaw", ns)
            local = {p.get("id"): r.get("id") + "_" + p.get("id")
                     for p in kl.findall("s:listOfParameters/s:parameter", ns)}
            rid = "v_" + r.get("id").rstrip("_")
            out.append("double %s = %s;" % (rid, self.expr(kl.find("m:math", ns), local)))
            for x in r.findall("s:listOfReactants/s:speciesReference", ns):
                flux[x.get("species")].append("- %s%s" % (x.get("stoichiometry", "1") + " * " if x.get("stoichiometry", "1") not in ("1", "1.0") else "", rid))
            for x in r.findall("s:listOfProducts/s:speciesReference", ns):
                flux[x.get("species")].append("+ %s%s" % (x.get("stoichiometry", "1") + " * " if x.get("stoichiometry", "1") not in ("1", "1.0") else "", rid))
        out.append("")
        ruled = {r.get("variable") for r in rate_rules}
        for s in cmts:
            if s.get("id") in ruled:
                continue
            terms = flux[s.get("id")]
            rhs = " ".join(terms).lstrip("+ ").strip() if terms else "0"
            if rhs.startswith("- "):
                rhs = "-" + rhs[2:]
            out.append("dxdt_%s = %s;" % (self.name(s.get("id")), rhs))
        for r in rate_rules:
            out.append("dxdt_%s = %s;" % (self.name(r.get("variable")), self.expr(r.find("m:math", ns))))
        return "\n".join(out) + "\n"


if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    rename = {}
    for a in sys.argv[3:]:
        if a.startswith("--rename="):
            rename = dict(kv.split("=") for kv in a.split("=", 1)[1].split(","))
    c = Converter(src, rename)
    open(dst, "w").write(c.convert("Converted from " + src, []))
