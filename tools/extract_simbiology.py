"""Extract reactions, rules, parameters and variants from a SimBiology project.

Usage: python tools/extract_simbiology.py <project.sbproj> <out.json>

A .sbproj file is a zip archive. Its simbiodata.mat holds the model as a
MATLAB handle object (class SimBiology.Model) serialised as one flat cell
array of every property value. This reads that cell array without MATLAB and
recovers, by their storage patterns:
  - reactions: the reaction string followed by its rate expression
  - rules: strings of the form "name = expression"
  - parameters: name, '', '', [], value, units
  - variants: cell arrays of {'parameter', name, 'Value', value}; the variant
    name is stored after its contents
Used to rebuild the T-cell engager model of Hosseini et al. (npj Syst Biol
Appl 2020;6:28) from the project in its Supplementary Software.
"""
import ast
import json
import re
import struct
import sys
import zipfile
import zlib

FMT = {1: 'b', 2: 'B', 3: 'h', 4: 'H', 5: 'i', 6: 'I', 7: 'f', 9: 'd', 12: 'q', 13: 'Q', 16: 'B', 17: 'H', 18: 'I'}
SZ = {'b': 1, 'B': 1, 'h': 2, 'H': 2, 'i': 4, 'I': 4, 'f': 4, 'd': 8, 'q': 8, 'Q': 8}


def decompress(mat):
    out, p = b'', 128
    while p < len(mat):
        t, n = struct.unpack_from('<II', mat, p)
        d = mat[p + 8:p + 8 + n]
        out += zlib.decompress(d) if t == 15 else struct.pack('<II', t, n) + d
        p += 8 + n
    return out


class Reader:
    def __init__(self, b):
        self.b = b

    def tag(self, p):
        t, n = struct.unpack_from('<II', self.b, p)
        if t >> 16:
            return t & 0xffff, t >> 16, p + 4, p + 8
        return t, n, p + 8, p + 8 + ((n + 7) // 8) * 8

    def vals(self, t, n, dp):
        f = FMT.get(t, 'B')
        return list(struct.unpack_from('<%d%s' % (n // SZ[f], f), self.b, dp))

    def children(self, p, end):
        items = []
        while p < end:
            t, n, dp, nxt = self.tag(p)
            items.append(self.matrix(dp, dp + n) if t == 14 else None)
            p = nxt
        return items

    def matrix(self, p, end):
        if end - p < 8:
            return None
        t, n, dp, p = self.tag(p); cls = self.vals(t, n, dp)[0] & 0xff
        t, n, dp, p = self.tag(p)                         # dims
        t, n, dp, p = self.tag(p)                         # name
        if cls == 1:
            return self.children(p, end)
        if cls in (2, 3):
            if cls == 3:
                t, n, dp, p = self.tag(p)
            t, n, dp, p = self.tag(p)
            t, n, dp, p = self.tag(p)
            return ('struct', self.children(p, end))
        if cls == 17:
            return None
        t, n, dp, p = self.tag(p)
        v = self.vals(t, n, dp)
        if cls == 4:
            return ''.join(chr(c) for c in v)
        return v if len(v) != 1 else v[0]


def main(project, out):
    mat = zipfile.ZipFile(project).read('simbiodata.mat')
    raw = decompress(mat)
    r = Reader(raw)
    # the handle object: 'handle', 'SimBiology.Model', then the cell array
    t, n, dp, nxt = r.tag(72)
    cells = r.matrix(dp, dp + n)

    reactions, rules, params, variants = [], [], {}, []
    ident = re.compile(r'[A-Za-z_][A-Za-z0-9_]*$')
    for i, c in enumerate(cells):
        if isinstance(c, str) and '->' in c and i + 1 < len(cells):
            rate = cells[i + 1]
            if isinstance(rate, str) and not rate.startswith('mw') and rate:
                reactions.append({'reaction': c, 'rate': rate})
        elif isinstance(c, str) and ' = ' in c and '->' not in c:
            rules.append(c)
        elif (isinstance(c, str) and ident.match(c) and i + 5 < len(cells)
              and cells[i + 1] == '' and cells[i + 2] == '' and cells[i + 3] == []
              and isinstance(cells[i + 4], (int, float))):
            params.setdefault(c, cells[i + 4])
        elif isinstance(c, list) and c and all(isinstance(x, list) and len(x) == 4 and x[0] == 'parameter' for x in c):
            variants.append({'values': {x[1]: x[3] for x in c}})
    # a variant's name is the first plain string after its contents
    k = 0
    for i, c in enumerate(cells):
        if isinstance(c, list) and c and all(isinstance(x, list) and len(x) == 4 and x[0] == 'parameter' for x in c):
            nm = next(x for x in cells[i + 1:i + 12] if isinstance(x, str) and x and not x.startswith(('mw', '<')) and '-' not in x[8:9])
            variants[k]['name'] = nm
            k += 1
    json.dump({'reactions': reactions, 'rules': rules, 'parameters': params, 'variants': variants},
              open(out, 'w'), indent=1)
    print(len(reactions), 'reactions,', len(rules), 'rules,', len(params), 'parameters,', len(variants), 'variants')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
