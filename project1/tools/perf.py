# coding: utf-8

from matplotlib.pylab import *
import json

with open("build/data.json") as f:
    data = json.loads(f.read())

def SI_unit(t):
    digits = log10(t)
    if digits >= 9:
        return "\SI{%.2f}{\second}" % (t / 1e9)
    if digits >= 6:
        return "\SI{%.2f}{\milli\second}" % (t / 1e6)
    if digits >= 3:
        return "\SI{%.2f}{\micro\second}" % (t / 1e3)
    return "\SI{%.2f}{\\nano\second}" % t

lines = {}

for run in data:
    n = run['n']
    if not n in lines:
        lines[n] = {}
    lines[n][run['algorithm']] = run['elapsed']

for n in sort(lines.keys()):
    runs = lines[n]

    lu = ""
    tri = ""

    if 'LU' in runs:
        lu = SI_unit(runs['LU'])
    if 'Tri' in runs:
        tri = SI_unit(runs['Tri'])

    print "%d & %s & %s \\\\" % (n, tri, lu)

