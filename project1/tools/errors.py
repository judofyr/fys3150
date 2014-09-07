from matplotlib.pylab import *
import json

with open("build/data.json") as f:
    data = json.loads(f.read())

def exact(x):
    return 1 - (1 - exp(-10))*x - exp(-10*x)

maxs = {}

for run in data:
    if run['algorithm'] == 'LU':
        continue

    if run['n'] > 1e5:
        continue

    xs = linspace(0, 1, run['n']+2)[1:-1]
    us = exact(xs)
    vs = run['values'][1:-1]
    eps = ma.log10(abs((vs-us)/us))
    plot(xs, eps, label="N=%d" % run['n'])
    maxs[run['n']] = nanmax(eps)

legend()
savefig("build/errors.png")

for n in sort(maxs.keys()):
    print "%d & %s \\\\" % (n, maxs[n])

