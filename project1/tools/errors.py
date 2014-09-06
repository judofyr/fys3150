from matplotlib.pylab import *
import json

with open("build/data.json") as f:
    data = json.loads(f.read())

def exact(x):
    return 1 - (1 - exp(-10))*x - exp(-10*x)

for run in data:
    if run['algorithm'] == 'LU':
        continue

    xs = linspace(0, 1, run['n'])
    us = exact(xs)
    vs = run['values']
    eps = log10(abs((vs-us)/us))
    plot(xs, eps, label="N=%d" % run['n'])

legend()
savefig("build/errors.png")



