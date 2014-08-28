from matplotlib.pylab import *
import json

with open("data.json") as f:
    data = json.loads(f.read())

def exact(x):
    return 1 - (1 - exp(-10))*x - exp(-10*x)

for run in data:
    if run['algorithm'] == 'LU':
        continue

    xs = linspace(0, 1, run['n'])
    us = exact(xs)
    vs = run['values']
    max_eps = float('-inf')
    max_u = 0
    for (u, v) in zip(us, vs):
        if u == 0:
            continue
        eps = log10(abs((v - u)/u))
        if eps > max_eps:
            max_eps = eps
            max_u = u
    print "N=%d eps=%f at %f" % (run['n'], max_eps, max_u)




