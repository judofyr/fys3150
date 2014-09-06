from matplotlib.pylab import *
import json

with open("build/data.json") as f:
    data = json.loads(f.read())

def exact(x):
    return 1 - (1 - exp(-10))*x - exp(-10*x)

for run in data:
    if run['algorithm'] == 'LU':
        continue
    if run['n'] > 1000:
        continue

    xs = linspace(0, 1, run['n'])
    ys = run['values']
    plot(xs, ys, label="N=%d" % run['n'])

xs = linspace(0, 1, 10000)
plot(xs, exact(xs), label="Exact")
legend()
savefig("build/values.png")


