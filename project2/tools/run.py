from matplotlib.pylab import *
from subprocess import Popen, PIPE
import sys

def float_line(line):
    return [float(x) for x in line.split()]

class System:
    def __init__(self, n_step, rho_max, w_r=None, algo='jacobi'):
        self.n_step = n_step
        self.rho_max = rho_max
        self.w_r = w_r
        self.x = linspace(0, rho_max, n_step)

        args = ["./build/schrodinger", algo, self.n_step, self.rho_max]
        if not self.w_r is None:
            args.append(self.w_r)
        args = map(str, args)
        self.popen = Popen(args, stdout=PIPE, stderr=PIPE)

    def wait(self):
        out, err = self.popen.communicate()
        self.states = []
        for line in out.splitlines():
            data = float_line(line)
            vec = zeros(self.n_step)
            val = data[0]
            vec[1:] = data[1:]
            n = sqrt(trapz(vec**2, self.x))
            self.states.append((val, vec/n))

        self.states.sort(key=lambda (val, vec): val)
        info = {}
        for line in err.splitlines():
            key, value = line.split('=')
            info[key] = value

        self.num_iter = 1
        if 'iterations' in info:
            self.num_iter = int(info['iterations'])
        if 'elapsed' in info:
            self.elapsed = int(info['elapsed'])

def plot_w_r():
    w_rs = [None, 0.01, 0.5, 1, 5]
    systems = [System(200, 7, w_r) for w_r in w_rs]

    for s in systems:
        s.wait()

    for i in range(3):
        title("lambda_%d" % (i + 1))
        for s in systems:
            if s.w_r is None:
                l = "Single electron"
            else:
                l = "w_r=%g" % s.w_r
            plot(s.x, s.states[i][1]**2, label=l)
        legend()
        savefig("build/l%d.png" % i)
        figure()

steps = [10,50,100,150,200,250,300]

def num_iter():
    systems = [System(n_step, 7) for n_step in steps]
    for s in systems:
        s.wait()
    for s in systems:
        print "%d & %d \\\\" % (s.n_step, s.num_iter)

def prec():
    s = System(280, 7)
    s.wait()
    exact = [3, 7, 11]
    for i in range(3):
        print "%d & %.5f \\\\" % (exact[i], s.states[i][0])

def SI_unit(t):
    digits = log10(t)
    if digits >= 9:
        return "\SI{%.2f}{\second}" % (t / 1e9)
    if digits >= 6:
        return "\SI{%.2f}{\milli\second}" % (t / 1e6)
    if digits >= 3:
        return "\SI{%.2f}{\micro\second}" % (t / 1e3)
    return "\SI{%.2f}{\\nano\second}" % t

def perf():
    systems = [[System(n_step, 7), System(n_step, 7, algo='arma')] for n_step in steps]
    for s1, s2 in systems:
        s1.wait()
        s2.wait()
        print "%d & %s & %s \\\\" % (s1.n_step, SI_unit(s1.elapsed), SI_unit(s2.elapsed))

for arg in sys.argv[1:]:
    eval(arg)()

