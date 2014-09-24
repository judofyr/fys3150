from matplotlib.pylab import *

files = [
    ('a', 'w_r = 0.01'),
    ('b', 'w_r = 1'),
    ('c', 'w_r = 5')
]

for x, title in files:
    filename = "build/wr%s.txt" % x
    data = [float(x) for x in open(filename)]
    plot(range(len(data)), data, label=title)

legend()
show()

