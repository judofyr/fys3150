import subprocess
import json
import sys

def run(*args):
    sys.stderr.write("[*] Running: %s\n" % (" ".join(args)))
    child = subprocess.Popen(args, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    (stdout, stderr) = child.communicate()
    values = [0.]
    for x in stdout.splitlines():
        values.append(float(x))
    values.append(0.)
    elapsed_start = stderr.find("elapsed: ") + 9
    elapsed_stop = stderr.find(" ", elapsed_start)
    elapsed = int(stderr[elapsed_start:elapsed_stop])
    return (values, elapsed)

data = []

for alg in ["LU", "Tri"]:
    for i in xrange(1, 7):
        n = 10**i
        if alg == "LU" and i >= 4:
            continue

        values, elapsed = run("./build/tridiagonal", alg, str(n))
        data.append(
            {
                'algorithm': alg,
                'n': n,
                'values': values,
                'elapsed': elapsed
            }
        )

print json.dumps(data)

