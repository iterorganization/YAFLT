
x = []
y = []
xlim = []
ylim = []
with open("test_1_accuracy.out", "r") as f:
    N = 0
    for line in f:
        if "Comparing" in line:
            N = int(line.split()[1])
            break
    for line in f:
        if line.startswith('Point'):
            sline = line.split()
            x.append(float(sline[1]))
            y.append(float(sline[2]))
        elif line.startswith('R goes from'):
            sline = line.split()
            xlim = [float(sline[3]), float(sline[5])]
        elif line.startswith('Z goes from'):
            sline = line.split()
            ylim = [float(sline[3]), float(sline[5])]

import matplotlib.pyplot as plt

plt.plot(x, y, '.')
plt.xlim(*xlim)
plt.ylim(*ylim)
plt.axis("equal")
plt.show()
