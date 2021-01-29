    import matplotlib as mpl
import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

file_name1 = sys.argv[1]

mpl.rcParams.update({'font.size': 32})

def read_file(file_name):
    ds, sz = defaultdict(list), defaultdict(list)
    with open(file_name) as f:
        f.readline()
        for line in f.readlines():
            l = line.split(',')
            density, size = map(float, l[2:4])
            year = int(l[1])
            ds[year].append(density)
            sz[year].append(size)
    x, d, s = [], [], []
    for year in sorted(ds.keys()):
        delta = 1.0 / (len(ds[year]) + 1)
        for i in range(len(ds[year])):
            x.append(year + delta*i)
            d.append(ds[year][i])
            s.append(sz[year][i])
    return x, d, s

x1, d1, s1 = read_file(file_name1)


fig, ax1 = plt.subplots()


line_ds1, = ax1.plot(x1, d1, 'r-', label='Dynamic', linewidth=3)
ax1.set_xlabel('Time')
ax1.set_xlim(xmin=min(x1), xmax=2018)
ax1.set_xbound(lower=min(x1), upper=2018)
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('Density')
ax1.tick_params('y')
# ax1.set_ylim([0,15])

# ax2 = ax1.twinx()

# ax2.set_xbound(lower=min(x1), upper=2018)

# # line_sz1, = ax2.plot(x1, s1, 'r-', label='size')


# ax2.set_ylabel('size', color='r')
# ax2.tick_params('y', colors='r')
# ax2.set_ylim([0,1500])
fig.tight_layout()

tmp_hand = []
if len(sys.argv) > 2:
    file_name2 = sys.argv[2]    
    x2, d2, s2 = read_file(file_name2)
    line_ds2, = ax1.plot(x2, d2, 'b--', label='ELS', linewidth=3)
    # line_sz2, = ax2.plot(x2, s2, 'k-', label=' s ize')
    # tmp_hand += [line_ds2,line_sz2]
    tmp_hand += [line_ds2]

# plt.legend(handles= [ line_ds1, line_sz1] + tmp_hand)
plt.legend(handles= [ line_ds1] + tmp_hand, loc=2, prop={'size':32})
plt.show()