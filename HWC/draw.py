import matplotlib as mpl
import sys
import numpy as np
import matplotlib.pyplot as plt


file_name1 = sys.argv[1]

def read_file(file_name):
    ds, sz = [], []
    with open(file_name) as f:
        f.readline()
        for line in f.readlines():
            # print line
            density, size = map(float, line.split(',')[2:4])
            # print density, size
            # r.append((density, size))
            ds.append(density)
            sz.append(size)
    return ds, sz

ds1, sz1 = read_file(file_name1)


fig, ax1 = plt.subplots()

x = range(len(ds1))
line_ds1, = ax1.plot(x, ds1, 'b-', label=file_name1+' density')
ax1.set_xlabel('updates/200')

# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('density', color='b')
ax1.tick_params('y', colors='b')
# ax1.set_ylim([0,15])


ax2 = ax1.twinx()
line_sz1, = ax2.plot(x, sz1, 'r-', label=file_name1+' size')

ax2.set_ylabel('size', color='r')
ax2.tick_params('y', colors='r')
# ax2.set_ylim([0,1500])
fig.tight_layout()

tmp_hand = []
if len(sys.argv) > 2:
    file_name2 = sys.argv[2]
    ds2, sz2 = read_file(file_name2)
    line_ds2, = ax1.plot(x, ds2, 'g-', label=file_name2+' density')
    line_sz2, = ax2.plot(x, sz2, 'k-', label=file_name2+' size')
    tmp_hand += [line_ds2,line_sz2]

plt.legend(handles= [ line_ds1, line_sz1] + tmp_hand)
plt.show()