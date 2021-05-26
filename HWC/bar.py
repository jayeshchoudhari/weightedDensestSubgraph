"""
========
Barchart
========

A bar plot with errorbars and height labels on individual bars
"""
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib as mpl


mpl.rcParams.update({'font.size': 32})


# dataset = {
#     'DBLP' : './save_result/dblp_{}_addremove_r4_d10',
#     'CiteULike' : './save_result/citeulike_{}_addremove_r4_d20',
#     'YouTube' : './save_result/youtube_{}_addremove_r4_d200',
# }


dataset = {
    'DBLP' : './save_result/dblp_{}_add_r4',
    'CiteULike' : './save_result/citeulike_{}_add_r4',
    'YouTube' : './save_result/youtube_{}_add_r4',
}


eps_map = {
    0.01 : '0_01',
    0.1 : '0_1',
    0.2 : '0_2',
    0.5 : '0_5',
}


eps = [0.01, 0.1, 0.2, 0.5]
# eps = [0.1, 0.2, 0.5]

N = len(eps)

def process(filename):
    ds, ts = [], []
    with open(filename) as f:
        f.readline()
        for line in f.readlines():
            l = line.split(',')
            density,  time = map(float, [l[2], l[-3]])
            ds.append(density)
            ts.append(time)    
        ts = ts[:-1]
    return np.average(ts), max(ds), np.average(ds)



result = {}

for e in eps:
    for k, v in dataset.items():
        name = v.format(eps_map[e])
        print name
        # avg_time, max_dens, avg_dens = process(name)
        # print avg_time, max_dens, avg_dens
        result[(e, k)] = process(name)

interest = 'avg_time'

i2pos = {
    'avg_time' : 0, 
    'max_dens' : 1, 
    'avg_dens':2
}

ylegend = {
    'avg_time' : "Microseconds", 
    'max_dens' : "Density", 
    'avg_dens': "Density"
}


pos = i2pos[interest]


fr = defaultdict(list)

alldata = ['DBLP', 'CiteULike', 'YouTube']

for k in alldata:
    for e in eps:
        fr[e].append(result[(e, k)][pos])



ind = np.arange(len(alldata))  # the x locations for the groups
width = 0.15      # the width of the bars

fig, ax = plt.subplots()

rect = []
cs = "rbgk"
for i in range(N):
    rect.append(ax.bar(ind+i*width, fr[eps[i]], width, color = cs[i]))


# rects1 = ax.bar(ind, men_means, width, color='r', yerr=men_std)

# women_means = (25, 32, 34, 20, 25)
# women_std = (3, 5, 2, 3, 3)
# rects2 = ax.bar(ind + width, women_means, width, color='y', yerr=women_std)

# add some text for labels, title and axes ticks
ax.set_ylabel( ylegend[interest])
# ax.set_title('Scores by group and gender')
ax.set_xticks(ind + (len(alldata))* width / 2)
ax.set_xticklabels(tuple(alldata))

# ax.set_ybound(lower=0, upper=200)
# ax.legend(   tuple(map(lambda x : x[0], rect)),    tuple(map(str, eps))  ,  loc=2, prop={'size':26})
ax.legend(   tuple(map(lambda x : x[0], rect)),    tuple(map(str, eps))  , prop={'size':26})

# ax.legend((rects1[0], rects2[0]), ('Men', 'Women'))


def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        # ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
        #         '%d' % int(height),
        #         ha='center', va='bottom')
for r in rect:
    autolabel(r)

plt.show()