import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from matplotlib.pyplot import MultipleLocator


def run():
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 20}

    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 18}

    light_red = "#C5645D"
    light_blue = '#4779AC'
    light_green = '#80A87C'
    light_yellow = '#D3A428'

    # Figure a
    coverage = ['7X',  '14X', '21X', '28X',  '35X']
    base_x =   [0,     10,    20,    30,     40]
    bwa_mem =  [31644, 59391, 88529, 115652, 146973]
    spring =   [24981, 38310, 46399, 50559,  60591]
    minicom =  [25747, 40943, 47775, 51216,  62403]
    pgrc =     [23492, 36725, 45163, 49366,  59632]

    plt.figure()
    ax = plt.subplot(131)
    ax.bar([x + 0 for x in base_x], bwa_mem, align='edge', width=1.5, color=light_red, label='BWA-MEM')
    ax.bar([x + 2.5 for x in base_x], spring,  align='edge', width=1.5, color=light_blue, label='SPRING')
    ax.bar([x + 4.5 for x in base_x], minicom, align='edge', width=1.5, color=light_green, label='Minicom')
    ax.bar([x + 6.5 for x in base_x], pgrc,    align='edge', width=1.5, color=light_yellow, label='PgRC')
    ax.set_xlabel("Coverage", fontdict=font2)
    ax.set_ylabel("Runtime(s)", fontdict=font2)
    ax.set_xticks([x + 3.5 for x in base_x])
    ax.set_xticklabels(coverage)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(labelsize=15)
    tick_labels = ax.get_xticklabels() + ax.get_yticklabels()
    [tl.set_fontname('Times New Roman') for tl in tick_labels]
    ax.legend(prop=font1)

    # Figure B
    # bwa_mem = [22580690796, 44389394657, 65878092845, 87258049888, 108613279707]
    # bwa_mem = [b / 1000000000 for b in bwa_mem]
    # spring =  [1046353920, 1269032960, 1500477440, 1724917760, 1939814400]
    # spring =  [s / 1000000000 for s in spring]
    # minicom = [1533306880, 1840066560, 2106398720, 2356889600, 2591240960]
    # minicom = [m / 1000000000 for m in minicom]
    # pgrc =    [1125829786, 1206283685, 1358332405, 1518122050, 1672546991]
    # pgrc =    [p / 1000000000 for p in pgrc]
    # ax = plt.subplot(222)
    # ax.bar([x + 0 for x in base_x], bwa_mem, align='edge', width=1.5, color=light_red, label='BWA-MEM')
    # ax.bar([x + 2.5 for x in base_x], spring,  align='edge', width=1.5, color=light_blue, label='SPRING')
    # ax.bar([x + 4.5 for x in base_x], minicom, align='edge', width=1.5, color=light_green, label='Minicom')
    # ax.bar([x + 6.5 for x in base_x], pgrc,    align='edge', width=1.5, color=light_yellow, label='PgRC')
    # ax.set_xlabel("Coverage", fontdict=font2)
    # ax.set_ylabel("Size(GB)", fontdict=font2)
    # ax.set_xticks([x + 3.5 for x in base_x])
    # ax.set_xticklabels(coverage)
    # ax.tick_params(axis='x', length=0)
    # ax.tick_params(labelsize=15)
    # tick_labels = ax.get_xticklabels() + ax.get_yticklabels()
    # [tl.set_fontname('Times New Roman') for tl in tick_labels]

    # Figure c
    bwa_mem = [624, 623, 622, 622, 622]
    spring =  [335, 287, 256, 233, 216]
    minicom = [339, 298, 270, 248, 218]
    pgrc =    [346, 292, 257, 233, 214]
    ax = plt.subplot(132)
    ax.bar([x + 0 for x in base_x], bwa_mem, align='edge', width=1.5, color=light_red, label='BWA-MEM')
    ax.bar([x + 2.5 for x in base_x], spring,  align='edge', width=1.5, color=light_blue, label='SPRING')
    ax.bar([x + 4.5 for x in base_x], minicom, align='edge', width=1.5, color=light_green, label='Minicom')
    ax.bar([x + 6.5 for x in base_x], pgrc,    align='edge', width=1.5, color=light_yellow, label='PgRC')
    ax.set_xlabel("Coverage", fontdict=font2)
    ax.set_ylabel("BWT-Extend/Read", fontdict=font2)
    ax.set_xticks([x + 3.5 for x in base_x])
    ax.set_xticklabels(coverage)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(labelsize=15)
    tick_labels = ax.get_xticklabels() + ax.get_yticklabels()
    [tl.set_fontname('Times New Roman') for tl in tick_labels]

    # Figure d
    bwa_mem = [29.11, 29.06, 29.02, 29.00, 28.99]
    spring =  [15.22, 12.96, 11.54, 10.54, 9.80]
    minicom = [15.81, 13.75, 12.43, 11.48, 10.72]
    pgrc =    [15.78, 13.28, 11.76, 10.72, 9.93]
    ax = plt.subplot(133)
    ax.bar([x + 0 for x in base_x], bwa_mem, align='edge', width=1.5, color=light_red, label='BWA-MEM')
    ax.bar([x + 2.5 for x in base_x], spring,  align='edge', width=1.5, color=light_blue, label='SPRING')
    ax.bar([x + 4.5 for x in base_x], minicom, align='edge', width=1.5, color=light_green, label='Minicom')
    ax.bar([x + 6.5 for x in base_x], pgrc,    align='edge', width=1.5, color=light_yellow, label='PgRC')
    ax.set_xlabel("Coverage", fontdict=font2)
    ax.set_ylabel("SAL/Read", fontdict=font2)
    ax.set_xticks([x + 3.5 for x in base_x])
    ax.set_xticklabels(coverage)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(labelsize=15)
    tick_labels = ax.get_xticklabels() + ax.get_yticklabels()
    [tl.set_fontname('Times New Roman') for tl in tick_labels]

    plt.show()


if __name__ == '__main__':
    run()
