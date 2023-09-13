import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from matplotlib.pyplot import MultipleLocator


def run():
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 21}

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
    spring =   [24981, 38310, 46399, 55309,  61367]
    minicom =  [25747, 40943, 47775, 57016,  62239]
    pgrc =     [24492, 38225, 46663, 56198,  60923]

    plt.figure(figsize=(20, 8), dpi=450)
    plt.subplots_adjust(left=0.06, right=0.99, bottom=0.1, top=0.95, wspace=0.3)
    ax = plt.subplot(131)
    ax.plot([x for x in base_x], bwa_mem, marker='s', markersize=10, markerfacecolor='none', linewidth=1, color=light_red, label='BWA-MEM')
    ax.plot([x for x in base_x], spring,  marker='^', markersize=10, markerfacecolor='none', linewidth=1, color=light_blue, label='CS(SPRING)')
    ax.plot([x for x in base_x], minicom, marker='o', markersize=10, markerfacecolor='none', linewidth=1, color=light_green, label='CS(Minicom)')
    ax.plot([x for x in base_x], pgrc,    marker='+', markersize=10, markerfacecolor='none', linewidth=1, color=light_yellow, label='CS(PgRC)')
    ax.set_xlabel("Coverage", fontdict=font2)
    ax.set_ylabel("Seeding time(s)", fontdict=font2)
    ax.set_xticks(base_x)
    ax.set_xticklabels(coverage)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(labelsize=15)
    tick_labels = ax.get_xticklabels() + ax.get_yticklabels()
    [tl.set_fontname('Times New Roman') for tl in tick_labels]
    ax.legend(prop=font2)
    ax.grid(True, axis='y', linestyle='-', color='gray', linewidth=0.3)
    ax.set_title('(a)', loc='left', fontdict=font1, pad=10)

    # Figure b
    bwa_mem = [624, 623, 622, 622, 622]
    spring =  [335, 287, 256, 233, 216]
    minicom = [339, 298, 270, 248, 218]
    pgrc =    [346, 292, 257, 233, 214]
    ax = plt.subplot(132)
    ax.plot([x for x in base_x], bwa_mem, marker='s', markersize=10, markerfacecolor='none', linewidth=1, color=light_red, label='BWA-MEM')
    ax.plot([x for x in base_x], spring,  marker='^', markersize=10, markerfacecolor='none', linewidth=1, color=light_blue, label='CS(SPRING)')
    ax.plot([x for x in base_x], minicom, marker='o', markersize=10, markerfacecolor='none', linewidth=1, color=light_green, label='CS(Minicom)')
    ax.plot([x for x in base_x], pgrc,    marker='+', markersize=10, markerfacecolor='none', linewidth=1, color=light_yellow, label='CS(PgRC)')
    ax.set_xlabel("Coverage", fontdict=font2)
    ax.set_ylabel("BWT-extend/Read", fontdict=font2)
    ax.set_xticks(base_x)
    ax.set_xticklabels(coverage)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(labelsize=15)
    tick_labels = ax.get_xticklabels() + ax.get_yticklabels()
    [tl.set_fontname('Times New Roman') for tl in tick_labels]
    ax.grid(True, axis='y', linestyle='-', color='gray', linewidth=0.3)
    ax.set_title('(b)', loc='left', fontdict=font1, pad=10)

    # Figure c
    bwa_mem = [29.11, 29.06, 29.02, 29.00, 28.99]
    spring =  [15.22, 12.96, 11.54, 10.54, 9.80]
    minicom = [15.81, 13.75, 12.43, 11.48, 10.72]
    pgrc =    [15.78, 13.28, 11.76, 10.72, 9.93]
    ax = plt.subplot(133)
    ax.plot([x for x in base_x], bwa_mem, marker='s', markersize=10, markerfacecolor='none', linewidth=1, color=light_red, label='BWA-MEM')
    ax.plot([x for x in base_x], spring,  marker='^', markersize=10, markerfacecolor='none', linewidth=1, color=light_blue, label='CS(SPRING)')
    ax.plot([x for x in base_x], minicom, marker='o', markersize=10, markerfacecolor='none', linewidth=1, color=light_green, label='CS(Minicom)')
    ax.plot([x for x in base_x], pgrc,    marker='+', markersize=10, markerfacecolor='none', linewidth=1, color=light_yellow, label='CS(PgRC)')
    ax.set_xlabel("Coverage", fontdict=font2)
    ax.set_ylabel("SAL/Read", fontdict=font2)
    ax.set_xticks(base_x)
    ax.set_xticklabels(coverage)
    ax.tick_params(axis='x', length=0)
    ax.tick_params(labelsize=15)
    tick_labels = ax.get_xticklabels() + ax.get_yticklabels()
    [tl.set_fontname('Times New Roman') for tl in tick_labels]
    ax.grid(True, axis='y', linestyle='-', color='gray', linewidth=0.3)
    ax.set_title('(c)', loc='left', fontdict=font1, pad=10)

    plt.savefig("./Figure 3.png")


if __name__ == '__main__':
    run()
