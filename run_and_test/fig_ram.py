import matplotlib.pyplot as plt
import numpy as np


def run():
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 21}

    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 18}

    light_red = "#e5c9c4"
    light_blue = '#b2d1d4'
    light_green = '#80A87C'
    light_grey = '#DEE0E1'
    light_yellow = '#D3A428'
    big_red = '#e83929'

    fig, ax1 = plt.subplots()
    fig.set_size_inches(6, 6)
    fig.set_dpi(450)
    fig.subplots_adjust(left=0.20, right=0.85, bottom=0.20, top=0.95)
    bwa_time, bwa_ram = 115652, 6
    mem2_time, mem2_ram = 77919, 27
    spring, minicom, pgrc = 55309, 57016, 56198
    comp_time, comp_ram = (spring + minicom + pgrc) / 3, 8
    ert_time, ert_ram = 41977, 74
    pmi_time, pmi_ram = 35295, 120

    ax1.bar(0, bwa_time, color=light_red, align='edge', width=2)
    ax1.bar(3, comp_time, color=light_green, align='edge', width=2)
    ax1.bar(6, mem2_time, color=light_blue, align='edge', width=2)
    ax1.bar(9, ert_time, color=light_grey, align='edge', width=2)
    ax1.bar(12, pmi_time, color=light_yellow, align='edge', width=2)

    ax1.set_ylabel("Runtime(s)", fontdict=font2)
    ax1.set_xticks([0, 3, 6, 10, 13])
    ax1.set_xticklabels(['FM-Index', 'Comp-Seed', 'BWA-MEM2', 'ERT', 'BWA-MEME'], rotation=45)
    ax1.tick_params(axis='x', length=0)
    ax1.tick_params(labelsize=15)
    [tl.set_fontname('Times New Roman') for tl in ax1.get_xticklabels()]
    [tl.set_fontname('Times New Roman') for tl in ax1.get_yticklabels()]
    ax1.grid(True, axis='y', linestyle='--', color='gray', linewidth=0.2)

    ax2 = ax1.twinx()
    ram_base = [1, 4, 7, 10, 13]
    ram_line = [bwa_ram, comp_ram, mem2_ram, ert_ram, pmi_ram]
    ax2.plot(ram_base, ram_line, color=big_red, marker='s', markersize=7, linewidth=2)

    ax2.set_ylabel("RAM(GB)", fontdict=font2, color=big_red)
    ax2.tick_params(axis='y', labelcolor=big_red)
    ax2.tick_params(labelsize=15)
    [tl.set_fontname('Times New Roman') for tl in ax2.get_yticklabels()]

    plt.savefig('image_ram.png')


if __name__ == '__main__':
    run()
