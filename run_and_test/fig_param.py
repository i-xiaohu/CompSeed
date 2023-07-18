import matplotlib.pyplot as plt


def fetch_time(fn: str):
    ret = 0
    with open(fn, 'r') as f:
        for line in f:
            if 'Seeding cost' in line:
                ret = float(line.strip().split()[2])
    return ret


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

    seed_param = [0.5, 1.0, 1.5, 2.0, 2.5][::-1]
    bwa_mem = [1866.11, 1619.26, 1537.84, 1536.88, 1484.28]
    spring = [834.19, 827.46, 801.87, 789.59, 778.9]
    minicom = [827.96, 854.33, 796.23, 796.63, 789.73]
    pgrc = [835.08, 838.4, 798.09, 799.28, 765.14]

    plt.figure(figsize=(5, 6), dpi=450)
    plt.subplots_adjust(left=0.20, right=0.90, bottom=0.1, top=0.95, wspace=0.3)
    ax = plt.subplot()
    ax.plot([x for x in seed_param], bwa_mem, marker='s', markersize=10, markerfacecolor='none', linewidth=1, color=light_red, label='BWA-MEM')
    ax.plot([x for x in seed_param], spring,  marker='^', markersize=10, markerfacecolor='none', linewidth=1, color=light_blue, label='SPRING')
    ax.plot([x for x in seed_param], minicom, marker='o', markersize=10, markerfacecolor='none', linewidth=1, color=light_green, label='Minicom')
    ax.plot([x for x in seed_param], pgrc,    marker='+', markersize=10, markerfacecolor='none', linewidth=1, color=light_yellow, label='PgRC')
    ax.set_xlabel("Re-seeding", fontdict=font2)
    ax.set_ylabel("Runtime(s)", fontdict=font2)
    ax.set_xticks(seed_param)
    ax.set_xticklabels(seed_param[::-1])
    ax.tick_params(labelsize=15)
    tick_labels = ax.get_xticklabels() + ax.get_yticklabels()
    [tl.set_fontname('Times New Roman') for tl in tick_labels]
    ax.grid(True, axis='y', linestyle='-', color='gray', linewidth=0.3)

    plt.savefig('./image_param.png')


def output_time():
    base_path = '/vol1/agis/ruanjue_group/jifahu/zsmem-experiments'
    bwa_mem = []
    for r in ['05', '10', '15', '20', '25']:
        _max = 0
        for c in ['spring', 'minicom', 'pgrc']:
            fn = '%s/%s/seed_log/elegans.r%s.bwa' % (base_path, c, r)
            _max = max(_max, fetch_time(fn))
        bwa_mem.append(_max)
    print('BWA-MEM', bwa_mem)

    for c in ['spring', 'minicom', 'pgrc']:
        comp_list = []
        for r in ['05', '10', '15', '20', '25']:
            fn = '%s/%s/seed_log/elegans.k%s.comp' % (base_path, c, r)
            comp_list.append(fetch_time(fn))
        print(c, comp_list)


if __name__ == '__main__':
    run()
    # output_time()
