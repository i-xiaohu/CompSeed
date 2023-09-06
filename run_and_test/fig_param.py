import matplotlib.pyplot as plt


def fetch_time(fn: str):
    ret = 0
    with open(fn, 'r') as f:
        for line in f:
            if 'Seeding cost' in line:
                ret = float(line.strip().split()[2])
    return ret


def run():
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 14}
    light_red = "#C5645D"
    light_blue = '#4779AC'
    light_green = '#80A87C'
    light_yellow = '#D3A428'

    seed_param = [0.5, 1.0, 1.5, 2.0, 2.5][::-1]
    bwa_mem = [140339, 121775, 115652, 115579, 111624]
    spring = [56758, 56300, 54559, 53723, 52996]
    minicom = [61731, 63697, 59366, 59395, 58881]
    pgrc = [58039, 58270, 55469, 55551, 53178]

    plt.figure(figsize=(5, 6), dpi=450)
    plt.subplots_adjust(left=0.20, right=0.90, bottom=0.1, top=0.95, wspace=0.3)
    ax = plt.subplot()
    ax.plot([x for x in seed_param], bwa_mem, marker='s', markersize=10, markerfacecolor='none', linewidth=1, color=light_red, label='BWA-MEM')
    ax.plot([x for x in seed_param], spring,  marker='^', markersize=10, markerfacecolor='none', linewidth=1, color=light_blue, label='CS(SPRING)')
    ax.plot([x for x in seed_param], minicom, marker='o', markersize=10, markerfacecolor='none', linewidth=1, color=light_green, label='CS(Minicom)')
    ax.plot([x for x in seed_param], pgrc,    marker='+', markersize=10, markerfacecolor='none', linewidth=1, color=light_yellow, label='CS(PgRC)')
    ax.set_xlabel("Re-seeding", fontdict=font2)
    ax.set_ylabel("Runtime(s)", fontdict=font2)
    ax.set_xticks(seed_param)
    ax.set_xticklabels(seed_param[::-1])
    ax.tick_params(labelsize=15)
    tick_labels = ax.get_xticklabels() + ax.get_yticklabels()
    [tl.set_fontname('Times New Roman') for tl in tick_labels]
    ax.grid(True, axis='y', linestyle='-', color='gray', linewidth=0.3)
    ax.legend(prop=font2)

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
