from prettytable import PrettyTable


def profile(fn:str):
    ret = 0
    with open(fn, 'r') as f:
        for line in f:
            if 'Seeding cost' in line:
                ret = float(line.strip().split()[2])
    return ret


def run():
    dataset = [
        'ecoli_SRR1562082_1',
        'elegans_SRR16905161',
        'chicken_SRR13537343_1',
        'human_ERP001775_s1_1',
        'human_ERP001775_s2_1',
        'human_ERP001775_s3_1',
        'human_ERP001775_s4_1',
        'human_ERP001775_s5_1',
        'human_ERR194146_1',
        'human_ERR194161_1',
        'human_ERR3239279_1',
        'human_SRR10965089_1'
    ]
    compressor = ['spring', 'minicom', 'pgrc']
    path_base = '/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/'
    header = ['Dataset', 'BWA', 'Spring', 'Minicom', 'PgRC']
    table = PrettyTable(header)
    for h in header:
        table.align[h] = 'r'
    table.align["Dataset"] = 'l'

    for d in dataset:
        row = [d]
        bwa = 0
        for c in compressor:
            bwa = max(bwa, profile('%s/%s/seed_log/%s.bwa' % (path_base, c, d)))
        row.append(str(bwa))
        for c in compressor:
            comp = profile('%s/%s/seed_log/%s.comp' % (path_base, c, d))
            row.append('%.0f (%.0f)' % (comp, 100 * comp / bwa))
        table.add_row(row)
        print('Dataset %s done' % d)
    print(table)


if __name__ == '__main__':
    run()
