from prettytable import PrettyTable


def fetch_time(fn: str):
    ret = 0
    with open(fn, 'r') as f:
        for line in f:
            if 'Seeding cost' in line:
                ret = float(line.strip().split()[2])
    return ret


def profile_time(dataset):
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
            bwa = max(bwa, fetch_time('%s/%s/seed_log/%s.bwa' % (path_base, c, d)))
        row.append('%.0f' % bwa)
        for c in compressor:
            comp = fetch_time('%s/%s/seed_log/%s.comp' % (path_base, c, d))
            row.append('%.0f (%.0f)' % (comp, 100 * comp / bwa))
        table.add_row(row)
    print(table)


def fetch_bwt(fn: str):
    ret = 0
    with open(fn, 'r') as f:
        for line in f:
            if 'BWT-extend:' in line:
                ret = int(line.strip().split()[1])
    return ret


def profile_bwt(dataset):
    compressor = ['spring', 'minicom', 'pgrc']
    path_base = '/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/'
    header = ['Dataset', 'BWA', 'Spring', 'Minicom', 'PgRC']
    table = PrettyTable(header)
    for h in header:
        table.align[h] = 'r'
    table.align["Dataset"] = 'l'

    for d in dataset:
        row = [d]
        bwa = -1
        for c in compressor:
            if bwa == -1:
                bwa = fetch_bwt('%s/%s/seed_log/%s.bwa' % (path_base, c, d))
            else:
                assert bwa == fetch_bwt('%s/%s/seed_log/%s.bwa' % (path_base, c, d))
        row.append(str(bwa))
        for c in compressor:
            comp = fetch_bwt('%s/%s/seed_log/%s.comp' % (path_base, c, d))
            row.append('%.0f (%.0f)' % (comp, 100 * comp / bwa))
        table.add_row(row)
    print(table)


def fetch_sal(fn: str):
    ret = 0
    with open(fn, 'r') as f:
        for line in f:
            if 'SA Lookup:' in line:
                ret = int(line.strip().split()[2])
    return ret


def profile_sal(dataset):
    compressor = ['spring', 'minicom', 'pgrc']
    path_base = '/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/'
    header = ['Dataset', 'BWA', 'Spring', 'Minicom', 'PgRC']
    table = PrettyTable(header)
    for h in header:
        table.align[h] = 'r'
    table.align["Dataset"] = 'l'

    for d in dataset:
        row = [d]
        bwa = -1
        for c in compressor:
            if bwa == -1:
                bwa = fetch_sal('%s/%s/seed_log/%s.bwa' % (path_base, c, d))
            else:
                assert bwa == fetch_sal('%s/%s/seed_log/%s.bwa' % (path_base, c, d))
        row.append(str(bwa))
        for c in compressor:
            comp = fetch_sal('%s/%s/seed_log/%s.comp' % (path_base, c, d))
            row.append('%.0f (%.0f)' % (comp, 100 * comp / bwa))
        table.add_row(row)
    print(table)


if __name__ == '__main__':
    all_dataset = [
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
    profile_time(all_dataset)
    profile_bwt(all_dataset)
    profile_sal(all_dataset)

