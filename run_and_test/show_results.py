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

    for (d, n) in dataset:
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


def fetch_sal(fn: str):
    ret = 0
    with open(fn, 'r') as f:
        for line in f:
            if 'SA Lookup:' in line:
                ret = int(line.strip().split()[2])
    return ret


def profile_bwt_sal(dataset):
    compressor = ['spring', 'minicom', 'pgrc']
    path_base = '/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/'
    header = ['Dataset', 'BWT', 'Spring', 'Minicom', 'PgRC', "SAL", "SPRING", "MINICOM", "PGRC"]
    table = PrettyTable(header)
    for h in header:
        table.align[h] = 'r'
    table.align["Dataset"] = 'l'

    for (d, n) in dataset:
        row = [d]
        bwa = -1
        for c in compressor:
            if bwa == -1:
                bwa = fetch_bwt('%s/%s/seed_log/%s.bwa' % (path_base, c, d))
            else:
                assert bwa == fetch_bwt('%s/%s/seed_log/%s.bwa' % (path_base, c, d))
        bwa /= n
        row.append('%.0f' % bwa)
        for c in compressor:
            comp = fetch_bwt('%s/%s/seed_log/%s.comp' % (path_base, c, d))
            comp /= n
            row.append('%.0f (%.0f)' % (comp, 100 * comp / bwa))

        bwa = -1
        for c in compressor:
            if bwa == -1:
                bwa = fetch_sal('%s/%s/seed_log/%s.bwa' % (path_base, c, d))
            else:
                assert bwa == fetch_sal('%s/%s/seed_log/%s.bwa' % (path_base, c, d))
        bwa /= n
        row.append('%.2f' % bwa)
        for c in compressor:
            comp = fetch_sal('%s/%s/seed_log/%s.comp' % (path_base, c, d))
            comp /= n
            row.append('%.2f (%.0f)' % (comp, 100 * comp / bwa))

        table.add_row(row)
    print(table)


if __name__ == '__main__':
    all_dataset = [
        ['ecoli_SRR1562082_1', 5825771],
        ['elegans_SRR16905161', 15046618],
        ['chicken_SRR13537343_1', 201583321],
        ['human_ERP001775_s1_1', 223571196],
        ['human_ERP001775_s2_1', 439498957],
        ['human_ERP001775_s3_1', 652258345],
        ['human_ERP001775_s4_1', 863941088],
        ['human_ERP001775_s5_1', 1075379007],
        ['human_ERR194146_1', 813180578],
        ['human_ERR194161_1', 843454257],
        ['human_ERR3239279_1', 420210145],
        ['human_SRR10965089_1', 376183716]
    ]
    profile_time(all_dataset)
    profile_bwt_sal(all_dataset)

