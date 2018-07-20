
def append_counts(path):
    with open(path, 'r') as infile:
        header = next(infile).split('\t')
        header_map = {i:j for i,j in zip(header, range(1,len(header) + 1))}
        for line in infile:
            line = line.split('\t')
            mir = line[header_map['mir']]
            mrna = line[header_map['bound']]
            counts = line[header_map['counts']]
            bound = "%s_%s" %(mir, mrna)
            dd[bound].append(int(counts))

def average_counts():
    for path in snu_files:
        append_counts(path)
    avg = {i:np.sum(dd[i])/3 for i in dd}
    return avg
avg = average_counts()