import subprocess

path = '/Users/mac9/Desktop/clash_jul2018/180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_3_1_comp_human_ebv_v28encode_mirbase_hybrids_ua.hyb.mirs.protein'
gtf_path = '/Users/mac9/Downloads/gencode.v28.annotation.gtf'
fa_path = '/Users/mac9/Downloads/gencode.v28.transcripts.fa'


def run_rna_fold(mir_seq, mrna_seq):

    folding = subprocess.check_output('echo "%s&%s" | RNAcofold --noPS' % (mir_seq, mrna_seq), shell=True).decode('utf8')
    folding_structure, free_energy = folding.split('\n')[1].split(' ', 1)
    mirna_half, mrna_half = folding_structure.split('&')
    free_energy = free_energy.lstrip('(').rstrip(')')
    
    return mirna_half, mrna_half, free_energy


def get_transcript_seq(fa_path):
    with open (fa_path) as fa:
        fa_dict = {}
        transcript, seq = '', ''
        
        for line in fa:
            if line[0] == '>':
                if seq:
                    fa_dict[transcript] = seq
                seq = ''
                transcript = line[1:].split('|')[0]
            else:
                seq+=line.strip()

    return fa_dict


def get_utrs(gtf_path):
    with open(gtf_path,'r') as gtf:
        dic = {}
        
        for line in gtf:

            # Transcripts with UTR and CDS
            #if 'protein_coding' in line:
            if not line.startswith('#'):
                _, _, feature, start, end, _, strand, _, attributes  = line.split('\t')
                start, end = int(start), int(end)

                if 'transcript_id' in attributes:
                    transcript_id = attributes.split('transcript_id')[1].split('"')[1] # Get rid of version info.. ENSG00000146648.2 > ENSG00000146648

                    if transcript_id not in dic:
                        dic[transcript_id] = {
                            'transcript_length': 0,
                            'strand': strand,
                            'utr_coords': []
                        }

                    if feature == 'exon':
                        exon_len = (end - start) + 1
                        dic[transcript_id]['transcript_length'] += exon_len

                    elif feature == 'start_codon':
                        dic[transcript_id]['start_codon'] = [start, end]

                    elif feature == 'stop_codon':
                        dic[transcript_id]['stop_codon'] = [start, end]

                    elif feature == 'UTR':
                        dic[transcript_id]['utr_coords'].append([start, end])


        for transcript in dic:
            trans = dic[transcript]
            if ('stop_codon' in trans and 'start_codon' in trans):
                start_codon = trans['start_codon']
                stop_codon = trans['stop_codon']
                strand = trans['strand']
                utr_coords = trans['utr_coords']

                if strand == '+':
                    dic[transcript]['UTR3'] = sum((last - first) + 1 for first, last in utr_coords if first >= stop_codon[0])
                    dic[transcript]['UTR5'] = sum((last - first) + 1 for first, last in utr_coords if last < start_codon[0])

                elif strand == '-':
                    dic[transcript]['UTR3'] = sum((last - first) + 1 for first, last in utr_coords if first <= stop_codon[1])
                    dic[transcript]['UTR5'] = sum((last - first) + 1 for first, last in utr_coords if last > start_codon[1])

                dic[transcript]['utr5range'] = (0, dic[transcript]['UTR5'])
                dic[transcript]['utr3range'] = (dic[transcript]['transcript_length'] - dic[transcript]['UTR3'], dic[transcript]['transcript_length'])
    return dic

utr_dict = get_utrs(gtf_path)
fasta_dict = get_transcript_seq(fa_path)

for path in ['/Users/mac9/Desktop/clash_jul2018/180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_3_1_comp_human_ebv_v28encode_mirbase_hybrids_ua.hyb.mirs.protein','/Users/mac9/Desktop/clash_jul2018/180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_1_1_comp_human_ebv_v28encode_mirbase_hybrids_ua.hyb.mirs.protein', '/Users/mac9/Desktop/clash_jul2018/180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_2_1_comp_human_ebv_v28encode_mirbase_hybrids_ua.hyb.mirs.protein']:
    with open(path) as infile, open("%s.annotated" % path, 'w') as outfile:
        
        header = [
            'read_id', 'sequence', '.', 'gene1', 'read1_start', 'read1_stop',
            'gene1_start', '    gene1_stop', 'gene1_align_score', 'gene2', 'read2_start', 
            'read2_stop', 'gene2_start', 'gene2_stop', 'gene2_align_score',
            'mir', 'mrna', 'transcript', 'mir_seq', 'mrna_seq', 'mir_fold',
            'mrna_fold', 'binding_energy', 'a1', 'm2_5' ,'m6_7', 'm8', 'm9_12', 'm13_17',
            'mrna_feature_start', 'mrna_feature_stop', 'count'
            ]
        outfile.write('%s\n' % '\t'.join(header))
        line_count = 0
        total_lines = subprocess.check_output('wc -l %s' % path, shell=True).decode('utf8').split()[0]
        
        print("total lines:", total_lines)
        
        for line in infile:
            line = line.strip().split('\t')
            mir_column, mrna_column = (3,9) if 'miR' in line[3] else (9,3)
            mir_gene = line[mir_column].split('_')[0]  
            transcript = line[mrna_column].split('_')[0]
            mrna_gene = line[mrna_column].split('_')[-4]  # Hugo symbol

            mrna_pos_start, mrna_pos_end = int(line[mrna_column + 3]), int(line[mrna_column + 4])

            mir_start, mir_stop = int(line[mir_column + 1]) - 1, int(line[mir_column + 2]) - 1
            mrna_start, mrna_stop = int(line[mrna_column + 1]) - 1, int(line[mrna_column + 2]) - 1
            sequence = line[1]
            mir_seq, mrna_seq = sequence[mir_start: mir_stop + 1], sequence[mrna_start: mrna_stop + 1] 

            mir, mrna, fe = run_rna_fold(mir_seq=mir_seq, mrna_seq=mrna_seq)
            mir_first_match = mir.find('(')
            mrna_first_match = mrna.rfind(')')
            mrna_start = mrna_first_match + mir_first_match
            count = line[0].split('_')[1]
            a = '0'
            
            if len(mrna_seq) > mrna_start + 1:
                if mrna_seq[mrna_start] == 'A':
                    a = '1'  # If 'A' is the base opposite mir position '1'.
            elif transcript in fasta_dict:
                mrna_align_start = int(line[mrna_column + 3]) - 1
                if len(fasta_dict[transcript]) >= mrna_start + mrna_align_start + 1:
                    fa_a = fasta_dict[transcript][mrna_start + mrna_align_start]
                    if fa_a == 'A':
                        a = '1'
                else:
                    a = 'unannotated'
            else:
                a = 'unknown'


            m8 = '0'
            if mir[7] == '(':
                m8 = '1'  # If bases at mir position '8' are complementary.

            mir2_5 = str(mir[1:5].count('('))
            mir6_7 = str(mir[5:7].count('('))
            mir9_12 = str(mir[8:12].count('('))
            mir13_17 = str(mir[12:17].count('('))



            startpos='n/a'
            endpos='n/a'

            if transcript in utr_dict and 'utr3range' in utr_dict[transcript] and 'utr5range' in utr_dict[transcript]:
                
                if utr_dict[transcript]['utr3range'][1] >= mrna_pos_start >= utr_dict[transcript]['utr3range'][0]:
                    startpos = "3'UTR"
                elif utr_dict[transcript]['utr5range'][0] <= mrna_pos_start <= utr_dict[transcript]['utr5range'][1]:
                    startpos = "5'UTR"
                elif utr_dict[transcript]['utr5range'][1] <= mrna_pos_start <= utr_dict[transcript]['utr3range'][0]:
                    startpos = "CDS"
                if utr_dict[transcript]['utr3range'][1] >= mrna_pos_end >= utr_dict[transcript]['utr3range'][0]:
                    endpos = "3'UTR"
                elif utr_dict[transcript]['utr5range'][0] <= mrna_pos_end <= utr_dict[transcript]['utr5range'][1]:
                    endpos = "5'UTR"
                elif utr_dict[transcript]['utr5range'][1] <= mrna_pos_end <= utr_dict[transcript]['utr3range'][0]:
                    endpos = "CDS"


            line = '\t'.join(line)
            to_append = '\t'.join([mir_gene, mrna_gene, transcript, mir_seq, mrna_seq, mir, mrna, fe, a, mir2_5, mir6_7, m8, mir9_12, mir13_17, startpos, endpos, count])
            outfile.write("%s\t%s\n" % (line, to_append))
            line_count+=1
            if line_count % 100 == 0:
                print(line_count, "lines processed.", total_lines, "lines total.")


