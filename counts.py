import pandas as pd
import glob
import numpy as np
from collections import defaultdict
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


script_name = sys.argv.pop(0)
paths = sys.argv  # hyb output file paths for replicate samples
samples = len(paths)

mirna_path = '/Users/mac9/Downloads/snu_ctl_mir_counts.tsv'  # microRNA counts from small fraction seq (not clash)
mrna_path = '/Users/mac9/Downloads/SNU719_gene_expression.tsv'  # mRNA TPMs from RNA-seq experiment
ago_paths = glob.glob("/Users/mac9/Desktop/clash_jul2018/mircounts/*hybrids.hyb.*")  # Files output from 'count_total_clash_mirs.py'


hyb_count_min = 5

dd = defaultdict(list)

def append_counts(path):
    '''Process a CLASH hyb file, appending mir-mrna binding counts to a defaultdict.'''

    with open(path, 'r') as infile:
        header = next(infile).strip('\n').split('\t')
        mir_col, mrna_col, counts_col = [header.index(i) for i in ['mir', 'mrna', 'count']]

        for line in infile:
            line = line.split('\t')
            mir, mrna, counts = line[mir_col], line[mrna_col], line[counts_col]
            bound = "%s_%s" % (mir, mrna)
            dd[bound].append(int(counts))


def average_counts(paths):
    '''Run append_counts for each sample, then return average counts.'''
    
    for path in paths:
        append_counts(path)

    avg = {i: np.sum(dd[i]) / len(paths) for i in dd}
    snucounts = pd.DataFrame.from_dict(avg, orient='index')
    snucounts.columns = ['counts']
    snucounts['mir'] = [i.rsplit('_',1)[0] for i in snucounts.index]
    snucounts['mrna'] = [i.rsplit('_',1)[1] for i in snucounts.index]
    snucounts = snucounts.set_index('mir')

    return snucounts


def mrna_expression(mrna_path):

    mrna = pd.read_table(mrna_path, index_col=0)
    mrna = np.mean(mrna,1)
    return dict(zip(mrna.index, mrna))

def mirna_expression(mirna_path):

    mir = pd.read_table(mirna_path, index_col=0)
    mir = np.mean(mir, 1)
    return dict(zip(mir.index, mir))

def total_ago_mirna(paths):

    new = defaultdict(list)
    samples = len(paths)
    assert samples > 0, "Need at least one path"
    if isinstance(paths, str):
        paths = [paths]
    for path in paths:
        with open(path) as infile:
            header = next(infile)
            for line in infile:
                mir, counts = line.strip('\n').split('\t')
                new[mir].append(int(counts))
    for key in new:
        new[key] = np.sum(new[key]) / samples
    
    return new
    

def plot_hyb_proportions(df, filter=hyb_count_min):
    '''df argument should be the snucounts df'''

    df = df[df['counts'] > filter]

    df = df[df['mrna'] != 'MTRNR2L12']  # weird potential artifact gene
    order = list(df.reset_index().groupby('mir')['counts'].sum().sort_values(ascending=False).index)
    border = [i for i in order if 'ebv' in i]
    horder = [i for i in order if 'hsa' in i]
    order = border + horder

    df = df.loc[order]
    plt.rc("axes", linewidth=2)
    fig, ax = plt.subplots(figsize=(36,12))
    enum_df = dict(zip(df.index.unique(),range(len(df.index.unique()))))
    for i,j in enumerate(df.index.unique()):
        print(j)
        z = df.loc[j]
        ind = [i] * z['counts'].size
        counts = z['counts']
        size = counts /10
        plt.scatter(ind, counts *100/np.sum(counts),s=size, linewidth=1, edgecolor='k',alpha=.5)

    plt.title("Unique targets of microRNAs",fontsize=35, fontweight='bold')
    plt.yticks(fontsize=20,fontweight='bold')
    plt.ylabel("Percent hybridization to each transcript",fontsize=22,fontweight='bold')
    ax.set_xticks(range(len(df.index.unique())))
    ax.set_xticklabels([i.split('miR-')[1] if 'miR' in i else i.split('_')[-1] for i in df.index.unique()], rotation=90, fontsize=8)
    plt.yticks(fontsize=22, fontweight='bold')
    plt.subplots_adjust(bottom=0.2)
    ax.grid(alpha=0.3)
    ax.set_axisbelow(True)
    #plt.yscale('log')
    plt.savefig("hyb_proportions.pdf")


def targets_per_mir_seaborn(df, mircount_dict, min_interactions_per=1):
    
    sns.set(rc={'figure.figsize':(15, 15) })
    sns.set_style("white")

    ebv = df.loc[[i for i in set(df.index) if 'ebv' in i and i in mircount_dict]]
    human = df.loc[[i for i in set(df.index) if 'hsa' in i and i in mircount_dict]]



    ebv = ebv[ebv['counts'] > min_interactions_per]
    human = human[human['counts'] > min_interactions_per]

    human_unique = human.groupby(level=0)['mrna'].size().sort_values() 

    ebv_unique = ebv.groupby(level=0)['mrna'].size().sort_values()
    human_counts = [mircount_dict[i] + 1  for i in human_unique.index] 
    ebv_counts = [mircount_dict[i] + 1 for i in ebv_unique.index]
    human_unique = list(human_unique) 
    ebv_unique = list(ebv_unique) 

    g = sns.JointGrid(np.log2(np.array(ebv_counts)+1), np.log2(np.array(ebv_unique)+1))
    sns.kdeplot(np.log2(np.array(human_counts)+1), ax=g.ax_marg_x, shade=True, color="gray", label="Human")
    sns.kdeplot(np.log2(np.array(human_unique)+1),vertical=True, ax=g.ax_marg_y, shade=True, color="gray")
    g.ax_joint.plot(np.log2(np.array(human_counts)+1), np.log2(np.array(human_unique)+1), "o", ms=4, color="0.80")

    sns.kdeplot(np.log2(np.array(ebv_counts)+1), ax=g.ax_marg_x, shade=True, color="#E55300",label='EBV')
    sns.kdeplot(np.log2(np.array(ebv_unique)+1), ax=g.ax_marg_y, vertical=True, shade=True, color="#E55300")
    g.ax_joint.plot(np.log2(np.array(ebv_counts)+1), np.log2(np.array(ebv_unique)+1), "o", ms=4, color="#E55300", label='EBV')

    g.ax_joint.set_xlabel(r'$\log_2(miRNA)$', fontweight='bold', fontsize=16)
    g.ax_joint.set_ylabel(r'$\log_2(Targets)$', fontweight='bold', fontsize=16)
    g.fig.set_tight_layout(tight=True)
    #g.ax_joint.grid()
    plt.legend()
    plt.savefig('targets_per_mir.pdf')


#def species_hyb_per_region(paths):

    

snucounts = average_counts(paths)
mrna = mrna_expression(mrna_path)
mirna = mirna_expression(mirna_path)
ago = total_ago_mirna(ago_paths)


def ago_bind_by_species(mirna_dict, ago_mirna_dict, lower_lim=50):
    '''Plots ago-complexed mir counts/ total mir counts for each species'''

    mirna_dict = {i: mirna_dict[i] + 1 for i in mirna_dict}  # To avoid division of 0
    ago_mirna_dict = {i: ago_mirna_dict[i] + 1 for i in ago_mirna_dict}

    ago_d = {i: ago_mirna_dict[i] / mirna_dict[i] for i in ago_mirna_dict if mirna_dict[i] > lower_lim}

    h_arr = [ago_d[i] for i in ago_d if 'hsa' in i]
    e_arr = [ago_d[i] for i in ago_d if 'ebv' in i]

    h_mean = np.mean(h_arr)
    e_mean = np.mean(e_arr)
    p = stats.ttest_ind(h_arr, e_arr)
    #plt.violinplot([h_arr, e_arr])
    #plt.scatter([0]*len(h_arr), h_arr, c='k')
    #plt.scatter([1]*len(e_arr), e_arr, c='r')
    plt.boxplot([h_arr, e_arr])
    #plt.boxplot(e_arr)

    print(lower_lim, len(ago_d), p)