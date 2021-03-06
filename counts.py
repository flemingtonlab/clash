import pandas as pd
import glob
import numpy as np
from collections import defaultdict
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import matplotlib.patches as mpatches

#script_name = sys.argv.pop(0)
paths = glob.glob("../etc/annotated_hyb_files/*")  # hyb output file paths for replicate samples
samples = len(paths)

mirna_path = '../etc/snu_ctl_mir_counts.tsv'  # microRNA counts from small fraction seq (not clash)
mrna_path = '../etc/SNU719_gene_expression.tsv'  # mRNA TPMs from RNA-seq experiment
ago_paths = glob.glob("../etc/hyb_total_mir_counts/*hybrids.hyb.*")  # Files output from 'count_total_clash_mirs.py'

hyb_count_min = 5


def get_species(string):
    if 'hsa' in string:
        return 'Human'
    else:
        return 'EBV'


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
    
    sns.set(rc={'figure.figsize':(15, 15)})
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
    g.ax_joint.plot(np.log2(np.array(human_counts)+1), np.log2(np.array(human_unique)+1), "o", ms=9, color="0.80", mec='k', mew=.3, alpha=.7)

    sns.kdeplot(np.log2(np.array(ebv_counts)+1), ax=g.ax_marg_x, shade=True, color="#E55300",label='EBV')
    sns.kdeplot(np.log2(np.array(ebv_unique)+1), ax=g.ax_marg_y, vertical=True, shade=True, color="#E55300")
    g.ax_joint.plot(np.log2(np.array(ebv_counts)+1), np.log2(np.array(ebv_unique)+1), "o", ms=9, color="#E55300", label='EBV', mec='k', mew=.3, alpha=.7)

    g.ax_joint.set_xlabel(r'$\log_2(miRNA)$', fontweight='bold', fontsize=16)
    g.ax_joint.set_ylabel(r'$\log_2(Targets)$', fontweight='bold', fontsize=16)
    g.fig.set_tight_layout(tight=True)
    g.ax_joint.grid()
    plt.legend()
    g.ax_joint.set_xlim([-1.3, 21])
    g.ax_joint.set_ylim([-1.3, 13.5])
    plt.savefig('targets_per_mir.pdf')


snucounts = average_counts(paths)
snucounts = snucounts[snucounts['mrna'] != 'MTRNR2L12']
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

    col_list=['grey', 'orange']
    col_list_palette = sns.xkcd_palette(col_list)
    sns.set_palette(col_list)
    fig, ax = plt.subplots(figsize=(8,12))
    sns.boxplot(data=[h_arr, e_arr], ax=ax, linewidth=2, width=.5)
    ax.set_ylabel(r'$\frac{[AGO-miRNA]}{miRNA}$' , fontsize=25, fontweight='bold')
    ax.tick_params(labelsize=20, which='both', length=4)
    ax.set_xticklabels(["Human", "EBV"])
    plt.yscale('log') 
    plt.title("miRNA-Argonaute Binding Efficiency", fontsize=25)
    print(lower_lim, len(ago_d), p)
    plt.tight_layout()
    plt.savefig("ago-bound_vs_total_byspecies.pdf")
    #return ax

def grouped(df, groupby2='m13_17'):

    h = df.groupby(['species', groupby2])['count'].sum().loc['Human']
    e = df.groupby(['species', groupby2])['count'].sum().loc['EBV']

    return h, e


def get_total_13_17(groupby='m13_17'):

    hd = defaultdict(list)
    ed = defaultdict(list)
    for i in paths:
        hyb = pd.read_table(i)
        hyb = hyb[hyb['mrna'] != '']
        counts = np.sum(hyb['count'])
        hyb['species'] = hyb['mir'].apply(get_species)
        h, e = grouped(hyb, groupby2=groupby)
        for key, val in h.items():
            hd[key].append(val/counts)
        for key, val in e.items():
            ed[key].append(val/counts)
    hmn = []
    for i in hd.keys():
        
        hmn.append(mean(hd[i]))
    hmn = sum(hmn)
    
    emn = []
    for i in ed.keys():

        emn.append(mean(ed[i]))
    emn = sum(emn)

    #hmn = np.mean(hd[relative_to]) 
    for key, values in hd.items(): 
        hd[key] = [100*value / hmn for value in values]
   
    # emn = np.mean(ed[relative_to])
    for key, values in ed.items(): 
        ed[key] = [100*value / emn for value in values]

    h_sem = {i: stats.sem(hd[i]) for i in hd}
    e_sem = {i: stats.sem(ed[i]) for i in ed}

    # Call plot function here
    return hd, ed

def get_individual_13_17(groupby='m13_17'):

    hd = defaultdict(list)
    ed = defaultdict(list)
    relative_to = 0

    for i in paths:
        hyb = pd.read_table(i)
        hyb = hyb[hyb['mrna'] != '']
        counts = np.sum(hyb['count'])
        hyb['species'] = hyb['mir'].apply(get_species)
        hyb = pd.DataFrame(hyb.groupby(['species','mir','mrna',groupby])['count'].size()).reset_index().set_index('species')
        for key, val in hyb.loc['Human'][groupby].value_counts().items():
            hd[key].append(val / counts)
        for key, val in hyb.loc['EBV'][groupby].value_counts().items():
            ed[key].append(val / counts)

    hmn = []
    for i in hd.keys():
        
        hmn.append(mean(hd[i]))
    hmn = sum(hmn)
    
    emn = []
    for i in ed.keys():

        emn.append(mean(ed[i]))
    emn = sum(hmn)

    #hmn = np.mean(hd[relative_to]) 
    for key, values in hd.items(): 
        hd[key] = [100*value / hmn for value in values]
   
    # emn = np.mean(ed[relative_to])
    for key, values in ed.items(): 
        ed[key] = [100*value / emn for value in values]

    h_sem = {i: stats.sem(hd[i]) for i in hd}
    e_sem = {i: stats.sem(ed[i]) for i in ed}

    return hd, ed

def kd_prep():
    hd = defaultdict(list)
    ed = defaultdict(list)
    samples = len(paths)
    for path in paths:
        hyb = pd.read_table(path)
        hyb = hyb[~ hyb['mrna'].str.startswith('MT')]
        hyb['a1'] = pd.to_numeric(hyb['a1'], errors='coerce')
        hyb['m2_5'] = pd.to_numeric(hyb['m2_5'], errors='coerce')
        hyb['m6_7'] = pd.to_numeric(hyb['m6_7'], errors='coerce')
        hyb['m8'] = pd.to_numeric(hyb['m8'], errors='coerce')
        hyb['seed_total'] = np.sum(hyb[['a1','m2_5','m6_7','m8']], 1)
        counts = np.sum(hyb['count'])

        hyb['species'] = hyb['mir'].apply(get_species)
        sm = []
        for i in hyb.index:
            if hyb.loc[i,'seed_total'] == 8:
                if hyb.loc[i, 'm13_17'] >= 3:
                    sm.append('8mer-supp')
                else:
                    sm.append('8mer')
            elif hyb.loc[i,'seed_total'] == 7:
                if hyb.loc[i,'m13_17'] >=3:
                    sm.append('7mer-supp')
                else:
                    sm.append('7mer')
            elif hyb.loc[i, 'seed_total'] == 6:
                if hyb.loc[i, 'm13_17'] >= 3:
                    sm.append('6mer-supp')
                else:
                    sm.append('6mer')
            else:
                if hyb.loc[i, 'm13_17'] >= 3:
                    sm.append('noseed-supp')
                else:    
                    sm.append('noseed')
        
        
            print(i)
        hyb['anno'] = sm
        hyb = pd.DataFrame(hyb.groupby(['species','mir','mrna','binding_energy','mrna_feature_stop', 'anno'])['count'].sum()).reset_index() 
        for i in hyb[hyb['species']=='Human'].index:
            hd["%s_%s_%s_%s_%s"% (hyb.loc[i, 'mir'], hyb.loc[i, 'mrna'], hyb.loc[i,'mrna_feature_stop'], hyb.loc[i, 'binding_energy'], hyb.loc[i, 'anno'])].append(hyb.loc[i,'count'])

        for i in hyb[hyb['species']=='EBV'].index:
            ed["%s_%s_%s_%s_%s"% (hyb.loc[i, 'mir'], hyb.loc[i, 'mrna'], hyb.loc[i,'mrna_feature_stop'], hyb.loc[i, 'binding_energy'], hyb.loc[i, 'anno'])].append(hyb.loc[i,'count'])
   
    h_seed = defaultdict(list)
    e_seed = defaultdict(list)    
    for i in hd:
        print(i)
        mir_i, mrna_i, feature, binding, seed_i = i.split('_')
        avg = sum(hd[i])/samples
        try:
            denom = mirna[mir_i] * mrna[mrna_i]
        except KeyError:
            print('i not found in mirna/mrna dict')
            continue
        if denom < 500:
            print('mirna/mrna is 0 - cannot be denominator..skipped')
            continue
        else:
            val = avg / denom
        h_seed[seed_i].append((mir_i, mrna_i, feature, binding, val))


    for i in ed:
        mir_i, mrna_i, feature, binding, seed_i = i.split('_')
        avg = sum(ed[i])/samples
        try:
            denom = mirna[mir_i] * mrna[mrna_i]
        except KeyError:
            print('i not found in mirna/mrna dict', i)
            continue
        if denom < 500:
            print('mirna/mrna is 0 - cannot be denominator..skipped', i )
            continue
        else:
            val = avg / denom
        e_seed[seed_i].append((mir_i, mrna_i, feature, binding, val))
        
    return e_seed, h_seed

def plot_kd(e_seed, h_seed):

    seed_order = ['noseed', 'noseed-supp', '6mer', '6mer-supp', '7mer',
                    '7mer-supp', '8mer', '8mer-supp']


    fig = plt.figure(figsize=(24,12))
    ax = plt.subplot(111)
    e_colors = {'CDS': "#c92400", "3'UTR": "#132ad8", "5'UTR": '#efb907', "n/a":"gray"}
    #e_colors = {'CDS': "#c92400", "3'UTR": "#132ad8", "5'UTR": '#efb907', "n/a":"gray"}

    for x, seed in enumerate(seed_order):
        color = [e_colors[i[-3]] for i in e_seed[seed]]
        be =[-float(i[-2]) for i in e_seed[seed]]
        kd = [i[-1] * 10000 for i in e_seed[seed]]
        plt.scatter([x+.1] * len(kd), be,  c=color,s=kd, alpha=.4)

        color = [e_colors[i[-3]] for i in h_seed[seed]]
        be =[-float(i[-2]) for i in h_seed[seed]]
        kd = [i[-1] * 10000 for i in h_seed[seed]]
        plt.scatter([x-.1] * len(kd), be,  c=color,s=kd, alpha=.4)

    ax.set_xticks(range(len(seed_order)))
    ax.set_xticklabels(seed_order, rotation=90 , fontsize=22, fontweight='bold')
    plt.ylabel(r'Binding Energy ($\mathbf{-(\Delta G)}$)', fontsize=22, fontweight='bold')
    plt.yticks(fontsize=22,fontweight ='bold')
    plt.grid(True,which="both",ls="--", alpha=0.7)
    ax.set_axisbelow(True)
    plt.title('microRNA Affinity for Target', fontsize=35, fontweight='bold')
    utr3 = mpatches.Patch(color='#132ad8', label="3' UTR")
    utr5 = mpatches.Patch(color='#efb907', label="5' UTR")
    cds = mpatches.Patch(color='#c92400', label="CDS")

    plt.legend(handles=[utr3, utr5, cds], fontsize=18, frameon=True, edgecolor='k', shadow=True, loc='upper right')
    plt.subplots_adjust(bottom=0.3)
    plt.savefig('../figs/binding_affinity_colors.png')
