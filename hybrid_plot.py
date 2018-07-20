import matplotlib.pyplot as plt
import pandas as pd
import glob
import numpy as np
from collections import defaultdict
import os
from matplotlib.gridspec import GridSpec
from scipy.interpolate import spline
import seaborn as sns


plt.close('all')
basedir = "/home/nate/Documents/groups/Erik"

snu_files = glob.glob("%s/annotated_clash/SNU*" % basedir)

dd = defaultdict(list)

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
snucounts = pd.DataFrame.from_dict(avg, orient='index')
snucounts.columns = ['counts']
snucounts['mir'] = [i.split('_')[0] for i in snucounts.index]
snucounts['mrna'] = [i.split('_')[1] for i in snucounts.index]
snucounts.index = snucounts['mir']

df = snucounts
bart_map = pd.read_csv("/home/nate/Documents/groups/Erik/star_barts_to_3p5pmap", index_col=0)
df['barts'] = bart_map
df = df.set_index(df['barts'])
df = df.drop(['barts','mir'],1)
df = df.loc[sorted(df.index.unique(),key =lambda x:int(x.split('BART')[1].split('-')[0]))]

mrna = pd.read_table("/home/nate/Documents/groups/Erik/SNU719_gene_expression.tsv", index_col=0)
mrna = np.mean(mrna,1)
mrna = dict(zip(mrna.index, mrna))

mir = pd.read_table("%s/SNUmirsctl.tsv" % basedir, index_col=0)
mir = np.mean(mir, 1)
mir.index = mir.index.map(lambda x:x.split(',')[1])
mir = dict(zip(mir.index, mir))

plt.rc("axes", linewidth=2)
fig, ax = plt.subplots(figsize=(24,12))
enum_df = dict(zip(df.index.unique(),range(len(df.index.unique()))))
for i,j in enumerate(df.index.unique()):
    print(j)
    z = df.loc[j]
    ind = [i] * z['counts'].size
    counts = z['counts']
    size = counts 
    plt.scatter(ind, counts*100/np.sum(counts),s=size, linewidth=1, edgecolor='k',alpha=.5)

plt.title("Cellular targets of microRNAs",fontsize=35, fontweight ='bold')
plt.yticks(fontsize=20,fontweight='bold')
plt.ylabel("Percent hybridization to each transcript",fontsize=22,fontweight='bold')
ax.set_xticks(range(len(df.index.unique())))
ax.set_xticklabels([i.split('miR-')[1] if 'miR' in i else i.split('_')[-1] for i in df.index.unique()],rotation=90,fontsize=22, fontweight='bold')


plt.title("Cellular targets of EBV microRNAs",fontsize=35, fontweight ='bold')
plt.yticks(fontsize=20,fontweight='bold')
plt.ylabel("Percent hybridization to each transcript",fontsize=22,fontweight='bold')
ax.set_xticks(range(len(df.index.unique())))
ax.set_xticklabels([i.split('miR-')[1] for i in df.index.unique()],rotation=90,fontsize=22, fontweight='bold')
plt.yticks(fontsize = 22, fontweight='bold')
plt.subplots_adjust(bottom=0.2)
ax.set_ylim([0,50])
ax.set_xlim(-1,39)
ax.grid(alpha=0.7)
ax.set_axisbelow(True)
def annotate(mir, rank, y, x_adjust):
    mir = 'ebv-miR-BART'+mir
    x = enum_df[mir]
    text = df.loc[mir].sort_values('counts')['mrna'].iloc[-(rank)]
    plt.annotate(text, (x + x_adjust, y),fontsize=22,fontweight = 'bold')

annotate('3-3p',1,40,1)
annotate('5-5p', 1, 24,1)
annotate('5-5p', 2, 15.75,1)
annotate('6-3p', 1, 34.5,.5)
annotate('7-3p', 1, 10,-4)
annotate('12',1,38,.4)
annotate('16',1,23,.7)
annotate('17-5p',1,19,.7)
annotate('17-3p',1,32,.4)
annotate('18-3p',1,40,.4)
plt.savefig("/home/nate/clash_presentation/cell_bart_targets.png")




# Seed plot
dd_seed = defaultdict(list)
def append_seed_counts(path):
    with open(path, 'r') as infile:
        header = next(infile).split('\t')
        header_map = {i:j for i,j in zip(header, range(1,len(header) + 1))}
        for line in infile:
            line = line.split('\t')
            mir = line[header_map['mir']]
            mrna = line[header_map['bound']]
            counts = line[header_map['counts']]
            seed = line[header_map['mirna-mrna_seed_match_type']]
            energy = line[header_map['Predicted_binding_energy']]
            region = line[header_map['region_aligned_to_start']]
            bound = "%s_%s_%s_%s_%s" %(seed, mir, mrna, energy, region)
            dd_seed[bound].append(int(counts))

def average_seed_counts():
    for path in snu_files:
        append_seed_counts(path)
    avg = {i:np.sum(dd_seed[i])/3 for i in dd_seed}
    return avg
avg_seed = average_seed_counts()
snucounts_seed = pd.DataFrame.from_dict(avg_seed, orient='index')
snucounts_seed.columns = ['counts']
snucounts_seed['seed'] = [i.split('_')[0] for i in snucounts_seed.index]
snucounts_seed['mir'] = [i.split('_')[1] for i in snucounts_seed.index]
snucounts_seed['mrna'] = [i.split('_')[2] for i in snucounts_seed.index]
snucounts_seed['energy'] = [i.split('_')[3] for i in snucounts_seed.index]
snucounts_seed['region'] = [i.split('_')[4] for i in snucounts_seed.index]
snucounts_seed = snucounts_seed.set_index('mir')
df = snucounts_seed
bart_map = pd.read_csv("/home/nate/Documents/groups/Erik/star_barts_to_3p5pmap", index_col=0)
df['barts'] = bart_map
df = df.set_index(df['barts'])
df = df.drop('barts',1)
df = df.loc[sorted(df.index.unique(),key =lambda x:int(x.split('BART')[1].split('-')[0]))]

fig2, ax2 = plt.subplots(figsize = (24, 12))
seeds = ['No seed match', "No seed match with 3' supplementary site", 'Offset 6mer', "Offset 6mer with 3' supplementary site",
                 "6mer site", "6mer site with 3' supplementary site", "7mer-A1 site", "7mer-A1 site with 3' supplementary site", "7mer-m8 site", 
                  "7mer-m8 site with 3' supplementary site", "Mismatch 8mer", "8mer site", "8mer site with 3' supplementary site"]

seed_short = ['No seed', "No seed (3'supp)", "Offset 6mer", "Offset 6mer (3'supp)", "6mer", "6mer (3'supp)", "7mer-A1",
                "7mer-A1 (3'supp)", "7mer-m8", "7mer-m8 (3'supp)", "Mismatch 8mer", "8mer", "8mer (3'supp)"]
infod = defaultdict(list)
reg_dd = defaultdict(dict)
for i, j in enumerate(seeds):
    print(j)
    z = df[df['seed'] == j]
    en = z['energy'].map(lambda x:-float(x))
    z['energy'] = en
    energy = pd.DataFrame(z.reset_index().groupby(['barts','mrna'])['energy'].mean()).reset_index().set_index('barts')

    z =  pd.DataFrame(z.reset_index().groupby(['barts','mrna'])['counts'].sum().reset_index()).set_index('barts')
    
    z['energy'] = energy['energy']
    colors = []
    colors_dict = {'CDS': '#AB4E68', "3'UTR": '#16DB93', "5'UTR":'#F29E4C'}
    gene_exp, mir_exp, counts, e, info, reg2 = [], [], [], [], [], []
    for mr, mi, count, en in zip(z['mrna'],z.index, z['counts'], z['energy']):
        if (mi in mir) and (mr in mrna) and mrna[mr] > 10 and count > 10:
            info.append("%s_%s"%(mi, mr))
            gene_exp.append(mrna[mr])
            mir_exp.append(mir[mi])
            counts.append(count*50000) 
            e.append(en)


            
   
    ind = [i] * len(gene_exp)
    counts = counts/(np.array(gene_exp) * np.array(mir_exp))
    
    info2=dict(zip(counts, info))
    infod[j] = info2
    plt.scatter(ind, e, s=counts, linewidth=1, edgecolor='k',alpha=.4)
ax2.set_xticks(range(len(seeds)))
ax2.set_xticklabels(seed_short, rotation=90, fontsize=22, fontweight='bold')
#plt.yscale('log')
plt.subplots_adjust(bottom=0.4)
plt.ylabel(r'Binding Energy ($\mathbf{-(\Delta G)}$)', fontsize=22, fontweight='bold')
plt.yticks(fontsize=22,fontweight ='bold')
#plt.ylabel(r'$\frac{hy}{\sqrt{miR x mRNA}}$',fontsize=40)
plt.grid(True,which="both",ls="--", alpha=0.7)
ax2.set_axisbelow(True)
plt.title("Affinity of EBV microRNAs to targeted mRNA", fontsize=35, fontweight='bold')
plt.annotate('IPO7', (11,30),fontsize=16, fontweight='bold')
plt.annotate('FBXO21', (10.8,23),fontsize=16, fontweight='bold')
plt.annotate('FNDC3A', (5.8,17),fontsize=16, fontweight='bold')
plt.annotate('BTBD1', (9.4,17),fontsize=16, fontweight='bold')
plt.annotate('TNRC6A', (9.3,20),fontsize=16, fontweight='bold')
plt.annotate('BTG2', (3.2,22),fontsize=16, fontweight='bold')
plt.annotate('PANK3', (1.2,15),fontsize=16, fontweight='bold')
plt.annotate('ZFPM1', (5.2,15),fontsize=16, fontweight='bold')
plt.annotate('SOX4', (10.35,14),fontsize=16, fontweight='bold')
plt.annotate('GLO1', (9.2,30),fontsize=16, fontweight='bold')
plt.savefig("/home/nate/clash_presentation/bart_mrna_seed_match.png")

#Done

stad = pd.read_csv('/home/nate/stad_protein_coding_volc.csv',index_col=0)
x=pd.DataFrame(df.groupby(['mrna','region'])['counts'].sum()).reset_index().set_index('mrna').sort_values('counts').reset_index()
dd = defaultdict(list)

for index,row in x.iterrows():
    dd[row['mrna']].append((row['counts'],row['region']))
d = defaultdict(float)
for mrna in dd:
    x = 0
    for tup in dd[mrna]:
        if tup[0] > x:
            x = tup[0]
            y = tup[1]
        d[mrna] = y
newd = defaultdict(float)
for match in infod.keys():
    for count, hybrid in infod[match].items():
        mrna = hybrid.split('_')[-1]
        if count > newd[mrna]:
            newd[mrna] = count
colors = []
for i in stad.index:
    if i in newd:
        if d[i] == 'CDS':
            colors.append('r')
        elif d[i] == "3'UTR":
            colors.append('g')
        elif d[i] == "5'UTR":
            colors.append('b')
        else:
            colors.append('pink')
            
    else:
        colors.append((1,1,1,1))
sizes = []
for i in stad.index:
    if i in newd:
        sizes.append(max(3,newd[i]/5))
    else:
        sizes.append(3)


#volcano?
fig3,ax3 = plt.subplots()
plt.scatter(stad['fc'],stad['p'],c=colors,s=sizes,alpha=.7)

newds = list(sorted([(j,i) for i,j in newd.items()]))



# Scatter plots
fig = plt.figure(figsize=(24,24))
st = pd.read_csv('/home/nate/Documents/groups/Erik/stad_ebv_tpm.csv', index_col=0)
mirse = pd.read_csv('/home/nate/Documents/groups/Erik/mirs_ebv_packaged.csv', index_col=0)
plot_these = ['IPO7', 'FBXO21', 'FNDC3A', 'BTBD1', 'JKAMP', 'NUP35', 'PCGF5', 'NCOA7', 'PITRM1', 'TNRC6A','BTG2', 'CD47']
corresponding_mirs =['ebv-miR-BART3-3p', 'ebv-miR-BART21-5p','ebv-miR-BART20-3p', 'ebv-miR-BART5-5p', 'ebv-miR-BART5-5p', 'ebv-miR-BART21-5p',
'ebv-miR-BART20-3p', 'ebv-miR-BART5-5p', 'ebv-miR-BART21-5p', 'ebv-miR-BART21-5p', 'ebv-miR-BART9-3p', 'ebv-miR-BART21-5p']
row=1
color_dict = {"3'UTR":'#87255B', 'CDS': '#D7263D', "5'UTR":'#FF570A'}
for ind, (gene,micro) in enumerate(zip(plot_these,corresponding_mirs)):
    sub = ind+1
    ax = plt.subplot(3,4,sub)
    color = color_dict[d[gene]] if d[gene] in color_dict else 'pink'
    
    plt.scatter(mirse.loc[micro], st.loc[gene],c=color, lw=1,s=30, edgecolor='k')
    plt.yticks(fontsize=18, fontweight='bold')
    plt.xticks(fontsize=18, fontweight='bold')
    plt.ylabel(gene, fontsize=18,fontweight='bold')
    plt.xlabel(micro, fontsize=18, fontweight='bold')
    plt.title(gene, fontsize=30, fontweight='bold')
    fit = np.polyfit(mirse.loc[micro], st.loc[gene], 1)
    plt.plot(mirse.loc[micro], fit[0]*mirse.loc[micro]+fit[1], c='k')

plt.tight_layout()



# GSEA

gsea_basedir = '/home/nate/Documents/groups/Erik/clash_gsea'
gsea_prefix = 'stad_gsea'

def gen_grp(cutoff=20):
    
    for region in ["5'UTR", "3'UTR", "CDS"]:

        nu = [i for  i in newds if d[i[1]] == region]
        region_text = region.replace("'","")
        with open("{base}/{region}_targets.grp".format(base=gsea_basedir, region=region_text), "w") as outfile:
            outfile.write(">%s\n" % region)
            for i in nu:
                if i[0] >= cutoff:
                    outfile.write("%s\n" %(i[1]))
        
        os.system("java -cp ~/Programs/gsea-3.0.jar -Xmx10000m "
                    "xtools.gsea.Gsea -res {base}/stad_gsea.txt -cls {base}/stad_gsea.cls "
                    "-gmx {region}_targets.grp -collapse false -out {base}/{region} -set_min 3".format(base=gsea_basedir, region=region_text))




input_colors = ["#0EAD69","#29335C","#DD6031"]
number_of_genes = 20502 # From tcga firebrowse dataset

def import_df(region):
    
    region_text = region.replace("'","")

    rand_folder_name = os.popen("ls -trlh {base}/{region}".format(base=gsea_basedir, region=region_text)).read().split('\n')[-2].split()[-1]
    path = "{}/{}/{}/{}_targets.grp.xls".format(gsea_basedir, region_text, rand_folder_name, region_text)
    df = pd.read_table(path, index_col=0, usecols=["RANK IN GENE LIST", "RUNNING ES"])
    df.loc[number_of_genes] = 0 # Ensure the lines end at 0
    df.iloc[0] = 0
    rank = np.array(df.index)
    es = np.array(df['RUNNING ES'])

    return rank, es

def plot(regions=("CDS","3'UTR", "5'UTR"), spline_points=1000):

    colors = []
    ax = plt.subplot(gs[:10,:]) # GridSpec

    num_regions = len(regions)
    yarr = np.zeros([num_regions, spline_points]) # Initialize an array
    xnew = np.linspace(0, number_of_genes, spline_points)

    for ind, region in enumerate(regions):
        print(region)
        rank, es = import_df(region=region)
        ynew = spline(rank, es, xnew, order=1, kind='smoothest')
        active_color = input_colors[ind]
        p = plt.plot(xnew, ynew, linewidth=3, label=region, c=active_color)
        colors.append(p[0].get_color())

        yarr[ind] = ynew  # Store the spline in an array for later

    for i in range(len(yarr)):
        sorted_curves = np.argsort(np.min(yarr, 1)) # Sorted indices of minimum values from low to high
        min_curve = sorted_curves[i]
        y0 = yarr[min_curve]
        c = colors[min_curve]
        plt.fill_between(xnew, y0, np.zeros(1000), color=c, alpha=.06) # Fill AUC

    for ind, region in enumerate(regions):
        active_color = input_colors[ind]

        rank, es = import_df(region)
        ax2 = plt.subplot(gs[12+ind,:], sharex=ax)
        ax2.axis('off')
        plt.plot(rank, [0]*len(rank),'|',ms=30, c=active_color)

 
    plt.xlim([-1, number_of_genes + 1])
    ax.set_xlabel('Rank', fontsize=23, fontweight='bold')
    ax.set_ylabel('ES', fontsize=23, fontweight='bold')
   # ax.set_title(pathway.replace('_',' '), fontsize=26, fontweight='bold')
    ax.set_title("EBV microRNA targets", fontsize=26, fontweight='bold')

    plt.sca(ax)
    plt.xticks(fontsize=20, fontweight='bold')
    plt.yticks(fontsize=20, fontweight='bold')
    ax.axhline(0, c='k', linewidth=4)
    ax.grid()
    ax.legend(loc='upper right',fancybox=True, fontsize=16)
    #ax.axvline(8000, ls='dashed', lw=2, c='k')

gen_grp(cutoff=15)
gs = GridSpec(16,1)     

plt.rc('axes', linewidth=2)
fig = plt.figure(figsize=(12,8)) # Initiate figure here in case I decide to loop through pathways and create subplots
plot()

plt.savefig("clash_gsea_stacked.png" ,dpi=300)  
plt.close('all')

#plt.subplots_adjust(hspace=0.3, wspace=0.3, left=.1, right=.9)
#plt.tight_layout()



#corr heatmaps

def heatplot(region, count_threshold,count_max_threshold=10000):
    
    nu = [i[1] for  i in newds if d[i[1]]==region and i[0]>=count_threshold and i[0]<count_max_threshold]

    for i in mirse.index:                                                          
        if 'BH' not in i and 'mir' not in i:
            nu.append(i)
    
    st2 = st.T
    mirse2 = mirse.copy() 
    mirse2 = mirse2.T
    st2[mirse2.columns] = mirse2 

    cc=st2.T.loc[nu].dropna().T.corr()

    c = sns.clustermap(cc,robust=True,cmap='RdBu_r',figsize=(36,36),linewidth=.1, linecolor='k',vmin=-1,vmax=1)
    plt.setp(c.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.setp(c.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.savefig("{region}.thresh_{thresh}_{max}.correlation_plot.png".format(region=region,thresh=count_threshold,max=count_max_threshold), dpi=300)

regions=("CDS","3'UTR", "5'UTR")
for i in regions:
    try:
        heatplot(region=i, count_threshold=50,count_max_threshold=1000000)
    except:
        print("too few meet threshold", i)

def heatplot_all(count_threshold, count_max_threshold=10000):
    
    nu = [i[1] for  i in newds if i[0]>=count_threshold and i[0]<count_max_threshold]

    for i in mirse.index:                                                          
        if 'BH' not in i and 'mir' not in i:
            nu.append(i)
    
    st2 = st.T
    mirse2 = mirse.copy() 
    mirse2 = mirse2.T
    st2[mirse2.columns] = mirse2 

    cc=st2.T.loc[nu].dropna().T.corr()

    c = sns.clustermap(cc,robust=True,cmap='RdBu_r',figsize=(36,36), method='ward', vmin = -1, vmax = 1, linewidth=1,linecolor='k')
    plt.setp(c.ax_heatmap.xaxis.get_majorticklabels(), rotation=90, fontsize=16)
    plt.setp(c.ax_heatmap.yaxis.get_majorticklabels(), rotation=0,fontsize=16)
    plt.savefig("/home/nate/clash_presentation/all.thresh_{thresh}_{max}.correlation_plot.png".format(thresh=count_threshold,max=count_max_threshold), dpi=300)
heatplot_all(150,100000)
heatplot_all(0,0.5)


stadmirs = pd.read_table("/home/nate/Documents/groups/Erik/human_ebv_CPM_STAD.tsv",index_col=0)
k = [i for i in stadmirs.index if 'ebv-miR-BART' in i]


fig = plt.figure(figsize=(24,12))
ax = plt.subplot(111)
plt.bar(range(31),np.sum(stadmirs.ix[k,-31:]),ec='k',color='#DD6031', label='EBV') #31 ebv+ patients all at end of df
plt.bar(range(31),[1000000]*31-np.sum(stadmirs.ix[k,-31:]),bottom=np.sum(stadmirs.ix[k,-31:]),ec='k',color='#29335C', label='Human')
ax.set_yticklabels([0,20,40,60,80,100,120], fontsize=22, fontweight='bold')
ax.set_xticklabels('')
ax.set_xlabel("Stomach adenocarcinoma patients", fontsize=22,fontweight='bold')
ax.set_ylabel("Percent of total microRNAs", fontsize=22, fontweight='bold')
plt.title('EBV and human microRNAs in stomach adenocarcinoma patients', fontsize=35, fontweight='bold')
ax.grid(alpha=0.7)
ax.set_axisbelow(True)
ax.legend(loc='best', fontsize=22)
plt.savefig("/home/nate/clash_presentation/ebv_counts_by_percent_tcga.png")


stadmirs = stadmirs.loc[k]
t = [i for i in stadmirs if i[13:15]=='01' ]
stadmirs_t = stadmirs[t] 
n = [i for i in stadmirs if i[13:15] == '11']
stadmirs_n = stadmirs[n]
fig = plt.figure(figsize=(24,12))
ax = plt.subplot(111)

plt.scatter(range(len(stadmirs_t.columns)),np.sum(stadmirs_t),s=40,linewidth=1, edgecolor='k')
ax.axvline(402, c='k',ls='--')
ax.axhline(19350, c='k',ls='--')
plt.fill_between([402,437], 19403, ax.get_ylim()[1],alpha=.2)
plt.xlim(0,437)
plt.ylim(0,401808)
plt.yticks(fontsize=22, fontweight='bold')
plt.xticks(fontsize=22, fontweight='bold')
ax.set_ylabel("EBV microRNA counts",fontsize=22,fontweight='bold')
ax.set_xlabel("Stomach adenocarcinoma patients", fontsize=22,fontweight='bold')
plt.title("EBV+ stomach adenocarcinomas", fontsize=35,fontweight='bold')
plt.annotate("7.1%", (405,380000), fontsize=22, fontweight='bold')
plt.savefig('/home/nate/clash_presentation/ebv_pos_stad.png')

fig = plt.figure(figsize=(24,12))
ax = plt.subplot(111)
plt.scatter(range(len(stadmirs_n.columns)),np.sum(stadmirs_n),s=40,linewidth=1, edgecolor='k')
plt.ylim(0,401808)
plt.yticks(fontsize=22, fontweight='bold')
plt.xticks(fontsize=22, fontweight='bold')
ax.axhline(19350, c='k',ls='--')
ax.set_ylabel("EBV microRNA counts",fontsize=22,fontweight='bold')
ax.set_xlabel("Stomach adenocarcinoma patients", fontsize=22,fontweight='bold')
plt.title("Adjacent normal tissue", fontsize=35,fontweight='bold')
plt.annotate("0%", (40,380000), fontsize=22, fontweight='bold')
plt.savefig('/home/nate/clash_presentation/ebv_neg_normals.png')