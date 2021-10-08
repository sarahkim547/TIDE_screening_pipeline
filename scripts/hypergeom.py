'''
Author: Mudra Hegde
Email: mhegde@broadinstitute.org
Performs hypergeometric analysis for genetic perturbation screens
Input: 1. Input file with column headers
       2. Chip File
       3. Option to include/exclude singletons in output file;Default: Singletons included
       4. Option to specify n% for mean p-val, LFC calculation
       5. Option of number of control guides to be included in a random set
       6. Option of minimum number of perturbations for a gene to be included in the volcano plot
       7. Option of maximum number of perturbations for a gene to be included in the volcano plot
       8. Option of number of genes to be labeled on the plot
edits made by Peter Du (Mar 2018)
v2.3.2 - added option to specify max best guides to average for "Average LFC" calculation; aesthetic changes
v2.3.3 - control genes start with "NO_CURRENT" or "NO_SITE" or "INACTIVE" or "ONE_INTERGENIC"
'''
import pandas as pd
import numpy as np
from scipy import stats
from math import log10
import csv, argparse, os, sys, re
from datetime import datetime
from decimal import Decimal
import networkx as nx
import matplotlib as mpl
mpl.use('Agg')
mpl.rc('pdf',fonttype=42)
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import warnings

warnings.filterwarnings("ignore")

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file',
        type=str,
        help='File containing sgRNA sequence and score(s) columns with headers',
        required=True)
    parser.add_argument("--output-dir", 
	type=str,
	required=True,
	help="Directory where result is written") 
    parser.add_argument('--chip-file',
        type=str,
        help='Chip file',
        required=True)
    parser.add_argument('--sing-pert',
        type=str,
        default='Y',
        help='Y to include genes with single perturbations in output file, N to exclude; Default: Y')
    parser.add_argument('--fraction',
        type=int,
        default=100,
        help='Fraction of perturbations to be included in mean LFC, pval calculation, 1 - 100; Default: 100')
    parser.add_argument('--ctrls-num',
        type=int, 
        default=4,
        help='Number of control guides to be included in a random set; default = 4')
    parser.add_argument('--min-pert',
        type=int,
        default=1,
        help='Minimum number of perturbations for a gene to be included in the volcano plot; default = 1')
    parser.add_argument('--max-pert',
        type=int,
        default=8,
        help='Maximum number of perturbations for a gene to be included in the volcano plot; default = 8')
    parser.add_argument('--label-num',
        type=int,
        default=3,
        help='Number of genes to be labeled on the plot; default = 3')
    parser.add_argument('--plot',
                        type=str,
                        default=False,
                        help='Generate plot; default = N')
    parser.add_argument('--max',
                        type=int,
                        default=-1,
                        help='Only consider the best N guides for all calculations; default = all')
    return parser

'''
Sorts a data frame on the column specified and re-indexes
Argument: Dataframe, column to sort on, direction of sorting
Return Value: Sorted and re-index dataframe
'''
def sort_reindex(df, col, direction):
    if direction == 'P' or direction == 'p':
        df = df.sort_values([col], ascending=False)
        df.index = range(0, len(df))
        #df['Rank'] = map(lambda x:x+1, range(len(df)))
    elif direction == 'N' or direction == 'n':
        df = df.sort_values([col], ascending=True)
        df.index = range(0, len(df))
    else:
        print 'Please enter a relevant direction; P for positive and N for negative'
        sys.exit(1)
    return df

def calc_hypergeom_scores(merged, st_in, ge, frac):
    grouped = merged.groupby('Gene Symbol')
    sps = list(st_in['Guide Sequence'])
    tot_sps = len(sps)
    sp_rank = {}
    for i in range(len(sps)):
        sp_rank[sps[i]] = i+1
    g_p_val = {}
    for g in ge:
        ge_df = grouped.get_group(g)
        ge_df = ge_df.drop_duplicates('Guide Sequence')
        gene = g
        guides = list(ge_df['Guide Sequence'])
        guide_list = ';'.join(guides)
        length_sps = len(guides)
        rank = list()
        rank_hash = {}
        for s in guides:
            r=sp_rank[s]
            rank.append(r)
            rank_hash[s] = r
        rank.sort()
        lfcs = list(ge_df['Score'])
        lfcs.sort()
        op_lfcs = ';'.join(str(round(l,2)) for l in lfcs)
        #op_lfcs = ';'.join([str(l) for l in lfcs])
        all_ranks = ';'.join(str(r) for r in rank)

        # math starts here
        if 0 < max_guide < length_sps:
            while len(lfcs) > max_guide:
                lfcs = lfcs[:-1]
        if 0 < max_guide < length_sps:
            while len(rank) > max_guide:
                rank = rank[:-1]
        p_values = [-log10(stats.hypergeom.pmf(rank.index(x)+1, tot_sps, length_sps, x)) for x in rank]
        p_values.sort(reverse=True)
        num_perts = int(frac/100.0*length_sps)
        avg_p_val = np.mean(p_values[0:num_perts])
        avg_lfc = np.mean(lfcs[0:num_perts])
        op_p_vals = ';'.join("%.2E" %Decimal(p) for p in p_values)
        g_p_val[g] = op_lfcs+'_'+str(avg_lfc)+'_'+op_p_vals+'_'+str(avg_p_val)+'_'+guide_list+'_'+all_ranks
    return g_p_val


def generate_chip(chip,num):
    print 'Generating temp chip file...'
    ctrls = chip[chip['Gene Symbol'].str.startswith('NO_SITE') |
                 chip['Gene Symbol'].str.startswith('NO_CURRENT') |
                 chip['Gene Symbol'].str.startswith('INACTIVE') |
                 chip['Gene Symbol'].str.startswith('ONE_INTERGENIC')]
    new_chip = chip[~chip['Gene Symbol'].str.startswith('NO_SITE') |
                    ~chip['Gene Symbol'].str.startswith('NO_CURRENT') |
                    ~chip['Gene Symbol'].str.startswith('INACTIVE') |
                    ~chip['Gene Symbol'].str.startswith('ONE_INTERGENIC')]
    new_chip = new_chip[['Guide Sequence','Gene Symbol']]
    random_assign = []
    for x in range(1,(len(ctrls)/num)+1):
        random_assign.extend([x]*num)
    ctrls.index = range(0,len(ctrls))
    ctrls = ctrls.ix[0:len(random_assign)-1,]
    ctrls['random_assign'] = random_assign
    new_ctrl = pd.DataFrame()
    for i in range(1,(len(ctrls)/num)+1):
        ra_df = ctrls[ctrls['random_assign'] == i]
        ra_df['Gene Symbol'] = 'NO_SITE_'+str(i)
        if len(ra_df) == num:
            new_ctrl = new_ctrl.append(ra_df)
    if len(new_ctrl) > 0:
        new_ctrl = new_ctrl[['Guide Sequence','Gene Symbol']]
        new_chip = new_chip.append(new_ctrl)
    return new_chip

def plot_volcano(outputfile,min_pert,max_pert,label_num,c):
    print 'Generating volcano plot...'
    output_df = pd.read_table(outputfile)
    output_df = output_df[(output_df['Number of perturbations']>=min_pert)&(output_df['Number of perturbations']<=max_pert)]
    ctrls = output_df[output_df['Gene Symbol'].str.contains('NO_SITE') |
                      output_df['Gene Symbol'].str.startswith('NO_CURRENT') |
                      output_df['Gene Symbol'].str.startswith('INACTIVE') |
                      output_df['Gene Symbol'].str.startswith('ONE_INTERGENIC')]
    fig,ax = plt.subplots()
    ax.scatter(output_df['Average LFC'],output_df['Average -log(p-values)'], color='0.5', alpha=0.8, edgecolor='', zorder=1, label='Genes')
    ax.scatter(ctrls['Average LFC'],ctrls['Average -log(p-values)'], color='0.2', alpha=0.8, edgecolor='', zorder=2, label='Control guides')
    ax.set_title(c,fontname="Helvetica",fontsize=12)
    ax.set_xlabel('Average fold change(log2)',fontname="Helvetica",fontsize=12)
    ax.set_ylabel('Average p-value(-log10)',fontname="Helvetica",fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    for tick in ax.get_xticklabels():
        tick.set_fontname('Helvetica')
        tick.set_fontsize(12)
    for tick in ax.get_yticklabels():
        tick.set_fontname('Helvetica')
        tick.set_fontsize(12)    
    output_df = output_df.sort_values(['Average LFC'],ascending=False)

    top_hits = output_df.head(label_num)
    for i,r in top_hits.iterrows():
        ax.scatter(r['Average LFC'],r['Average -log(p-values)'], c='red', alpha=0.9, edgecolor='', zorder=3)
    repel_labels(ax, top_hits['Average LFC'], top_hits['Average -log(p-values)'], top_hits['Gene Symbol'], k=0.2)

    bottom_hits = output_df.tail(label_num)
    for i,r in bottom_hits.iterrows():
        ax.scatter(r['Average LFC'], r['Average -log(p-values)'], c='red', alpha=0.9, edgecolor='', zorder=3)
        # ax.annotate(r['Gene Symbol'],
        #             xy=(r['Average LFC'], r['Average -log(p-values)']), xycoords='data',
        #             xytext=(5, -1), textcoords='offset points', fontSize=10
        #             )
    repel_labels(ax, bottom_hits['Average LFC'], bottom_hits['Average -log(p-values)'], bottom_hits['Gene Symbol'], k=0.2)
    #plt.xticks([-15, -10, -5, 0, 5, 10, 15])
    plt.legend(loc='lower right')
    fig.savefig(o_folder+'/'+c+'.pdf')
    return 1

# k is repulsion factor
def repel_labels(ax, x, y, labels, k=0.01):
    G = nx.DiGraph()
    data_nodes = []
    init_pos = {}
    for xi, yi, label in zip(x, y, labels):
        data_str = 'data_{0}'.format(label)
        G.add_node(data_str)
        G.add_node(label)
        G.add_edge(label, data_str)
        data_nodes.append(data_str)
        init_pos[data_str] = (xi, yi)
        init_pos[label] = (xi, yi)

    pos = nx.spring_layout(G, pos=init_pos, fixed=data_nodes, k=k)

    # undo spring_layout's rescaling
    pos_after = np.vstack([pos[d] for d in data_nodes])
    pos_before = np.vstack([init_pos[d] for d in data_nodes])
    scale, shift_x = np.polyfit(pos_after[:,0], pos_before[:,0], 1)
    scale, shift_y = np.polyfit(pos_after[:,1], pos_before[:,1], 1)
    shift = np.array([shift_x, shift_y])
    for key, val in pos.items():
        pos[key] = (val*scale) + shift

    for label, data_str in G.edges():
        ax.annotate(label,
                    xy=pos[data_str], xycoords='data',
                    xytext=pos[label], textcoords='data',)
                    # arrowprops=dict(arrowstyle="-",
                    #                 shrinkA=0, shrinkB=0,
                    #                 connectionstyle="arc3",
                    #                 color='k'))
    # expand limits
    # all_pos = np.vstack(pos.values())
    # x_span, y_span = np.ptp(all_pos, axis=0)
    # mins = np.min(all_pos-x_span*0.15, 0)
    # maxs = np.max(all_pos+y_span*0.15, 0)
    # ax.set_xlim([mins[0], maxs[0]])
    # ax.set_ylim([mins[1], maxs[1]])


if __name__ == '__main__':
    args = get_parser().parse_args()
    plot = args.plot
    max_guide = args.max
    inputfile = args.input_file
    input_df = pd.read_table(inputfile)
    cols = list(input_df.columns)[1:]
    cols = [re.sub('[\s+\.]', '_', x) for x in cols]
    cols.insert(0,'Guide Sequence')
    input_df.columns = cols
    cols_iter = cols[1:]
    ref = pd.read_table(args.chip_file)
    ref_colnames = list(ref.columns)
    ref_colnames[0:2] = ['Guide Sequence', 'Gene Symbol']
    ref.columns = ref_colnames
    include = args.sing_pert
    frac = args.fraction
    num = args.ctrls_num
    min_pert = args.min_pert
    max_pert = args.max_pert
    label_num = args.label_num
    inputname = inputfile.split('/')[-1]
    ref = generate_chip(ref,num)
    o_folder = args.output_dir #os.path.dirname(inputfile)+'/volcano_'+inputname[:-4]+'_'+str(datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
    if not os.path.exists(o_folder):
        os.makedirs(o_folder)
    print '\nOutput will be written to ' + o_folder
    ref.to_csv(o_folder+'/temp_chip_file.txt',sep='\t',index=False)
    counter = 1
    for ci,c in enumerate(cols_iter):
        print 'Analyzing ... ' + str(counter) + '/' + str(len(cols_iter))
        outputfile = o_folder+'/'+c+'_volcano_frac'+str(frac)
        if max_guide > 0:
            outputfile = outputfile + '_max' + str(max_guide)
        outputfile = outputfile + '.txt'
        st_in = input_df[['Guide Sequence',c]]
        st_in = st_in.rename(columns={c:'Score'})
        merged = pd.merge(st_in, ref, on='Guide Sequence')
        ge=list(merged['Gene Symbol'])
        ge=list(set(ge))

        merged = sort_reindex(merged, 'Score', 'P')
        st_in = sort_reindex(st_in, 'Score', 'P')
        g_p_val_P = calc_hypergeom_scores(merged, st_in, ge, frac)
        merged = sort_reindex(merged, 'Score', 'N')
        st_in = sort_reindex(st_in, 'Score', 'N')
        g_p_val_N = calc_hypergeom_scores(merged, st_in, ge, frac)

        with open(outputfile,'w') as o:
            w = csv.writer(o, delimiter='\t', lineterminator='\n')
            w.writerow(('Gene Symbol', 'Average LFC', 'Average -log(p-values)', 'Number of perturbations', 'Perturbations', 'Individual LFCs', 'Ascending ranks','Individual ascending -log(p-values)', 'Descending ranks', 'Individual descending -log(p-values)'))
            for g in ge:
                in_lfcs, avg_lfc_p, p_p_vals, avg_p_val, guide_list, ranks_p = g_p_val_P[g].split('_')
                in_lfcs, avg_lfc_n, n_p_vals, avg_n_val, guide_list, ranks_n = g_p_val_N[g].split('_')

                if float(avg_p_val) > float(avg_n_val):
                    avg_p_val = float(avg_p_val)
                    avg_lfc = float(avg_lfc_p)
                else:
                    avg_p_val = float(avg_n_val)
                    avg_lfc = float(avg_lfc_n)
                if include == 'N':
                    if len(guide_list.split(';')) != 1:
                        w.writerow((g, avg_lfc, avg_p_val, len(in_lfcs.split(';')), guide_list, in_lfcs, ranks_n, n_p_vals, ranks_p, p_p_vals))
                else:
                    w.writerow((g, avg_lfc, avg_p_val, len(in_lfcs.split(';')), guide_list, in_lfcs, ranks_n, n_p_vals, ranks_p, p_p_vals))
        if plot:
            val = plot_volcano(outputfile,min_pert,max_pert,label_num,c)
        counter += 1
    with open(o_folder+'/README.txt','w') as o:
        w = csv.writer(o,delimiter='\t')
        w.writerow((['Code Version: 2.3']))
        w.writerow((['Input file:'+inputfile]))
        w.writerow((['Chip file:'+args.chip_file]))
        w.writerow((['Output folder:'+o_folder]))
        w.writerow((['Fraction of perturbations for mean p-value, LFC:'+str(frac)]))
                    

