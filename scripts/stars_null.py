'''
Creates null distribution for genetic perturbation screens at a user-defined threshold
Input: 1. Input file
       2. Chip File
       3. Threshold
       4. Number of iterations
       5. Use first perturbation value (Y/N)
       6. Max number of guides to use (best N)
edited by Peter Du (Apr 2018)
Updated by TD to take output arg instead of doing datetime
'''

import pandas as pd
from random import shuffle
from scipy import stats
import sys, csv, os, argparse
from math import log, log10

def getParser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=True)
    required = parser.add_argument_group('required')
    required.add_argument('--input-file',
        type=str,
        help='File containing spacer sequence',
        required=True)
    required.add_argument('--chip-file',
        type=str,
        help='Chip file',
        required=True)
    required.add_argument('--thr',
        type=int,
        help='Threshold',
        required=True)
    required.add_argument('--num-ite',
        type=int,
        help='Number of random iterations',
        required=True)
    parser.add_argument('--use-first-pert',
        type=str,
        default='Y',
        help='Y/N for value of first perturbation to be used in final STARS score calculation for a gene; default=Y')
    parser.add_argument('--max',
        type=int,
        default=-1,
        help='Use only the best N guides to generate the null file')
    parser.add_argument("--output-dir", 
	type=str, 
	help="output directory where results will be written")
    return parser

'''
Sort a data frame on the column specified and re-indexes
Argument: Dataframe, column to sort on, direction of sorting
Return Value: Sorted and re-index dataframe
'''
def sort_reindex(df, col, direction):
    if direction == 'P' or direction == 'p':
        df = df.sort_values([col], ascending=False)
        df.index = range(0, len(df))
    elif direction == 'N' or direction == 'n':
        df = df.sort_values([col], ascending=True)
        df.index = range(0, len(df))
    else:
        print 'Please enter a relevant direction; P for positive and N for negative'
        sys.exit(1)
    return df


if __name__ == '__main__':
    args = getParser().parse_args()
    inputfile = args.input_file
    max_pert = args.max
    st_in = pd.read_table(inputfile)
    st_in = pd.DataFrame(st_in.iloc[:,0])
    st_in.columns = ['Spacer Sequence']
    st_in['Score'] = range(1,len(st_in)+1)
    ref = pd.read_table(args.chip_file)
    ref_colnames = list(ref.columns)
    ref_colnames[0:2] = ['Spacer Sequence', 'Gene']
    ref.columns = ref_colnames
    direction = 'N'
    thr = args.thr
    num_ite = args.num_ite
    first_pert = args.use_first_pert
    out_dir = args.output_dir
    #os.path.dirname(inputfile)
    #print '\nOutput will be written to ' + out_dir
    outputfile= out_dir + '/Null_STARSOutput8_thr' + str(thr)
    if max_pert > 0:
        outputfile = outputfile + '_max' + str(max_pert)
    outputfile = outputfile + '.txt'
    with open(outputfile,'w') as o:
        w = csv.writer(o,delimiter='\t',lineterminator='\n')
        for rand in range(num_ite):
            print str(rand + 1) + '/' + str(num_ite)
            scores = list(st_in['Score'])
            shuffle(scores)
            st_in['Score'] = scores
            merged = pd.merge(st_in, ref, on='Spacer Sequence')
            merged = sort_reindex(merged, 'Score', direction)
            st_in = sort_reindex(st_in, 'Score', direction)
            grouped=merged.groupby('Gene')
            sps = list(st_in['Spacer Sequence'])
            tot_sps = len(sps)

            sp_rank = {}
            for i in range(tot_sps):
                sp_rank[sps[i]] = i+1
            thr_rank = (thr*tot_sps)/100
            ge=list(merged['Gene'])
            ge=list(set(ge))
            i=1
            print "Analyzing ..."
            for i,g in enumerate(ge):
                ge_df=grouped.get_group(g)
                gene = g
                spacers=list(ge_df['Spacer Sequence'])
                length_sps=len(spacers)
                rank=list()
                for s in spacers:
                    r=sp_rank[s]
                    rank.append(r)
                rank.sort()
                if first_pert == 'N' or first_pert == 'n':
                    rank=rank[1:]
                    j=2
                    index_add = 2
                elif first_pert == 'Y' or first_pert == 'y':
                    index_add = 1
                    j=1
                else:
                    print >> sys.stderr,'Please specify Y/y for value of first perturbation to be included else specify N/n.'
                    sys.exit(1)
                indices = [r for r in rank if r <= thr_rank]
                # remove until only best N indices
                indices.sort()
                if 0 < max_pert < length_sps:
                    while len(indices) > max_pert:
                        indices.pop(-1)
                if len(indices) == 0:
                    continue
                else:
                    sc_array = [stats.binom.pmf(indices.index(x)+index_add, length_sps, x/float(tot_sps)) for x in indices]
                    spacer_rank = sc_array.index(min(sc_array))+index_add
                    if min(sc_array) != 0:
                        log_score = -log10(min(sc_array))
                    else:
                        log_score = -log10(1E-10)
                    w.writerow([log_score])
