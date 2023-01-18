# Combines the output of run-rmsFinder-Redondo-Salvo.sh
import glob
import re
import sys
import calculateAvoidanceKmerGroupsAcrossDataset as ca
import csv

def combineFiles(directory, pattern):
    '''combines all files matching a pattern''' 
    files = glob.glob(directory+pattern)    
    all_entries = []
    for filename in files:
        plasmidname = re.sub('_%s' % pattern, '', re.sub('.*/', '', filename))        
        with open(filename, 'r') as f:
            _ = f.readline # ignore first line
            new_entries = [x.strip('\n') for x in f.readlines()]
            if new_entries[0]!='':
                new_entries = [','.join([plasmidname, x]) for x in new_entries]
                _ = [all_entries.append(x) for x in new_entries[1:]]
    return(all_entries)

output_dir="/well/shaw/users/amu125/projects/restriction-sites/data/plasmids/rmsFinder/"

rms_header = 'plasmid,sequence,contig,contig_description,pos_MT,pos_RE,prot_MT,prot_RE,hit_MT,sim_MT,hit_RE,sim_RE'
mt_header = 'plasmid,qseqid,sseqid,pident,length,qlen,evalue,coverage_threshold_met,target,similarity,hit_type,n_REBASE_hits,identicalTarget'
re_header = mt_header

rms_files = combineFiles(output_dir, '*RMS.csv')
print('done RMS')
mt_files = combineFiles(output_dir, '*MT.csv')
print('done MT')
re_files = combineFiles(output_dir, '*RE.csv')
print('done RE')
with open(output_dir+'/RMS_combined.csv', 'w') as f:
    f.write(rms_header+'\n')
    for rms in rms_files:
        f.write('%s\n' % rms)

with open(output_dir+'/MT_combined.csv', 'w') as f:
    f.write(mt_header+'\n')
    for mt in mt_files:
        split_mt = list(csv.reader([mt]))[0] 
        print(split_mt)
        target = split_mt[8]
        target_list = [re.sub(',', '', x) for x in target.split()]  # check for more than one target
        for target in target_list:
            possible_targets = ca.possibleSequences(target)
            for possible_target in possible_targets:
                possible_target_string = ','.join([','.join(split_mt[0:8]), possible_target, ','.join(split_mt[9:])])
                f.write('%s\n' % possible_target_string)

with open(output_dir+'/RE_combined.csv', 'w') as f:
    f.write(re_header+'\n')
    for mt in re_files:
        f.write('%s\n' % re)


