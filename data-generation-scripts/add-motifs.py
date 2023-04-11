import generateKmerTargetDBlevels as gk 

with open('/well/shaw/users/amu125/projects/restriction-sites/db-outputs/motif-counts-useful.csv', 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i>0:
                for pos in gk.possibleSequences(line.split(',')[0]):
                        print(pos+','+line.strip('\n'))
