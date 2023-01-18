import generateKmerTargetDB as gk

print('motif,length,ambiguity,n.genomes,n.species,species')

species_dict = {}
with open('/well/shaw/users/amu125/projects/restriction-sites/db-outputs/species-motifs.csv', 'r') as f:
	for line in f.readlines():
		species, motif = line.strip().split(',')
		if motif in species_dict.keys():
			species_dict[motif].append(species)
		else:
			species_dict[motif] = [species] 

with open('/well/shaw/users/amu125/projects/restriction-sites/db-outputs/motif-counts.csv', 'r') as f:
	for line in f.readlines():
		count, motif = line.strip().split(',')
		motif_length = len(motif)	
		ambiguity_count = len(gk.possibleSequences(motif))
		species_count = len(species_dict[motif])
		species_list = '-'.join(species_dict[motif])
		print(motif+','+str(motif_length)+','+str(ambiguity_count)+','+str(count)+','+str(species_count)+','+species_list)
