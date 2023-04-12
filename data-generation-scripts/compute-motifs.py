import generateKmerTargetDBlevels as gk

print('motif,length,ambiguity,n.genomes,n.species,species')

species_dict = {}
with open('../data/species-motifs.csv', 'r') as f:
	for i, line in enumerate(f.readlines()):
		if i==0:
			continue
		species, motif = line.strip().split(',')
		if motif in species_dict.keys():
			species_dict[motif].append(species)
		else:
			species_dict[motif] = [species] 

# dict of possible sequences for each k
possible_sequences = {4:0, 5:0, 6:0, 7:0, 8:0, 9:0}

for motif in species_dict.keys():
	possible_sequences[len(motif)] += len(gk.possibleSequences(motif))

print(possible_sequences)


