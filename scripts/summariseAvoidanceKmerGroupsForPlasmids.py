
import generateKmerTargetDBlevels as gk
import calculateAvoidanceKmerGroupsAcrossDatasetv2 as ca
import pandas as pd
import pickle
import argparse

def get_options():
	parser = argparse.ArgumentParser(description='Summarise output for Redondo-Salvo plasmids by k-mer groupings.',
                                  prog='summariseAvoidanceKmerGroupsForPlasmids')
	parser.add_argument('--outdir', help='output directory', required=True)
	parser.add_argument('--targetdict', help='target dictionary', required=True)
	parser.add_argument('--k', help='k-mer length', required=True)
	parser.add_argument('--n', help='subsampling', required=True)
	parser.add_argument('--inclusive', help='whether kmer groups are inclusive', action='store_true', default=False)
	return parser.parse_args()



#target_dict = '/well/shaw/users/amu125/projects/restriction-sites/db-outputs/5-species-target-dict.pkl'
#inclusive = True

#with open(target_dict, 'rb') as handle:
#        species_kmer_dict = pickle.load(handle)
#print(species_kmer_dict['Escherichia_coli'])

#id_2_species_dict = {}
#with open('/well/shaw/users/amu125/data/redondo-salvo/id_2_species_filtered.csv') as f:
#	for line in f.readlines():
#		plasmid, species = line.strip().split(',')
#		id_2_species_dict[plasmid] = species

# Strip out non-E. faecalis
#efaec_dict = {k:v for k, v in id_2_species_dict.items() if v=='Escherichia_coli'}
#id_2_species_dict = efaec_dict

def summariseResults(big_df, output_file, species_dict, inclusive=False):
	with open(output_file, 'w') as output_f:
		output_f.write('plasmid,species,kmer_category,n,rank,score\n')
		for each_species in set(id_2_species_dict.values()):
			groups_dict = {group: [] for group in gk.KMER_CATEGORIES}
			kmer_dict = species_dict[each_species]
			if inclusive==False:
				for group in groups_dict.keys():
					for kmer in kmer_dict[group]:
						groups_dict[group].append(kmer.lower()) # need lower because of rmes results
			elif inclusive==True:
				_ = [groups_dict['Palindromic'].append(kmer.lower()) for kmer in kmer_dict['Palindromic']] # need lower due to rmes results
				cumulative_kmers = kmer_dict['Species']
				_ = [groups_dict['Species'].append(kmer.lower()) for kmer in cumulative_kmers]
				for group in ['Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom']:
					cumulative_kmers += kmer_dict[group]
					_ = [groups_dict[group].append(kmer.lower()) for kmer in list(cumulative_kmers)] # need lower due to rmes 
			print(groups_dict)
			local_dict = {k:v for k, v in id_2_species_dict.items() if v==each_species}
			for plasmid, species in local_dict.items():
				small_df = big_df[big_df['plasmid'] == plasmid] 
			
				print(plasmid, species)
				results = ca.summariseRmes(small_df, groups_dict)
				for i, category in enumerate(gk.KMER_CATEGORIES):
					output_f.write(','.join([str(x) for x in [plasmid,species,category,results[i][0], results[i][1], results[i][2]]])+'\n')	

id_2_species_dict = {}
with open('/well/shaw/users/amu125/data/redondo-salvo/id_2_species_filtered.csv') as f:
	for line in f.readlines():
		plasmid, species = line.strip().split(',')
		id_2_species_dict[plasmid] = species



def main():


	args = get_options()
	target_dict = args.targetdict
	with open(target_dict, 'rb') as handle:
		species_kmer_dict = pickle.load(handle)
	output_dir = args.outdir
	inclusive = args.inclusive
	subsampling = str(args.n)
	k = str(args.k)
#big_df = pd.read_csv('/well/shaw/users/amu125/projects/restriction-sites/data/plasmids/rmes/redondo-salvo-k-6.csv')
#summariseResults(big_df, 'redondo-salvo-summary-n-0.csv', inclusive=True)

	print('SUBSAMPLING:'+str(subsampling))
	big_df = pd.read_csv('/well/shaw/users/amu125/projects/restriction-sites/data/plasmids/rmes/redondo-salvo-k-'+k+'-n-'+subsampling+'.csv')
	summariseResults(big_df, output_dir+'/'+'redondo-salvo-summary-k-'+k+'-n-'+str(subsampling)+'-'+str(inclusive)+'.csv', species_kmer_dict, inclusive=inclusive)

if __name__=="__main__":
	main()
