#!/bin/bash
#$ -cwd -V
#$ -N panEsc -j y
#$ -q short.qc
#$ -pe shmem 4

echo "****************************************************"
echo "SGE job ID: "$JOBID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started."
echo "TIME: "`date`
echo "****************************************************"


# Requires
# CONDA
module load Anaconda3/2020.11
eval "$(conda shell.bash hook)"
conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env

# And that the refseq files have been downloaded to $db_dir and are there as .fna (fasta)
database_name=Escherichia coli
short_name=$(python /well/shaw/users/amu125/programs/scratch/rescomp/makeShortName.py $database_name)
base_dir=/Users/Liam/Dropbox/_Projects/Restriction-sites/R-M-and-plasmids/scripts
db_dir=$base_dir/Database_init
prepare_out_dir=$base_dir/prepare_out
annotate_out_dir=$base_dir/annotate_out
pangenome_out_dir=$base_dir/pangenome_out
core_out_dir=$base_dir/core_genome
echo "****************************************************"
echo "TIME: "`date`
echo "PanACoTA: Preparing the genomes..."
PanACoTA prepare --norefseq -o $prepare_out_dir -d $db_dir --threads 0 
echo "****************************************************"
echo ""
echo "****************************************************"
echo "TIME: "`date`
echo "PanACoTA: Annotating the genomes..."
PanACoTA annotate --info $prepare_out_dir/LSTINFO-NA-filtered-0.0001_0.06.txt -n $short_name --prodigal --threads 0 --force -r $annotate_out_dir
echo "****************************************************"
echo " "
echo "****************************************************"
echo "TIME: "`date`
echo "PanACoTA: Producing pangenome..."
PanACoTA pangenome -l $annotate_out_dir/LSTINFO-LSTINFO-NA-filtered-0.0001_0.06.lst -n $short_name -d $annotate_out_dir/Proteins -o $pangenome_out_dir 
echo "****************************************************"
echo " "
echo "****************************************************"
echo "TIME: "`date`
echo "PanACoTA: Producing core genome..."
PanACoTA corepers -p $pangenome_out_dir/PanGenome-"$short_name".All.prt-clust-0.8-mode1.lst -o $core_out_dir -t 0.99 -X
echo "****************************************************"
echo " "
echo "****************************************************"
echo "TIME: "`date`
echo "PanACoTA-output: Making core genome fastas..."
python /well/shaw/users/amu125/programs/scratch/rescomp/makeCoreGenomeFastas-with-plasmid-category.py --core_genome $core_out_dir/PersGenome_PanGenome-"$short_name".All.prt-clust-0.8-mode1.lst-all_0.99-mixed.lst --gene_dir $annotate_out_dir/Genes --out_dir $core_out_dir 
python /well/shaw/users/amu125/programs/scratch/rescomp/statsFastaFilesDir.py $core_out_dir/$short_name > $base_dir/genome_component_stats.txt
echo "Done."
echo "****************************************************"
echo " "
echo "****************************************************"
echo "TIME: "`date`
echo "RMES: Running RMES on the fastas..."
for N in 0 2500 5000 10000 50000 100000; # make the subsampling dirs
do
    mkdir -p $base_dir/rmes-4-$N
    mkdir -p $base_dir/rmes-5-$N
    mkdir -p $base_dir/rmes-6-$N
done
echo "Made the directories"
# Currently this all runs sequentially: suggest parallelising or splitting into array jobs to speed up for future
cd $core_out_dir
echo "Changed directory to" $core_out_dir
if [ -f job-ids-rmes.txt ];
then
    rm job-ids-rmes.txt
fi
for genome in $(ls *core *accessory_chrom *accessory_plas );
do
    for N in 0 2500 5000 10000 50000 100000;
    do    
        for k in $(seq 4 6);
        do
            echo $genome,$k,$N >> job-ids-rmes.txt
        done
    done
done 

# Now use this generated file to submit a default array job
total_N=$(wc $core_out_dir/job-ids-rmes.txt -l | cut -d " " -f 1)
echo "Expected number of files from RMES jobs:" $total_N
# Check total number of jobs isnt >70k (scheduler limit is 75k)
if [ $total_N -gt 70000 ];
then
    python /well/shaw/users/amu125/programs/scratch/rescomp/make_default_job.py "genome=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 1); k=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 2); N=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 3); conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env; echo "\$genome"; python /well/shaw/users/amu125/programs/scratch/rescomp/runRMES-with-subsampling.py --fasta  "\$genome" --output $base_dir/rmes-"\$k"-"\$N"/"\$genome".csv --bases_subsample "\$N" --k "\$k"" --array 1-70000 > rmes-job-1.sh
    qsub rmes-job-1.sh -sync y
    python /well/shaw/users/amu125/programs/scratch/rescomp/make_default_job.py "genome=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 1); k=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 2); N=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 3); conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env; echo "\$genome"; python /well/shaw/users/amu125/programs/scratch/rescomp/runRMES-with-subsampling.py --fasta  "\$genome" --output $base_dir/rmes-"\$k"-"\$N"/"\$genome".csv --bases_subsample "\$N" --k "\$k"" --array 70001-$total_N > rmes-job-2.sh
    qsub rmes-job-2.sh -sync y
fi
if [ $total_N -lt 70000 ];
then
    python /well/shaw/users/amu125/programs/scratch/rescomp/make_default_job.py "genome=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 1); k=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 2); N=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 3); conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env; echo "\$genome"; python /well/shaw/users/amu125/programs/scratch/rescomp/runRMES-with-subsampling.py --fasta  "\$genome" --output $base_dir/rmes-"\$k"-"\$N"/"\$genome".csv --bases_subsample "\$N" --k "\$k"" --array 1-$total_N > rmes-job.sh
    qsub rmes-job.sh
fi
echo "TIME: "`date`
echo "Submitting RMES job on core/accessory genomes..."
qsub rmes-job.sh -sync y

# Want to ensure all those jobs have finished
echo "Waiting for those jobs to finish..."
# sleep 120 # Will give them two minutes
sleep 60 # give the jobs a minute headstart before checking
echo "Expecting $total_N files in total"
n_rmes_output=0
while [ $n_rmes_output -lt $total_N ]
do 
    sleep 10
    n_rmes_output=$(find $base_dir/rmes-*/ -name "*.csv" | wc -l | cut -d " " -f 1)
    echo "Found $n_rmes_output so far..."
done
echo "Found $n_rmes_output out of $total_N expected files"
cd $base_dir
echo "TIME: "`date`
echo "Merging the RMES output..."
# Merge the RMES outputs by accession into files for each k,N combination
python /well/shaw/users/amu125/programs/scratch/rescomp/make_default_job.py "conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env; for k in 4 5 6; do  for N in 0 2500 5000 10000 50000 100000;  do  python  /well/shaw/users/amu125/programs/scratch/rescomp/mergeRMES.py  --dir $base_dir/rmes-\$k-\$N --suffix ".csv" --output $base_dir/rmes-k-\$k-\$N-all.csv; python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/rmes-\$k-\$N;  done; done" > merge-rmes-job.sh
qsub merge-rmes-job.sh
echo "Done"
echo "****************************************************"
echo "TIME: "`date`
echo "****************************************************"

# Move panacota files which we want to keep
echo "Moving files..."
mv $prepare_out_dir/LSTINFO-NA-filtered-0.0001_0.06.txt $base_dir/panacota-filtered-genomes.txt -f
python /well/shaw/users/amu125/programs/scratch/rescomp/getPanacotaGenomeFTP.py --input $base_dir/panacota-filtered-genomes.txt --output $base_dir/genome-ftps.txt
mv $pangenome_out_dir/PanGenome-"$short_name".All.prt-clust-0.8-mode1.lst $base_dir/pangenome.lst -f
mv $core_out_dir/PersGenome_PanGenome-"$short_name".All.prt-clust-0.8-mode1.lst-all_0.99-mixed.lst $base_dir/persistent_genome_t0.99_X.lst -f
# Remove directories if they exist
if [ -d $base_dir/Genes ]; then
    rm -r $base_dir/Genes
fi
if [ -d $base_dir/Proteins ]; then
    rm -r $base_dir/Proteins
fi
mv $annotate_out_dir/Genes $base_dir/ -f
mv $annotate_out_dir/Proteins $base_dir/ -f

# rm $base_dir/job-*.sh
echo "****************************************************"
echo "TIME: "`date`
echo "Running rmsFinder..."
cd $base_dir
if [ -f $base_dir/Proteins/"$short_name".All.prt ];
then
    rm $base_dir/Proteins/"$short_name".All.prt
fi

ls $base_dir/Proteins/*.prt | rev | cut -d "/" -f 1 | rev > Proteins/proteomes.txt
N_proteomes=$(wc -l Proteins/proteomes.txt | cut -d " " -f 1)
mkdir -p $base_dir/rms-results 
cd $base_dir/rms-results
python /well/shaw/users/amu125/programs/scratch/rescomp/make_default_job.py "conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env; genome=\$(sed -n "\${SGE_TASK_ID}"p $base_dir/Proteins/proteomes.txt); python /well/shaw/users/amu125/programs/rmsFinder/rmsFinder.py $base_dir/Proteins/"\$genome" --panacotafasta --output $base_dir/rms-results/"\$genome" --hmm tesson --db all" --name rmsFinder --array "1-"$N_proteomes"" > job-rmsFinder.sh
qsub job-rmsFinder.sh -sync y
cd $base_dir
echo "****************************************************"
echo "TIME: "`date`

# Remove all the temporary files produced by panacota along the way
echo "Removing files not needed at TIME: "`date`
python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $prepare_out_dir
python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $annotate_out_dir
python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $pangenome_out_dir
python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $core_out_dir
python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/core_genomes
python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/Genes
echo "Expecting $N_proteomes RMS files in total"
sleep 120 # give the rmsFinder jobs a 2-minute headstart before checking
n_rms_output=0
while [ $n_rms_output -lt $N_proteomes ]
do 
    sleep 30
    n_rms_output=$(ls $base_dir/rms-results/*RMS.csv | wc -l | cut -d " " -f 1)
    echo "Found $n_rms_output so far..."
done
echo "Found $n_rms_output out of $N_proteomes expected files"
echo "Merging the rmsFinder output..."
# Merge the rmsFinder RMS outputs
python /well/shaw/users/amu125/programs/scratch/rescomp/mergeRMS.py --dir $base_dir/rms-results --suffix "RMS" --out $base_dir/rmsFinder-RMS-merged.csv
python /well/shaw/users/amu125/programs/scratch/rescomp/mergeRMS.py --dir $base_dir/rms-results --suffix "MT" --out $base_dir/rmsFinder-MT-merged.csv
python /well/shaw/users/amu125/programs/scratch/rescomp/mergeRMS.py --dir $base_dir/rms-results --suffix "RE" --out $base_dir/rmsFinder-RE-merged.csv
sleep 20
echo "Removing final files no longer needed..."
python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/rms-results
python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/Proteins
echo "Removing the initial fasta files..."
python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/Database_init
python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $core_out_dir
echo "Finished at TIME:"`date`
echo "****************************************************"
