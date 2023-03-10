######
import sys
import os
# Variables go here
base_dir = os.getcwd()
database_name = sys.argv[1]
job_name = 'pan'+database_name[0:3]
######
print('#!/bin/bash')
print('#$ -cwd -V')
print('#$ -N %s -j y' % job_name)
print('#$ -q short.qc')
print('#$ -pe shmem 4')
print('')
print('echo "****************************************************"')
print('echo "SGE job ID: "$JOBID')
print('echo "SGE task ID: "$SGE_TASK_ID')
print('echo "Run on host: "`hostname`')
print('echo "Operating system: "`uname -s`')
print('echo "Username: "`whoami`')
print('echo "Started."')
print('echo "TIME: "`date`')
print('echo "****************************************************"')
print('')
print('')
print('# Requires')
print('# CONDA')
print('module load Anaconda3/2020.11')
print('eval "$(conda shell.bash hook)"')
print('conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env')
print('')
print('# And that the refseq files have been downloaded to $db_dir and are there as .fna (fasta)')
print('database_name=%s' % database_name)
print('short_name=$(python /well/shaw/users/amu125/programs/scratch/rescomp/makeShortName.py $database_name)')
print('base_dir=%s' % base_dir)
print('db_dir=$base_dir/Database_init')
print('prepare_out_dir=$base_dir/prepare_out')
print('annotate_out_dir=$base_dir/annotate_out')
print('pangenome_out_dir=$base_dir/pangenome_out')
print('core_out_dir=$base_dir/core_genome')
print('echo "****************************************************"')
print('echo "TIME: "`date`')
print('echo "PanACoTA: Preparing the genomes..."')
print('PanACoTA prepare --norefseq -o $prepare_out_dir -d $db_dir --threads 0 ')
print('echo "****************************************************"')
print('echo ""')
print('echo "****************************************************"')
print('echo "TIME: "`date`')
print('echo "PanACoTA: Annotating the genomes..."')
print('PanACoTA annotate --info $prepare_out_dir/LSTINFO-NA-filtered-0.0001_0.06.txt -n $short_name --prodigal --threads 0 --force -r $annotate_out_dir')
print('echo "****************************************************"')
print('echo " "')
print('echo "****************************************************"')
print('echo "TIME: "`date`')
print('echo "PanACoTA: Producing pangenome..."')
print('PanACoTA pangenome -l $annotate_out_dir/LSTINFO-LSTINFO-NA-filtered-0.0001_0.06.lst -n $short_name -d $annotate_out_dir/Proteins -o $pangenome_out_dir ')
print('echo "****************************************************"')
print('echo " "')
print('echo "****************************************************"')
print('echo "TIME: "`date`')
print('echo "PanACoTA: Producing core genome..."')
print('PanACoTA corepers -p $pangenome_out_dir/PanGenome-"$short_name".All.prt-clust-0.8-mode1.lst -o $core_out_dir -t 0.99 -X')
print('echo "****************************************************"')
print('echo " "')
print('echo "****************************************************"')
print('echo "TIME: "`date`')
print('echo "PanACoTA-output: Making core genome fastas..."')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/makeCoreGenomeFastas-with-plasmid-category.py --core_genome $core_out_dir/PersGenome_PanGenome-"$short_name".All.prt-clust-0.8-mode1.lst-all_0.99-mixed.lst --gene_dir $annotate_out_dir/Genes --out_dir $core_out_dir ')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/statsFastaFilesDir.py $core_out_dir/$short_name > $base_dir/genome_component_stats.txt')
print('echo "Done."')
print('echo "****************************************************"')
print('echo " "')
print('echo "****************************************************"')
print('echo "TIME: "`date`')
print('echo "RMES: Running RMES on the fastas..."')
print('for N in 0 2500 5000 10000 50000 100000; # make the subsampling dirs')
print('do')
print('    mkdir -p $base_dir/rmes-4-$N')
print('    mkdir -p $base_dir/rmes-5-$N')
print('    mkdir -p $base_dir/rmes-6-$N')
print('done')
print('echo "Made the directories"')
print('# Currently this all runs sequentially: suggest parallelising or splitting into array jobs to speed up for future')
print('cd $core_out_dir')
print('echo "Changed directory to" $core_out_dir')
print('if [ -f job-ids-rmes.txt ];')
print('then')
print('    rm job-ids-rmes.txt')
print('fi')
print('for genome in $(ls *core *accessory_chrom *accessory_plas );')
print('do')
print('    for N in 0 2500 5000 10000 50000 100000;')
print('    do    ')
print('        for k in $(seq 4 6);')
print('        do')
print('            echo $genome,$k,$N >> job-ids-rmes.txt')
print('        done')
print('    done')
print('done ')
print('')
print('# Now use this generated file to submit a default array job')
print('total_N=$(wc $core_out_dir/job-ids-rmes.txt -l | cut -d " " -f 1)')
print('echo "Expected number of files from RMES jobs:" $total_N')
print('# Check total number of jobs isnt >70k (scheduler limit is 75k)')
print('if [ $total_N -gt 70000 ];')
print('then')
print('    python /well/shaw/users/amu125/programs/scratch/rescomp/make_default_job.py "genome=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 1); k=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 2); N=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 3); conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env; echo "\$genome"; python /well/shaw/users/amu125/programs/scratch/rescomp/runRMES-with-subsampling.py --fasta  "\$genome" --output $base_dir/rmes-"\$k"-"\$N"/"\$genome".csv --bases_subsample "\$N" --k "\$k"" --array 1-70000 > rmes-job-1.sh')
print('    qsub rmes-job-1.sh -sync y')
print('    python /well/shaw/users/amu125/programs/scratch/rescomp/make_default_job.py "genome=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 1); k=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 2); N=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 3); conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env; echo "\$genome"; python /well/shaw/users/amu125/programs/scratch/rescomp/runRMES-with-subsampling.py --fasta  "\$genome" --output $base_dir/rmes-"\$k"-"\$N"/"\$genome".csv --bases_subsample "\$N" --k "\$k"" --array 70001-$total_N > rmes-job-2.sh')
print('    qsub rmes-job-2.sh -sync y')
print('fi')
print('if [ $total_N -lt 70000 ];')
print('then')
print('    python /well/shaw/users/amu125/programs/scratch/rescomp/make_default_job.py "genome=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 1); k=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 2); N=\$(sed -n "\${SGE_TASK_ID}"p $core_out_dir/job-ids-rmes.txt | cut -d "," -f 3); conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env; echo "\$genome"; python /well/shaw/users/amu125/programs/scratch/rescomp/runRMES-with-subsampling.py --fasta  "\$genome" --output $base_dir/rmes-"\$k"-"\$N"/"\$genome".csv --bases_subsample "\$N" --k "\$k"" --array 1-$total_N > rmes-job.sh')
print('    qsub rmes-job.sh')
print('fi')
print('echo "TIME: "`date`')
print('echo "Submitting RMES job on core/accessory genomes..."')
print('qsub rmes-job.sh -sync y')
print('')
print('# Want to ensure all those jobs have finished')
print('echo "Waiting for those jobs to finish..."')
print('# sleep 120 # Will give them two minutes')
print('sleep 60 # give the jobs a minute headstart before checking')
print('echo "Expecting $total_N files in total"')
print('n_rmes_output=0')
print('while [ $n_rmes_output -lt $total_N ]')
print('do ')
print('    sleep 10')
print('    n_rmes_output=$(find $base_dir/rmes-*/ -name "*.csv" | wc -l | cut -d " " -f 1)')
print('    echo "Found $n_rmes_output so far..."')
print('done')
print('echo "Found $n_rmes_output out of $total_N expected files"')
print('cd $base_dir')
print('echo "TIME: "`date`')
print('echo "Merging the RMES output..."')
print('# Merge the RMES outputs by accession into files for each k,N combination')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/make_default_job.py "conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env; for k in 4 5 6; do  for N in 0 2500 5000 10000 50000 100000;  do  python  /well/shaw/users/amu125/programs/scratch/rescomp/mergeRMES.py  --dir $base_dir/rmes-\$k-\$N --suffix \".csv\" --output $base_dir/rmes-k-\$k-\$N-all.csv; python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/rmes-\$k-\$N;  done; done" > merge-rmes-job.sh')
print('qsub merge-rmes-job.sh')
print('echo "Done"')
print('echo "****************************************************"')
print('echo "TIME: "`date`')
print('echo "****************************************************"')
print('')
print('# Move panacota files which we want to keep')
print('echo "Moving files..."')
print('mv $prepare_out_dir/LSTINFO-NA-filtered-0.0001_0.06.txt $base_dir/panacota-filtered-genomes.txt -f')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/getPanacotaGenomeFTP.py --input $base_dir/panacota-filtered-genomes.txt --output $base_dir/genome-ftps.txt')
print('mv $pangenome_out_dir/PanGenome-"$short_name".All.prt-clust-0.8-mode1.lst $base_dir/pangenome.lst -f')
print('mv $core_out_dir/PersGenome_PanGenome-"$short_name".All.prt-clust-0.8-mode1.lst-all_0.99-mixed.lst $base_dir/persistent_genome_t0.99_X.lst -f')
print('# Remove directories if they exist')
print('if [ -d $base_dir/Genes ]; then')
print('    rm -r $base_dir/Genes')
print('fi')
print('if [ -d $base_dir/Proteins ]; then')
print('    rm -r $base_dir/Proteins')
print('fi')
print('mv $annotate_out_dir/Genes $base_dir/ -f')
print('mv $annotate_out_dir/Proteins $base_dir/ -f')
print('')
print('# rm $base_dir/job-*.sh')
print('echo "****************************************************"')
print('echo "TIME: "`date`')
print('echo "Running rmsFinder..."')
print('cd $base_dir')
print('if [ -f $base_dir/Proteins/"$short_name".All.prt ];')
print('then')
print('    rm $base_dir/Proteins/"$short_name".All.prt')
print('fi')
print('')
print('ls $base_dir/Proteins/*.prt | rev | cut -d "/" -f 1 | rev > Proteins/proteomes.txt')
print('N_proteomes=$(wc -l Proteins/proteomes.txt | cut -d " " -f 1)')
print('mkdir -p $base_dir/rms-results ')
print('cd $base_dir/rms-results')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/make_default_job.py "conda activate /well/shaw/users/amu125/miniconda3/envs/panacota_rmes_env; genome=\$(sed -n "\${SGE_TASK_ID}"p $base_dir/Proteins/proteomes.txt); python /well/shaw/users/amu125/programs/rmsFinder/rmsFinder.py $base_dir/Proteins/"\$genome" --panacotafasta --output $base_dir/rms-results/"\$genome" --hmm tesson --db all" --name rmsFinder --array "1-"$N_proteomes"" > job-rmsFinder.sh')
print('qsub job-rmsFinder.sh -sync y')
print('cd $base_dir')
print('echo "****************************************************"')
print('echo "TIME: "`date`')
print('')
print('# Remove all the temporary files produced by panacota along the way')
print('echo "Removing files not needed at TIME: "`date`')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $prepare_out_dir')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $annotate_out_dir')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $pangenome_out_dir')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $core_out_dir')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/core_genomes')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/Genes')
print('echo "Expecting $N_proteomes RMS files in total"')
print('sleep 120 # give the rmsFinder jobs a 2-minute headstart before checking')
print('n_rms_output=0')
print('while [ $n_rms_output -lt $N_proteomes ]')
print('do ')
print('    sleep 30')
print('    n_rms_output=$(ls $base_dir/rms-results/*RMS.csv | wc -l | cut -d " " -f 1)')
print('    echo "Found $n_rms_output so far..."')
print('done')
print('echo "Found $n_rms_output out of $N_proteomes expected files"')
print('echo "Merging the rmsFinder output..."')
print('# Merge the rmsFinder RMS outputs')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/mergeRMS.py --dir $base_dir/rms-results --suffix "RMS" --out $base_dir/rmsFinder-RMS-merged.csv')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/mergeRMS.py --dir $base_dir/rms-results --suffix "MT" --out $base_dir/rmsFinder-MT-merged.csv')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/mergeRMS.py --dir $base_dir/rms-results --suffix "RE" --out $base_dir/rmsFinder-RE-merged.csv')
print('sleep 20')
print('echo "Removing final files no longer needed..."')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/rms-results')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/Proteins')
print('echo "Removing the initial fasta files..."')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $base_dir/Database_init')
print('python /well/shaw/users/amu125/programs/scratch/rescomp/rmDir.py --dir $core_out_dir')
print('echo "Finished at TIME:"`date`')
print('echo "****************************************************"')
