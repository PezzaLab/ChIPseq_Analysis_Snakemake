################
# Note:
# Please replace the words between curly braces with the aproppiate value 
# (curly braces should disappear as well)
################

# load modules
ml slurm python/3.10.2 pandas/1.4.2


# Dry run summary ("mock" run, no files are created)
cd /s/Pezza/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio -np --quiet
# Dry run extended version, text printed to the screen is also
# saved on file 'dry_run.txt'
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio -np | tee dry_run.txt


# DAG (A plot to see the dependencies between the rules)
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio -np --rulegraph | \
dot -Tsvg > dag.svg


# Real run (This will actually create the files3)
cd /s/pezzar-lab/{library_name} && \
sbatch --mail-user agustin-carbajal@omrf.org --job-name snkmk --mem 400 --wrap \
"snakemake --profile Config/Profiles/slurm_quio --notemp"
