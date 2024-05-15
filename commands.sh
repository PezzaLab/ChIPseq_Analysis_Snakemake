################
# Note:
# To do each step copy all the lines (together) from that step
# (up to the next step and paste them on the command line. Steps 2 to 4
#  (included) are optional, 1 and 4 are mandatory.
################

# 1) load modules
ml slurm python/3.10.2 pandas/1.4.2

# 2) OPTIONAL. Dry run summary ("mock" run, no files are created)
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio -np --quiet

# 2.bis) OPTIONAL. Dry run extended version, text printed to the screen
# is also saved on file 'dry_run.txt'
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio -np | tee dry_run.txt

# 3) OPTIONAL. DAG (A plot to see the dependencies between the rules).
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio -np --rulegraph | \
dot -Tsvg > dag.svg


# 4) Real run (This will actually create the files3)
cd /s/pezzar-lab/{library_name} && \
sbatch --mail-type END,FAIL --job-name snkmk --mem 400 --wrap \
"snakemake --profile Config/Profiles/slurm_quio --notemp"

# -----------
# Info on snakemake pipeline used to process this library:
# Commit hash:
# {sha}
# Commit date:
# {commit_date}
