################
# Note:
# To do each step copy all the lines (together) from that step
# (up to the next step and paste them on the command line. Steps 2 to 4
#  (included) are optional, 1 and 4 are mandatory.
################

# 1) load modules
module purge && ml slurm
source /hpc-prj/pezza/conda/bin/activate
conda activate snakemake

# 2) OPTIONAL. Dry run summary ("mock" run, no files are created)
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio -np --quiet

# 2.bis) OPTIONAL. Dry run extended version, text printed to the screen
# is also saved on file 'dry_run.txt'
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio -np | tee dry_run.txt

# 3) OPTIONAL. DAG (A plot to see the dependencies between the rules).
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio --rulegraph | \
dot -Tsvg > rulegraph.svg


# 4) Real run (This will actually create the files)
cd /s/pezzar-lab/{library_name} && \
# Generate new terminal session
tmux new -s {library_name}
	# If already created one and want to re-join, use following command
tmux attach -t {library_name}
	# run command
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio --notemp

# ---------
# Helpful commands:
## Check running jobs
squeue --me -o %i%.60k

## Find log files
find -regex '.*log'

## See all available tmux sessions
tmux ls

## Kill a specific tmux session
tmux kill-session -t <session_name>

# -----------
# Info on snakemake pipeline used to process this library:
# Commit hash:
# {sha}
# Commit date:
# {commit_date}
