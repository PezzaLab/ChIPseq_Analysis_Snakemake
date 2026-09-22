################################################################################
# ChIP-seq / CUT&RUN Snakemake Pipeline Execution Commands
#
# HOW TO RUN:
# To do each step, copy all lines (together) from that step (up to the next step)
# and paste them into the command line.
# Steps 2 to 3 (included) are optional; steps 1 and 4 are mandatory.
#
# PIPELINE EXECUTION MODES:
# The pipeline supports two target execution modes:
#
# 1) 'basic' mode: Standard ChIP-seq / CUT&RUN Analysis
#    - Scope: delivers qctrl files, peaks and bigwigs (regular, doublestranded bigwigs).
#    - To run: Replace target 'full' with 'basic' in the commands below
#      (e.g., snakemake ... basic).
#
# 2) 'full' mode: Complete / Meiotic Hotspot & Single-Strand Analysis (Default)
#    - Scope: Everything in 'basic' mode PLUS hotspot analysis, single strand bigwigs, hotspot aggregate profiles and heatmaps and a couple more things.
#    - To run: Keep target 'full' in the commands below (e.g., snakemake ... full).
################################################################################

# 1) Load modules
module purge && ml slurm python/3.14.7

# 2) OPTIONAL. Dry run summary ("mock" run, no files are created)
# Note: target is 'full' by default; replace with 'basic' for basic analysis mode.
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio_repeat_10 -np --quiet full

# 2.bis) OPTIONAL. Dry run extended version, text printed to the screen
# is also saved on file 'dry_run.txt'
# Note: target is 'full' by default; replace with 'basic' for basic analysis mode.
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio_repeat_10 -np full | tee dry_run.txt

# 3) OPTIONAL. DAG (A plot to see the dependencies between the rules).
# Note: target is 'full' by default; replace with 'basic' for basic analysis mode.
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio_repeat_10 --rulegraph full | \
dot -Tsvg > rulegraph.svg


# 4) Real run (This will actually create the files)
cd /s/pezzar-lab/{library_name}
# Generate new terminal session
tmux new -s {library_name}
	# If you already created the session and want to re-join, use:
	# tmux attach -t {library_name}
	# Run command inside tmux (replace 'full' with 'basic' for basic analysis mode):
cd /s/pezzar-lab/{library_name} && \
snakemake --profile Config/Profiles/slurm_quio_repeat_10 --notemp full

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

##########################
### BACKUP IN DROPBOX ####
##########################
# Note: 'Results/checksums.tsv' is automatically generated at the end of the Snakemake
# pipeline run (cataloging MD5 checksums, sizes, and timestamps for all deliverables
# under Results/ exported below). It will be backed up automatically with the first command.
# To regenerate on demand:
#   bash workflow/Scripts/generate_rclone_checksums.sh .
# or:
#   snakemake --profile Config/Profiles/slurm_quio_repeat_10 checksums

# Load rclone
ml rclone

# Check what we are leaving behind:
rclone copy --dry-run /s/pezzar-lab/{library_name} dropboxOMRF:Bioinformatics/Libraries/{library_name} \
--min-size 50M 2>&1 | grep -E -v '.*bam|bw|fq\.gz|matrix|fastq\.gz.|\.sra*'

# Do backup
rclone copy /s/pezzar-lab/{library_name} dropboxOMRF:Bioinformatics/Libraries/{library_name} \
--max-size 50M \
--filter '+ .snakemake/slurm_logs/**' \
--filter '- .*' \
--filter '- .*/' \
--filter '- ~*' \
--filter '- *fastp.html' \
--filter '- not_copy*' \
--filter '- *.bai' \
--filter '- *.bam' \
--filter '- *.bw' \
--filter '- *.filename' \
--filter '- *.matrix' \
--filter '- *.fq.gz' \
--filter '- *.fastq.gz' \
--filter '- *.sra' \
--filter '- *.homer_anotated.tsv' \
--filter '- *.gappedPeak' \
--filter '- Intersect_HSs_plus_minus_2000_bp/' \
--filter '+ **blacklist**' \
--filter '- Peaks/**'

# Copy bigwigs
## Check if it is copying anything other than a .bw
rclone --dry-run copy /s/pezzar-lab/{library_name} dropboxOMRF:Bioinformatics/Libraries/{library_name} \
--include '*.bw' 2>&1 | grep -v '.*\.bw'

## Copy
rclone copy /s/pezzar-lab/{library_name} dropboxOMRF:Bioinformatics/Libraries/{library_name} \
--include '*.bw'
