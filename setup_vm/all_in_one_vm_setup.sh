#!/usr/bin/env bash
set -euo pipefail

echo "=============================================================================="
echo ">>> [1/5] Setting up SLURM & Munge on the VM..."
echo "=============================================================================="
sudo apt-get install -y -qq slurmd slurmctld slurm-client munge lmod rsync git

# Configure Munge key if missing
if [ ! -s /etc/munge/munge.key ]; then
    sudo dd if=/dev/urandom of=/etc/munge/munge.key bs=1 count=1024
fi
sudo chown munge:munge /etc/munge/munge.key
sudo chmod 400 /etc/munge/munge.key
sudo systemctl enable munge
sudo systemctl restart munge

# Configure SLURM
HOSTNAME_VM=$(hostname -s)
NUM_CPUS=$(nproc)
TOTAL_MEM=$(free -m | awk '/^Mem:/{print $2}')
REAL_MEM=$(( TOTAL_MEM > 1000 ? TOTAL_MEM - 500 : TOTAL_MEM ))

sudo mkdir -p /etc/slurm /var/spool/slurmd /var/spool/slurmctld /var/log/slurm
sudo tee /etc/slurm/slurm.conf > /dev/null <<SLURM_CONF
ClusterName=local_hpc
SlurmctldHost=${HOSTNAME_VM}
MpiDefault=none
ProctrackType=proctrack/cgroup
ReturnToService=2
SlurmctldPidFile=/run/slurmctld.pid
SlurmctldPort=6817
SlurmdPidFile=/run/slurmd.pid
SlurmdPort=6818
SlurmdSpoolDir=/var/spool/slurmd
SlurmUser=slurm
StateSaveLocation=/var/spool/slurmctld
SwitchType=switch/none
TaskPlugin=task/affinity

# Timers
InactiveLimit=0
KillWait=30
MinJobAge=300
SlurmctldTimeout=120
SlurmdTimeout=300
Waittime=0

# Scheduling
SchedulerType=sched/backfill
SelectType=select/cons_tres
SelectTypeParameters=CR_Core_Memory

# Logging
SlurmctldDebug=info
SlurmctldLogFile=/var/log/slurm/slurmctld.log
SlurmdDebug=info
SlurmdLogFile=/var/log/slurm/slurmd.log
JobAcctGatherType=jobacct_gather/linux

# Node & Partition (matches pezzar-lab serial partition)
NodeName=${HOSTNAME_VM} CPUs=${NUM_CPUS} RealMemory=${REAL_MEM} State=UNKNOWN
PartitionName=serial Nodes=${HOSTNAME_VM} Default=YES MaxTime=INFINITE State=UP
SLURM_CONF

sudo chown -R slurm:slurm /var/spool/slurmd /var/spool/slurmctld /var/log/slurm /etc/slurm
sudo systemctl enable slurmd slurmctld
sudo systemctl restart slurmctld
sudo systemctl restart slurmd
sleep 2
sudo scontrol update NodeName="${HOSTNAME_VM}" State=RESUME 2>/dev/null || true

echo ""
echo "SLURM Cluster Status:"
sinfo

echo "=============================================================================="
echo ">>> [2/5] Setting up Lmod Environment Modules..."
echo "=============================================================================="
MODULE_DIR="${HOME}/modulefiles"
mkdir -p "${MODULE_DIR}"

create_mod() {
    local name="$1"
    local ver="$2"
    local bin="${3:-/usr/local/bin}"
    mkdir -p "${MODULE_DIR}/${name}"
    cat <<MOD_EOF > "${MODULE_DIR}/${name}/${ver}.lua"
whatis("Mock module for ${name} ${ver}")
prepend_path("PATH", "${bin}")
MOD_EOF
}

create_mod "bamutil" "1.0.15"
create_mod "bedtools" "2.30.0"
create_mod "bioconductor" "3.14"
create_mod "bwa" "0.7.15"
create_mod "cutadapt" "3.7"
create_mod "deeptools" "3.4.3"
create_mod "fastp" "0.23.2"
create_mod "homer" "5.1"
create_mod "macs2" "2.2.7.1"
create_mod "multiqc" "1.15"
create_mod "phantompeakqualtools" "1.2"
create_mod "picard" "2.21.2"
create_mod "R" "4.1.2-mkl"
create_mod "samtools" "1.14"
create_mod "sambamba" "0.8.2"

grep -qxF "module use ${MODULE_DIR}" "${HOME}/.bashrc" || echo "module use ${MODULE_DIR}" >> "${HOME}/.bashrc"
echo "Lmod modules ready in ${MODULE_DIR}"

echo "=============================================================================="
echo ">>> [3/5] Setting up Target Snakemake Python Environment..."
echo "=============================================================================="
CONDA_EXE="${HOME}/miniforge3/bin/conda"
MAMBA_EXE="${HOME}/miniforge3/bin/mamba"

if [ -x "${MAMBA_EXE}" ]; then
    SOLVER="${MAMBA_EXE}"
else
    SOLVER="${CONDA_EXE}"
fi

# Create target environment named snakemake_target if not exists
if ! "${CONDA_EXE}" env list | grep -q "snakemake_target"; then
    echo "Creating 'snakemake_target' conda environment..."
    ${SOLVER} create -y -n snakemake_target -c conda-forge -c bioconda \
        python=3.12 pandas numpy pyyaml gitpython snakemake
fi

# Install the exact SLURM executor plugin version (2.8.0)
"${HOME}/miniforge3/envs/snakemake_target/bin/pip" install --quiet \
    snakemake-executor-plugin-slurm==2.8.0 \
    snakemake-executor-plugin-slurm-jobstep==0.6.1

echo "Target environment packages:"
"${HOME}/miniforge3/envs/snakemake_target/bin/python" --version
"${HOME}/miniforge3/envs/snakemake_target/bin/snakemake" --version
"${HOME}/miniforge3/envs/snakemake_target/bin/pip" list | grep -i "snakemake-executor"

echo "=============================================================================="
echo ">>> [4/5] Pipeline Sync Status..."
echo "=============================================================================="
echo "Pipeline located at: ${HOME}/projects/ChIPseq_snakemake_pipeline"

echo "=============================================================================="
echo ">>> [5/5] Setup Complete!"
echo "=============================================================================="
