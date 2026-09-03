#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Single-Node SLURM & Lmod Environment Setup for Linux VM
# Mimics HPC environment for testing the ChIP-seq / CUT&RUN Snakemake Pipeline
# ==============================================================================

echo ">>> [1/4] Installing SLURM, Munge, and Lmod..."
sudo apt-get update
sudo apt-get install -y slurmd slurmctld slurm-client munge lmod

# Generate munge key if missing
if [ ! -f /etc/munge/munge.key ]; then
    sudo create-munge-key
fi
sudo systemctl enable munge
sudo systemctl restart munge

echo ">>> [2/4] Configuring Single-Node SLURM..."
HOSTNAME_VM=$(hostname -s)
NUM_CPUS=$(nproc)
TOTAL_MEM=$(free -m | awk '/^Mem:/{print $2}')
REAL_MEM=$(( TOTAL_MEM > 1000 ? TOTAL_MEM - 500 : TOTAL_MEM ))

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

# TIMERS
InactiveLimit=0
KillWait=30
MinJobAge=300
SlurmctldTimeout=120
SlurmdTimeout=300
Waittime=0

# SCHEDULING
SchedulerType=sched/backfill
SelectType=select/cons_tres
SelectTypeParameters=CR_Core_Memory

# LOGGING & ACCOUNTING
SlurmctldDebug=info
SlurmctldLogFile=/var/log/slurm/slurmctld.log
SlurmdDebug=info
SlurmdLogFile=/var/log/slurm/slurmd.log
JobAcctGatherType=jobacct_gather/linux

# COMPUTE NODES & PARTITIONS
NodeName=${HOSTNAME_VM} CPUs=${NUM_CPUS} RealMemory=${REAL_MEM} State=UNKNOWN
PartitionName=serial Nodes=${HOSTNAME_VM} Default=YES MaxTime=INFINITE State=UP
SLURM_CONF

sudo mkdir -p /var/spool/slurmd /var/spool/slurmctld /var/log/slurm
sudo chown -R slurm:slurm /var/spool/slurmd /var/spool/slurmctld /var/log/slurm

echo ">>> [3/4] Starting SLURM daemons..."
sudo systemctl enable slurmd slurmctld
sudo systemctl restart slurmctld
sudo systemctl restart slurmd
sudo scontrol update NodeName="${HOSTNAME_VM}" State=RESUME || true

echo ">>> [4/4] Verifying SLURM status..."
sinfo
squeue

echo "=============================================================================="
echo "SLURM setup complete! You can now submit test jobs via 'sbatch' or 'sinfo'."
echo "=============================================================================="
