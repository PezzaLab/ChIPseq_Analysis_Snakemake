# Local VM Setup Guide (HPC Emulation)

This directory contains helper scripts to configure your Linux VM (Ubuntu/Debian) to mimic the HPC environment.

## 1. Set Up Single-Node SLURM
Run the SLURM setup script:
```bash
bash setup_vm/setup_slurm_vm.sh
```
This installs `slurmctld`, `slurmd`, and `munge`, configuring a single compute partition named `serial`.
Verify with:
```bash
sinfo
squeue
```

## 2. Set Up Environment Modules (Lmod)
Run the module generation script:
```bash
bash setup_vm/setup_modules_vm.sh
```
Then enable the module path:
```bash
module use ~/modulefiles
module avail
```
This makes all 15 modules required by `Config/config.yaml` (`samtools/1.14`, `bwa/0.7.15`, `macs2/2.2.7.1`, etc.) loadable via `ml <module>` or by Snakemake's `use-envmodules: true`.

## 3. Python & Snakemake Environment
On the VM, install the target Python dependencies from `target_env_info/pip_freeze.txt`:
```bash
pip install -r target_env_info/pip_freeze.txt
```
