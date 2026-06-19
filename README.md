# GTRspmix
**GTRspmix** is a protein mixture model with multiple exchangeability matrices (GTRs) and profiles.  
Profiles are clustered into several groups. Then each group has own GTR matrix.

Details are shown in our paper. Please cite it if you find this useful.  

> **Ryo Harada _et. al._ (2026). GTRspmix: Capturing Heterogeneity of Exchangeability Across Sites to Improve Protein Phylogenetics.** _bioRxiv_. [doi.org/10.64898/2026.06.18.729217](https://doi.org/10.64898/2026.06.18.729217)

**Note**: We are planning to update scripts to use real EM instead of approximation, and implementing GTRspmix estimation in IQ-TREE.

---


## Table of Contents
1. [Overview](#1-overview)
2. [Dependencies & Installation](#2-dependencies--installation)
3. [Usage Guide](#3-usage-guide)
    - [Preparation](#31-preparation)
    - [Choose your running mode](#32-choose-your-running-mode)
    - [Select Parameters to be Optimized](#33-select-parameters-to-be-optimized)
4. [Core Output Files](#4-core-output-files)
5. [Command Examples](#5-command-examples)
6. [Licence](#6-licence)


## 1. Overview
This script iteratively optimizes GTRspmix model parameters by repeating specific steps (EM approximation). Use `--opt-gtr` and/or `--opt-profile-FO`/`--opt-profile-F` to enable the corresponding optimization phases.

Each iteration consists of two main phases. Each phase is only executed if its corresponding flag is enabled.

| Phase | Step | Optimized params | Fixed params | Required flags |
|:--:|:--:|:--:|:--:|:--:|
| GTR Optimization | 1 | $\theta$ | GTR, Profile | `--opt-gtr` |
|  | 2 | GTR | $\theta$, Profile | `--opt-gtr` |
| Profile Optimization | 3 | $\theta$ | GTR, Profile | `--opt-profile-FO` or `--opt-profile-F` |
|  | 4 | Profile | $\theta$, GTR | `--opt-profile-FO` or `--opt-profile-F` |

- $\theta$ is a parameter set including branch lengths, weights, and rates.
- If only `--opt-gtr` is set: Step 3 and 4 are skipped. The script only optimizes GTR matrices and global parameters.
- If only `--opt-profile-FO` or `--opt-profile-F` is set: Step 1 and 2 are skipped. The script optimizes only frequency profiles and global parameters.
- If `--opt-gtr` and either `--opt-profile-FO` or `--opt-profile-F` are set: All steps (1–4) are executed in each iteration.

**Profile optimization is experimental and may be unstable.** 

For optimizing GTR and profile parameters, soft-partitioning is used to approximate EM algorithm.
We generate a sub-alignment for each profile cluster or profile class. 
Each site is allocated round($p \times b$) times to the sub-alignment, where $p$ is a posterior probability for the class and $b$ is a scaling factor.


---

## 2. Dependencies & Installation
### 2.1. Download GTRspmix files
Please download the source code from the latest release in this repository.
Unzip the downloaded archive and navigate to the directory:

### 2.2. Software Requirements
GTRspmix script has been tested with following specific versions:
- Python v3.11.15
- [IQ-TREE](https://github.com/iqtree/iqtree3) v2.3.4/v3.0.1/v3.1.1 
- [gotree](https://github.com/evolbioinfo/gotree) v0.4.5 

### 2.3.a. Installation via Singularity (Recommended)
[Singularity](https://github.com/sylabs/singularity) is highly recommended for HPC clusters to ensure the reproducibility of the environment. The Singularity image includes both IQ-TREE and gotree.
```
sudo singularity build gtrspmix.sif singularity.def
singularity run gtrspmix.sif -h
```

### 2.3.b. Manual Installation
If you prefer not to use Singulairty, you have to install IQ-TREE and gotree on your system.  
Then python requirements will be installed by 
```
pip install -r requirements.txt
```
If IQ-TREE and gotree are not available in your system `$PATH`, you must specify their absolute paths using the following flags during execution:
- `--iqtree /path/to/iqtree3`
- `--gotree /path/to/gotree`

---

## 3. Usage Guide

### 3.1. Preparation
In addition to alignment and guide tree, GTRspmix requires an initial profile mixture model.
- Custom Profiles: We recommend using [MaMMAL](https://github.com/TheBrownLab/mammal) or [MEOW](https://github.com/RogerLab/meow) to generate starting frequency profiles.

Or

- Empirical Models: Standard empirical models (C60 series and SXXCYY series) are pre-installed. You can specify these names (e.g., `--nexus-few C10 --nexus-many C60`, `--nexus C60`, and `--model S10C60`) without providing external files.

### 3.2. Choose your running mode
Select a mode based on your starting input.

| Mode | Input Requirements | Use Case |
|:--:|:--:|:--|
| FromScratch<br> (SPPC) | `--nexus-few` & <br> `--nexus-many` | Cluster a large profile set (e.g. MEOW60) into fewer groups (e.g. 10). | 
| FromScratch<br> (K-means) | `--nexus` & <br> `-km` | Directly cluster profiles using K-means. | 
| ReStart | `--nexus` & <br> `--json` | Resume an interrupted run using a model and clustering JSON from previous run. |
| PreDefined | `--model` | Fine-tune an empirical GTRspmix model (e.g., S10pfamC60). |

### 3.3. Select Parameters to be Optimized
You can control which parameters are optimized using the following flags:

- `--opt-gtr`: Optimize GTRs.  
- `--opt-profile-FO`: Optimize frequency profiles by ML estimation. (Note: This function is **Experimental**)
- `--opt-profile-F`: Optimize frequency profiles by weighted observed frequency. (Note: This function is **Experimental**)
- `--opt-gtr` and either `--opt-profile-FO` or `--opt-profile-F`: Performs full optimization of GTRspmix model parameters.

---

## 4. Core Output Files
| File | Description |
|:--:|:--:|
| `model_best.nex` | The final optimized model. |
| `GTRspmix_maker.log` | The main execution log. Check this file to monitor the optimization progress, specifically the $lnL$ values for each iteration. | 
| `d_cluster.json` | Cluster mapping information. This JSON file records which Profile classes belong to which GTR cluster. Essential for restarting runs or manually inspecting the model structure. | 
| `model_X.nex` | The model file generated at iteration $X$. Essential for restarting runs. | 


---


## 5. Command Examples
Detailed flags can be found via `gtrspmix.py --help`. Below are typical command examples.

### Optimize 10 GTRs with MEOW60 profiles (SPPC clustering)
```bash
gtrspmix.py \
--opt-gtr \
-s alignment.fasta \
-t guide.treefile \
--nexus-few meow_10.nex \
--nexus-many meow_60.nex \
-m-gtr20 ELM \
-m-rate G4 \
--scale-gtr 10 \
-me-theta 0.01 \
-me-gtr 0.99 \
-me 10 \
-nt 8 \
-mem 100G \
-o GTRspmix_out
```

### Optimize 10 GTRs with MEOW60 profiles (Kmeans clustering)
```bash
gtrspmix.py \
--opt-gtr \
-s alignment.fasta \
-te guide.treefile \
--kmeans 10 \
--nexus meow_60.nex \
-m-gtr20 ELM \
-m-rate G4 \
--scale-gtr 10 \
-me-theta 0.01 \
-me-gtr 0.99 \
-me 10 \
-nt 8 \
-mem 100G \
-o GTRspmix_out
```

### Optimize 10 GTRs with MEOW60 profiles (Restarting)
When you restart runs, please copy latest nexus file and d_cluster.json file. Then specify new output directory.
```bash
gtrspmix.py \
--opt-gtr \
-s alignment.fasta \
-te guide.treefile \
--json d_cluster.json \
--nexus model_best.nex \
-m-rate G4 \
--scale-gtr 10 \
-me-theta 0.01 \
-me-gtr 0.99 \
-me 10 \
-nt 8 \
-mem 100G \
-o GTRspmix_out_restart
```

### Fine-ture S10pfamC60 model
```bash
gtrspmix.py \
--opt-gtr \
-s alignment.fasta \
-te guide.treefile \
--model S10pfamC60 \
-m-rate G4 \
--scale-gtr 10 \
-me-theta 0.01 \
-me-gtr 0.99 \
-me 10 \
-nt 8 \
-mem 100G \
-o GTRspmix_out
```

### Optimize 10 GTRs and 60 profiles (SPPC)
```bash
gtrspmix.py \
--opt-gtr \
--opt-profile-FO \
-s alignment.fasta \
-t guide.treefile \
--nexus-few meow_10.nex \
--nexus-many meow_60.nex \
-m-gtr20 ELM \
-m-rate G4 \
--scale-gtr 10 \
--scale-profile 100 \
-me-theta 0.01 \
-me-gtr 0.99 \
-me-pro 0.01 \
-me 10 \
-nt 8 \
-mem 100G \
-o GTRspmix_out
```

### Optimize only profiles
```bash
gtrspmix.py \
--opt-profile-FO \
-s alignment.fasta \
-te guide.treefile \
--nexus meow_60.nex \
-m-gtr20 ELM \
-m-rate G4 \
--scale-profile 100 \
-me-theta 0.01 \
-me-pro 0.01 \
-me 10 \
-nt 8 \
-mem 100G \
-o GTRspmix_out
```

---

## 6. Licence
This project is licensed under the GPL-3.0 License. See the `LICENSE` file for details.






