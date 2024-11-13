## Table of Contents

1. [Activate environment for Brainomics day 3](#1-activate-environment-for-brainomics-day-3)

2. [Explore the course material and start the notebook](#2-explore-the-course-material-and-start-the-notebook)

3. [Solutions to the exercises](#3-solutions-to-the-exercises)

---------------

<br> </br>

## 1. Activate environment for Brainomics day 3

* You should already have cloned the Brainomics repository in the previous days. If you haven't done it yet, you can do it by running the following command:

```bash
git clone https://github.com/BrainOmicsCourse/2024.git
```

For day 3 we have a new container for the specific tasks of spatial transcriptomics. We will use the folder `day3` as the home of the container. The __path to the container__ will be shared with you in the chat and on the whiteboard. Open the terminal and execute the following commands, after modifying the path to the shared folder:

```bash
module load singularity/3.8.5
singularity exec -B /group/brainomics/ -H $HOME/BrainOmics2024/3_Day3  /group/brainomics/Container/courses_brainomics2_day3-0.0.1.sif /bin/bash
```

This should open a new terminal inside the container. You need to initialize the conda environment by running the following commands:

```bash
conda init
source ~/.bashrc
conda-env list # you should see a list of environments, including day3_env
conda activate day3_env
```

Now we activated our environment so we can start Jupyter Lab:

```bash
jupyter lab
```

Copy the URL from the terminal (including the token) and paste it in your browser. You should now see the Jupyter Lab interface.

<br> </br>

## 2. Explore the course material and start the notebook

On day 3 of Brainomics we will explore methods for spatial transcriptomics data analysis. Within the folder of day 3 you will find:

* `spatial_day3.ipynb`, the Jupyter notebook on which we will work today.
* two scripts `.py`, with additional code for the preprocessing of some of the data used in the notebook.
* folders containing the scripts of two of the packages will be using today: `SpaGE` and `Bansky`. You can find updates of the packages in the respective repositories: [SpaGE](https://github.com/tabdelaal/SpaGE) and [Bansky](https://github.com/prabhakarlab/Banksy).
* folders containing utility functions
* a folder containing the notebook with the solutions to the exercise, both as `.ipynb` and compiled in html format.

You can open the notebook as we will be working mainly on it today. The notebook is divided in 4 sections:

* Data loading and explorative data analysis
* Prediction of spatial patterns for unmeasured genes using scRNA-seq data with SpaGE
* Imputation of unmeasured genes with Tangram
* Spatial clustering with Bansky

## 3. Solutions to the exercises

Solutions to the exercise are in the shared folder, you can copy them in day 3 folder by copying this command in the terminal:

```bash
rsync -cvr /group/brainomics/InputData/day3/Compiled $HOME/BrainOmics2024/3_Day3
```

You should now see the Compiled folder in the 3_Day3 folder. There's a compiled .ipynb and .html file with the solutions to the exercises.