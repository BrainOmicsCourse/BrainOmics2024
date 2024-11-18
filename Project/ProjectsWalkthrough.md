# Projects walkthroughs

## Table of Contents

1. [Activate environment](#Activate-environment)

2. [Project notebooks](#Project-notebooks)

---------------

# Activate environment

In order to work with the notebooks for the project, you can activate the environment created for day 2. You can do this by running the following command in your terminal:

```bash
module load singularity/3.8.5
singularity exec -B /group/brainomics -H $HOME/BrainOmics2024/Project /group/brainomics/Container/Day2v2.sif jupyter lab -y --ip=0.0.0.0 --port=8888 --notebook-dir=$HOME/BrainOmics2024/Project
```

You can now copy the URL from the terminal (including the token) and paste it in your browser. You should now see the Jupyter Lab interface. 

# Project notebooks overview

The idea of the project is to apply what you've learned from day 1 and day 2 hands-on session on a new dataset. Specifically we will work with the Braun dataset, a single-cell RNA-seq dataset of developing human brain. You can have a look at the reference paper [here](https://doi.org/10.1126/science.adf1226). The data have already been downloaded and you can start playing with them with the Step1 notebook.

In the project folder you will find 2 notebooks:

* `BraunDataset_Step1.ipynb`, covering topics of day 1
* `BraunDataset_Step2.ipynb` covering partially the topics of day 2. You can explore more by yourself with the help of the trainers and of the hands-on session notebooks.

