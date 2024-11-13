## Table of Contents

1. [Activate environment for Brainomics day 2](#1-activate-environment-for-brainomics-day-2)

2. [Explore the course material and start the notebook](#2-explore-the-course-material-and-start-the-notebook)

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
singularity exec -B /group/brainomics/ -H $HOME/BrainOmics2024/2_Day2  /group/brainomics/Container/Day2v2.sif jupyter lab -y --ip=0.0.0.0 --port=8888 --notebook-dir=$HOME/BrainOmics2024/2_Day2
```

Copy the URL from the terminal (including the token) and paste it in your browser. You should now see the Jupyter Lab interface.

<br> </br>

## 2. Explore the course material and start the notebook

On day 2 of Brainomics we will explore advanced methods for single-cell transcriptomics data analysis, including metacells, trajectory inference and transcription factor activity computation. Within the folder of day 2 you will find 4 notebooks:

* `1_Start.ipynb`: dedicated to stringent QC filter, multiplets detection and the creation of metacells
* `2_Palantir.ipynb`: dedicated to metacells processing and trajectory detection
* `3_Trajectories`: dedicated to trajectories refinement and gene patterns characterization along the trajectories themselves
* `4_Decoupler`: dedicated to the computation of transcription factor activity in combination with previously inferred trajectories

Additionally there will be:

* `utils`: a folder containing scripts with utility functions and files
* `Resources.md`: a markdown file with references to the libraries and data used in the notebooks
* `Compiled`: a folder containing the notebooks with outputs both in `ipynb` and `html` format 
