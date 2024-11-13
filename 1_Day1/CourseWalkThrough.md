__SECTIONS__

- [1. Clone the course material](#1-clone-the-course-material)
    + [Connect to the VDI and login to linux VM](#connect-to-the-vdi-and-login-to-linux-vm)
    + [Clone the repo](#clone-the-repo)
- [2. Launch the analysis container](#2-launch-the-analysis-container)
  * [Load Singularity Module](#load-singularity-module)
  * [Launch the container](#launch-the-container)
- [3. Explore the course material in jupyter](#3-explore-the-course-material-in-jupyter)
- [4. First notebook from raw data to dimensionality reduction](#4-first-notebook-from-raw-data-to-dimensionality-reduction)
- [5. Second notebook cluster characterization](#5-second-notebook-cluster-characterization)
---------------

<br> </br>

## 1. Clone the course material

* As a first step, we will clone the repo containing the material for this course day.
* After the VM login, we will open a terminal window and, from the course folder, clone the repo.


#### Connect to the VDI and login to linux VM


#### Clone the repo

Clone the repo:

```bash
cd $HOME
git clone https://github.com/BrainOmicsCourse/BrainOmics2024.git
```
--------

<br> </br>

## 2. Launch the analysis container

### Load Singularity Module

```bash
module load singularity/3.8.5
```
### Launch the container

```bash
singularity exec -B <PATH/TO/THE/SHARED/FOLDER> -H $HOME/BrainOmics2024/1_Day1/ \
    <PATH/TO/THE/SHARED/FOLDER>/Container/Day1_container.sif \
    jupyter lab -y --ip=0.0.0.0 --port=8888 --notebook-dir=$HOME/BrainOmics2024/1_Day1/
```

---------

<br> </br>

## 3. Explore the course material in jupyter

* Our tutorial will focus on key steps of scRNASeq analysis taking as example a publicly available dataset: scRNASeq profiling of human fetal cortex at mid-gestation.
* Beyond the jupyter notebooks, in this repo you can find the following material:

#### Dataset information

* [__Dataset information__](Resources.md)
* [__'AssembleAdata'__](Compiled/0_AssembleAdata.html) assembly html file: it reports the first step of the analytical workflow, to assemble the starting AnnData. While we will not perform this part in practice, you can look at the workflow by opening in the web browser.  

#### [Helper file](HelperFunctions/Day1Helper.py)
* It contains custom functions that are used in the notebooks.
* Keeping custom function in a separate .py file makes the flow more organized and the functions easier to maintain.


#### [Resource file](Resources.md)
It contains links to relevant external resources and material.

-------------

<br> </br>

## 4. First notebook from raw data to dimensionality reduction

* In the __Day1__ folder, you find the [__'1_FiltNormBatch'__](1_FiltNormBatch.ipynb) notebook structured in _markdown cells_ that contains information and comments and _code cells_ with the actual code chunks.
* Starting from the assembled anndata, we will work through the following analytical steps: quality control; filtering; normalization; feature selection; cluster identification.
* If you get lost, in the 'Compiled' folder you will find the [__notebook version__](Compiled/1_FiltNormBatch.ipynb) compiled with the outcome of each code chunk


-----------

<br> </br>

## 5. Second notebook cluster characterization

* Always in the __Day1__ folder, you find the [__'2_Clusters'__](2_Clusters.ipynb) notebook.
* Starting from the output of the previous notebook, we will characterize clusters and check distribution and levels of genes of interest.
* Also in this case, the [__compiled version__](Compiled/2_Clusters.ipynb) of the notebook is available in the 'Compiled' folder

-----------

<br> </br>
