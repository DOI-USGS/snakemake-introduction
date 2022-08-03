# Issue 0: Installing Snakemake
In order to walk through this learning module, we will create a new Python environment and install Snakemake and other required packages.

Please make sure you have Anaconda installed on your computer before proceeding with the following steps.

If you are working in Windows, you will need to open your Anaconda prompt to follow the steps below. If you are working in Linux or macOS, you can just open your command line.

## Install using environment.yaml
You can create a Conda environment for this tutorial with all the required packages in one step by running the following command:
`
conda env create -f environment.yaml
`

Note: You may have issues installing with the environment.yaml file. This could be due to the fact that Conda sometimes has trouble installing Snakemake. If you received any errors when creating the environment, skip down to the **Create environment and install packages manually** section.

If you did not receive any errors when creating the environment, test that the Conda environment was created successfully by activating it: `conda activate snakemake-tutorial`

If that worked, let's also test out your installation of Snakemake by running `snakemake --help`. You should see the documentation for Snakemake printed to your console.

If you are unable to activate your Conda environment, you will need to start over with the manual installation. Proceed to the **Create environment and install packages manually** section.

If you activated your Conda environment, but Snakemake failed to run successfully, delete any failed installations of Snakemake by running `conda remove snakemake`, and then pick up with **Step 3** in the **Create environment and install packages manually** section.

## Create environment and install packages manually
You can manually create your Conda environment and install the required packages by following the steps below.

1. Create a new Conda environment: `conda create -n snakemake-tutorial python=3.8`

2. Activate the environment: `conda activate snakemake-tutorial`

3. Install Mamba, a fast and robust replacement for the Conda package manager that is better able to handle the installation of Snakemake: `conda install -c conda-forge mamba`. *Note: We are only installing Mamba in our snakemake-tutorial Conda environment, but it can be used to install other python packages as well. If you would like to use Mamba to install other python packages, consider installing it to your base Conda environment.*

4. Install Snakemake (please [read the docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for more details on installation). If you are working in a Mac or Linux OS, you can install the full version of Snakemake: `mamba install -c conda-forge -c bioconda snakemake`. If you are working in Windows, some of the dependencies in the full version of Snakemake will not work for you. You can install the minimal version of Snakemake:`mamba install -c conda-forge bioconda::snakemake-minimal`

5. Test out your installation of Snakemake by running `snakemake --help`. You should see the documentation for snakemake printed to your console.

6. Install any additional packages needed for this tutorial (listed in the environment.yaml file) using Conda or Mamba.




**Last Resort**: If neither of the above options work, you can try installing Snakemake or snakemake-minimal with the appropriate Conda command:

`
conda install -c conda-forge bioconda::snakemake
conda install -c conda-forge bioconda::snakemake-minimal
`
