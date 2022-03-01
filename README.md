# Machine Learning for Drug Discovery

This repository contains all code and examples used in the 
ACS *in focus* book [Machine Learning for Drug Discovery](https://pubs.acs.org/doi/book/10.1021/acsinfocus.7e5017).

# Using this repo

All code was prepared to be executed in a Python environment.
For the reader's convenience, we provide a [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment file
that can be used to create a python environment with previously 
selected and tested packages.

To install Conda in your system, please refer to [Conda's installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

After Conda ins installed and ready, navigate to the `conda_env` folder and execute the following command to create the environment:

`conda env create -f environment_ML_Env.yml`

Now activate your environment with:

`conda activate environment_ML_Env`

And you will be ready to initialize the Jupyter notebook form the parent folder using:

`jupyter notebook`
