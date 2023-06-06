![yeast](./yeast.png)

TECAN-reader based yeast growth measurements in R
=================================================

# Installation

We recommend using [this yaml file](tecan.yaml) to create a dedicated Conda environment with all necessary dependencies (as explained below).

1. Install the [Conda](https://docs.conda.io/) package/environment management system (if you already have Conda skip to step 2):

   On MacOS, run:
   ```
   $ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
   $ sh Miniconda3-latest-MacOSX-x86_64.sh
   ```
   On Linux, run:
   ```
   $ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   $ sh Miniconda3-latest-Linux-x86_64.sh
   ```

   **IMPORTANT:** If in doubt, respond with "yes" to the following question during installation: "Do you wish the installer to initialize Miniconda3 by running conda init?". In this case Conda will modify your shell scripts (*~/.bashrc* or *~/.bash_profile*) to initialize Miniconda3 on startup. Ensure that any future modifications to your *$PATH* variable in your shell scripts occur **before** this code to initialize Miniconda3.

2. Clone the TECAN_analysis GitHub repository:
   ```
   $ git clone https://github.com/lehner-lab/TECAN_analysis.git
   ```

3. Create the tecan Conda environment and activate it:
   ```
   $ conda env create -f TECAN_analysis/tecan.yaml
   $ conda activate tecan
   ```

4. Start an R session and install the growthrates package:
   ```
   > install.packages("growthrates")
   ```
   
# Usage
   ```
   $ conda activate tecan
   $ TECAN_analysis/TECAN_384_well.R -h
   ```
   
# Demo
   ```
   $ conda activate tecan
   $ TECAN_analysis/TECAN_384_well.R TECAN_analysis/TECAN_growths.xlsx
   ```

# Dependencies

Requires the following R packages (see [tecan.yaml](tecan.yaml)):
* readxl
* growthrates
* argparser
* beeswarm
* scales
* data.table
  
Although care has been taken to make the code distribution-independent, it is possible that some of the scripts only work on Unix/MacOS systems, and may need to be modified in order to run on Windows systems.
