<img src="https://github.com/lehner-lab/TECAN_analysis/blob/main/yeast.png" width="500">


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
   $ TECAN_analysis/tecan.R -h
   usage: tecan.R [--] [--help] [--opts OPTS] [--method METHOD]
          [--parameter PARAMETER] [--wells WELLS] [--deadThreshold
          DEADTHRESHOLD] [--lagThreshold LAGTHRESHOLD] [--outputPrefix
          OUTPUTPREFIX] [--designPath DESIGNPATH] [--plotWidth PLOTWIDTH]
          [--plotHeight PLOTHEIGHT] [--ODThreshold ODTHRESHOLD] excel_path
   
   TECAN-reader based yeast growth measurements and summary plots
   
   positional arguments:
     excel_path           Path to the Excel file
   
   flags:
     -h, --help           show this help message and exit
   
   optional arguments:
     -x, --opts           RDS file containing argument values
     -m, --method         Maximum growth rate method [default: heuristic]
     -p, --parameter      h-parameter; number of consecutive time points
                          to evaluate maximum growth rate [default: 15]
     -w, --wells          Comma-separated list of well ids [default: all]
     -d, --deadThreshold  Growth rate threshold for dead variants
                          [default: 0.05]
     -l, --lagThreshold   Lag time threshold (in hours) for problematic
                          variants [default: 48]
     -o, --outputPrefix   Output path prefix [default: no output file;
                          print results to stdout]
     --designPath         Path to the plain text file with 'Well',
                          'Plasmid' and 'Well class' columns [optional].
                          'Well' column is required whereas 'Plasmid' and
                          'Well class' columns are optional; valid 'Well
                          class' column values are: 'Sample', 'Positive
                          control', 'Negative control', 'Blank' and
                          'Discarded'; plasmid-specific growth curves are
                          produced if 'Plasmid' column is supplied;
                          plasmid-specific boxplots are produced if
                          'Plasmid' and 'Well class' columns are supplied.
     --plotWidth          Plot width in inches [default: 8]
     --plotHeight         Plot height in inches [default: 8]
     -O, --ODThreshold    Minimum optical density required to escape
                          'deadThreshold' [default: 0]
   ```
   
# Demo
   ```
   $ conda activate tecan
   $ TECAN_analysis/tecan.R TECAN_analysis/TECAN_growths.xlsx
   ```

# Dependencies

Requires the following R packages (see [tecan.yaml](tecan.yaml)):
* readxl
* growthrates
* argparser
* beeswarm
* scales
* data.table
* ggplot2
* reshape2
  
Although care has been taken to make the code distribution-independent, it is possible that some of the scripts only work on Unix/MacOS systems, and may need to be modified in order to run on Windows systems.
