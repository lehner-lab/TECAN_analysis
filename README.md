![yeast](./yeast.png)

TECAN-reader based yeast growth measurements in R
=================================================

Maxi's [script](.TECAN_384_well.R) and some [toy data](.TECAN_growths.xlsx) for calculating maximum yeast growth rates from 384-well TECAN measurements, in R.

```
R --version
```

The current R version is also shown when opening RStudio or the R Console.

Script requires the following R packages: [**`readxl`**](https://cran.r-project.org/web/packages/readxl/index.html), [**`growthrates`**](https://cran.r-project.org/web/packages/growthrates/index.html), [**`beeswarm`**](https://cran.r-project.org/web/packages/beeswarm/index.html), [**`scales`**](https://cran.r-project.org/web/packages/scales/index.html).

Although care has been taken to make the code distribution-independent, it is possible that some of the scripts only work on Unix/MacOS systems, and may need to be modified in order to run on Windows systems.
