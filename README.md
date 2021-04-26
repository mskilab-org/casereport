# casereport
Genome graph based case reports

# For developers
## Manifest of files
| File      | Description |
| ----------- | ----------- |
| wgs.report.rmd | the rmd report to be knitted |
| wgs.run.R | the main R script to prepare data and plots |
| wgs.casereport.task | new task file |
| config.R | the R script of preset colors and other static global variables |
| utils.R | the R script of utility functions |
| flow.deploy | deploy file |
| run.sh | bash script to start running the module |


## What to do
In the rmd file organize sections of contents, fill in with text and
figures. Figures should be pre-made within the main R script as PNG
files and inserted using knitr::include_graphics().

