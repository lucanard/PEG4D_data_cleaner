# PEG4D_data_cleaner
This package is intended for the PEG4D users to easily clean their datasets from unwanted and unnecessary peaks in their data. In particular, it focuses on the elimination of false positives peaks coming from different source, like matrix effects, laboratory contaminants, column bleedings, or wrongly assigned double peaks.

## How the package works:
The input file for the following package is the data.matrix obtained from the chromatof software after peak detection and alignment. The matrix needs to be saved as csv file and imported as:
*read.csv("file.directory", header = T, stringsAsFactors = F)*
In the first step, the package eliminates all the redundant columns of the files, keeping only the Peak, first retention time, second retention time, the quantification mass, and the peak areas per each sample across all the compounds.
The next step is to eliminate all the variables present in the blank samples, below the (selectable) limit of detection.
It also eliminates all the variables known as column bleedings.
In the further step, all the variables having a similar retention time and identified under the same name are grouped, within a time span limit (selectable).
In the last step, it eliminates redundant peaks having a maximum value (maximum among all the samples) below a threshold limit (selectable). We suggest to use 50000, because injecting several standards we observed lack of linearity below this limit for all of them.

### how to install the package
This package has 3 other dependencies:
"stringr", "chron" and "lubridate".
You also need devtools to install this package
First install such packages
*install.packages(c("chron", "lubridate", "stringr", "devtools"))*
Then call devtools package and install the GCxGC_data_cleaner package
*library(devtools)
install_github("GCxGC.Leco.analyzer", repo = "lucanard/PEG4D_data_cleaner")
library(GCxGC.Leco.analyzer)*
