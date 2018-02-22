# PEG4D_data_cleaner
This package is intended to the PEG4D users to easily clean their datasets from unwanted and unnecessary peaks in their data. In particular, it focuses on the elimination of false positives peaks coming from different source, like matrix effects, laboratory contaminants, column bleedings, or wrongly assigned double peaks.

### how to install the package
This package has 3 other dependencies:
"stringr", "chron" and "lubridate".
You also need devtools to install this package
First install such packages
*install.packages(c("chron", "lubridate", "stringr", "devtools"))*
Then call devtools package and install the GCxGC_data_cleaner package
*library(devtools)
*
