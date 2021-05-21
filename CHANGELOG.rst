=========
Changelog
=========


Version 5.1.1 (20210521)
========================
- Added a new entry to the namelist to tune the number of previous days to use to compute average temperature and suppress melt (originally 10 days in v5.1.0).
- Modified the source code accordingly, as well as made some housekeeping (e.g., added a2dTaC_MeanDaysSuppressMelt as new state variable and T_SuppressMeltDays as new output in place of T_10Days)
- Added processing time as new global attribute of output nc files. This attribute reports the timestamp when the output file was created. 

Version 5.1.0 (20210326)
========================
- First public release of S3M 5.1.0

