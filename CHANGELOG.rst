=========
Changelog
=========

Version 5.3.4 (20250206)
========================
- modified configure.sh to support deployments in a docker

Version 5.3.3 (20241113)
========================
- correct bug in restart procedure wrt loading ice thickness

Version 5.3.2 (20241031)
========================
- Performed some minor changes to S3M_Module_Tools_Generic.f90, function filter_array_unique, to neglect -9999.0 while counting glacier IDs (this significantly speeds up that routine)

Version 5.3.1 (20240917)
========================
- Changed license to EUPL. 

Version 5.3.0 (20231215)
========================
- Improved handling of input, assimilation, and restard data to allow multiple, parallel runs on the same input data. The improved approach generates temporary copies of inpit files using process PID to avoid problems with concurrent reading of the same netCDF file.

Version 5.2.0 (20230208)
========================
- Improved handling of restard data, with a more rigorous check on the existence of both file and NC layers to avoid spurious results
- Extended number of output layers in the basic output mode to also include cumulative daily melt and snowfall at all hours, and LWC
- Modified reset of cumulative daily melt and snowfall in S3M_Module_Phys_Snow

Version 5.1.2 (20210720)
========================
- Added a provision in S3M_Main.f90 to shift the frequency at which output file are dumped by iDtData_Forcing compared to iDtData_Output. In the typical scenario where iDtData_Output is 3600 s (1 h), this convenient feature makes sure that output files are saved 1 h before iDtData_Output (this guarantees that the useful 11PM file are saved and can be used as a restart).

Version 5.1.1 (20210521)
========================
- Added a new entry to the namelist to tune the number of previous days to use to compute average temperature and suppress melt (originally 10 days in v5.1.0).
- Modified the source code accordingly, as well as made some housekeeping (e.g., added a2dTaC_MeanDaysSuppressMelt as new state variable and T_SuppressMeltDays as new output in place of T_10Days)
- Added processing time as new global attribute of output nc files. This attribute reports the timestamp when the output file was created. 

Version 5.1.0 (20210326)
========================
- First public release of S3M 5.1.0

