
# This file is ignored when building the package. It contains scripts for development purposes.
#
# Use this script to check that the package can be installed without errors on Windows, Linux and Mac.
if ("devtools" %in% installed.packages()) install.packages("devtools")
devtools::check_rhub(
  platforms = c("windows-x86_64-devel", "windows-x86_64-oldrel", "windows-x86_64-patched", "windows-x86_64-release")
) # Currently (2020-02-17) it only checks on Windows due to system requirement (rhub issue, not ours).
