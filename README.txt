The detail of the implementation is written in the report.
Please refer to "report.pdf" for more details.

In order to integrate the code into PBRT system (only version 2 is supported),
please put "medianCutEnvironmentLight.cpp" and "medianCutEnvironmentLight.cpp.h" under <pbrt_folder>/src/shapes/.
(the 2 other files are for comparison only)

Then change the original <pbrt_folder>/src/core/api.cpp to
"api.cpp" in this repository.
