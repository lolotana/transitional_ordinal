# Autoregressive_transitional_ordinal
R simulation code to replicate simulation study

This is a simplified & commented R code used for the simulation part of the manuscript
"Autoregressive transitional ordinal model to test for treatment effect in neurological trials with complex endpoints".
Please refer to the manuscript for details on the simulation study.

The R code only runs if provided with patient data formatted accordingly to the given instructions.
The patient data is not publicly available. Interested researcher may apply for data access to emsci.org, which is 
usually granted for research-only purposes.

The R code is set up to run on a simplified version of the whole simulation. Change settings (original settings commented out) under the "set simulation settings" code chunck to run full simulation. The latter is very (!) time-consuming, you may want to run it on a multi-core server and/or parallelize it.

