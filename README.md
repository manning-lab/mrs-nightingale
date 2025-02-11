This pipeline is hardcoded for particular data, but may hopefully be useful as a reference.

`mrs.R` does the following:
* Harmonizes metabolite names between the two input files
* Creates histograms of each metabolite to see if the distributions are normal
* Inverse-normal transforms each metabolite (because the distributions were seen to be not normal)
* Collects multiple vectors of metabolite weights from the 2nd input file. (There is a different vector of weights to be used depending on a sample's age group / sex group / T2D endpoint type / analysis method. Scores are calculated for every possible combination of characteristics, for every sample -- the correct score for each sample can be chosen later. The score calculation is computationally fast enough that code simplicity was chosen over efficiency.)
* Each score vector is multiplied against the metabolite matrix. This results in an output file with columns of scores.
