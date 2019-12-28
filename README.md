# Enhancer-detection

This file provides basic information to run "Enhancer_detection.pl" software, including input/output files.

Input files:
------------
1) CE7002_Dmel.fa (e.g Eve_stripe_3/7 enhancer region from D. melanogaster)
2) CE7002_Dpse.fa (e.g Orthologous intergenic region from D. pseudoobscura)
3) Dmel_background (Background frequency score from D. melanogaster)
4) Dpse_background (D. pseudoobscura background frequency score)
5) Dpse_null_distribution (Background score distribution from D. pseudoobscura) 

Output files:
-------------

1) plots/CE7002_Dpse_max.xls - the profile contains the mixed metric score of each sliding windows.

2) graphs/CE7002_Dpse_500.pdf - the scanning of Eve_stripe_3/7 enhancer against the orthologous intergenic region from D. pseudoobscura.

3) CE7002-Dpse.res - output of the scanning process
		     "score" - global maximum mixed metrix score
		     "max" - window number with the global maximum score
		     "ends" - the set of consecutively merged windows that exceeds the threshold value

Running Instructions:
--------------------

Make sure that the input files are exists in the current directory

You need to provide your current directory PATH name in the second line, (e.g) "usr lib /HOME DIRECTORY/lib/perl5/site_perl/5.8.8/" into "usr lib /Home/Jone/lib/perl5/site_perl/5.8.8/"

run command: perl Enhancer_detection.pl Enhancer_id(e.g CE7002) Species(e.g Dpse)

The output file will be created in the current directory named "CE7002-Dpse.res" and the profile output file will be directory to "plots" directory (will be created by user). Using R plot the pdf file will be created in the graphs directory.
