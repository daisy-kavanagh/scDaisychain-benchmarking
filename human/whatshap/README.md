Readbased phasing was performed with whatshap_run.sh.
The script compare_phase_concordance.py was then used to compare concordance of within phase set SNP assignment from whatshap and the scDaisychain results.
In order to further assess the impact of discordance between the 2 methods, the scDaisychain produced vcfs were corrected such that:
1) The majority orientation between scDaisychain and whatshap was selected by choosing the H1/H2 orientation that maximimises concordance.
2) Any discordant SNPs that were discordant to the majority orientation in the scDaisychain vcf were flipped to match their whatshap assignment.

The levels of concordance per phase set and the impact on Xi expression was then assessed using analyse_whatshap_comparison.ipynb

