
-- codes to characterize the relationship between gene transcription divergence and PEI number across species

1. System requirements

(1) centOS 7
(2) python 2.7
(3) R 3.5
(4) R pacakage MatchIt 3.0.2

2. Installation guide

No installation needed
Directly run the following codes one after another:

(1) 01.prepare_formatch.py
(2) 02.formatch.R
(3) 03.prepare_forplot.new.py
(4) 04.plot_scatter.R

3. Demo

Input:
divergence_time.txt
merge_codingtable_forcluster.final.tpm.nochicken.nomouse.normed.filtered.xls
merge_Enum.nochicken.nomouse.filtered.xls
orth_genes.nochicken.nomouse.txt
sample_infor.rest.txt
sus_scrofa.Tau.7tissues.xls

Output:
moreEnumandControl_cor.correct.filtered.pdf

4. Instructions for use

Change the hard-coded paths in the code to your own paths where you place the input and output files
