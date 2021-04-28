
-- codes to characterize gene coexpression associated with TAD

1. System requirements

(1) PBS (Protable Batch System)
(2) centOS 7
(3) python 2.7
(4) R 3.4

2. Installation guide

No installation needed
Directly run the following codes one after another:

(1) 01.convert_mRNA_forcor.py
(2) 02.getall_cors.py
(3) 03.get_TSSfromgtf.py
(4) 04.filter_hasexprs.py
(5) 05.get_sepnumber_group
(6) 06.prepare_forfinalplot.py, 06.prepareall_forfinalplot.py

3. Demo

Input:
20700_protein_coding.exp.gene.name.tpm.rmPAD.0.1.txt
sus_scrofa.analysis.gtf.zip
total.combined.domaincalls
TSS_PG.txt

Output:
forplot_PG_BC1.txt

4. Instructions for use

Change the hard-coded paths in the code to your own paths where you place the input and output files
