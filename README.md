# PACO (PAthway COmparator)
Shiny R web-app for comparing perturbed pathways associated to different phenotypes.

<b>URL:</b>
<a href="https://alpha.dmi.unict.it/shiny/users/gmicale/PathwayComparator/">https://alpha.dmi.unict.it/shiny/users/gmicale/PathwayComparator/</a>.

<b>Input data:</b>

Two types of text files can be uploaded by the user for pathway comparison:
- PHENSIM simulation file
- Custom file

Example files are available in the "Data" folder.

Currently supported organisms:
- Human
- Mouse

More organisms will be supported in the future.

<br/>

<b>PHENSIM simulation file:</b>

PHENSIM (PHENotype SIMulator) (<a href="https://phensim.tech/"></a>) is a tool developed to simulate the de-regulation of pathways biological elements, as a result of the over- or under-expression of user-specified molecules (e.g. genes or miRNAs).

De-regulation of pathway nodes is measured by an activity score, which can be positive or negative, denoting a biological element which is more or less active than normal condition.

After a new simulation has been launched and complete, PHENSIM simulation file can be downloaded as follows:
- Go to the "Simulations" panel on the left sidebar
- Click on "Show simulation" ('eye' icon on the right) to view more details about the simulation
- Go to the "Download results" box at the bottom of the page and click on "Download raw results"

Scores associated to nodes are the values "Activity Score" column.

Example:

<pre>
# Pathway Id	Pathway Name	Node Id	Node Name	Is Endpoint	Is Direct Target	Activity Score	P-Value	Adjusted P-Value	Log-Probabilities (Activation, Inhibition, Others)	Pathway Activity Score	Pathway p-value	Pathway Adjusted p-value	Pathway Log-Probabilities (Activation, Inhibition, Others)	Direct Targets	Average Node Perturbation	Average Pathway Perturbation
R-HSA-198753	ERK/MAPK targets	6197	RPS6KA3	Yes	No	0.0	0.9980000000000008	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9	0.0	0.5760000000000004	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9		0.0	-2.4006391023997053E-5
R-HSA-198753	ERK/MAPK targets	6196	RPS6KA2	No	No	0.0	0.9960000000000008	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9	0.0	0.5760000000000004	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9		0.0	-2.4006391023997053E-5
R-HSA-198753	ERK/MAPK targets	hsa-miR-199b-3p	hsa-miR-199b-3p	No	No	0.0	0.9970000000000008	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9	0.0	0.5760000000000004	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9		-6.314636151692764E-13	-2.4006391023997053E-5
</pre>

<br/>

<b>Custom file:</b>

The user can provide the app with any custom text file containing a score for each biological element.
The text file must contain three columns, "Organism","Node" and "Score, indicating for each biological entity (gene, miRNAs, compound) the organism to which it belongs or refers and a perturbation score.
Organism column is the common name of one of the supported organisms (e.g. 'Human', 'Mouse').
Node column specifies the id of the biological element, which must be: 
- The Entrez ID for genes
- The miRBase entry name for miRNAs
- The ChEBI id, the KEGG COMPOUND id or the KEGG GLYCAN id for chemical compounds
- The KEGG DRUG id for drugs
Perturbation score can be any real number (positive, negative or zero).

Example:

<pre>
Organism	Node	Score
Human	6197	-3.74
Human	6196	0
Human	hsa-miR-199b-3p	0
Human	hsa-miR-128-3p	2.5
Human	hsa-miR-214-3p	0
Human	chebi:43474	0
</pre>

<br/>

<b>References:</b>

- Alaimo S, Rapicavoli RV, Marceca GP, La Ferlita A, Serebrennikova OB, et al. (2021). <i>PHENSIM: Phenotype Simulator.</i> PLOS Computational Biology 17(6): e1009069. <a href="https://doi.org/10.1371/journal.pcbi.1009069">https://doi.org/10.1371/journal.pcbi.1009069</a>

