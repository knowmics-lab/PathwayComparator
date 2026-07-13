# PACO (PAthway COmparator)
Shiny R web-app for comparing perturbed pathways associated to different phenotypes.

<b>URL:</b>
<a href="https://paco.dioncogen.eu/">https://paco.dioncogen.eu/</a>.

<h2>Input data</h2>

Three types of text files can be uploaded by the user for pathway comparison:
- Custom file
- MITHrIL perturbation file
- PHENSIM simulation file

Example files are available in the "Data" folder.

Currently supported organisms:
- Human
- Mouse
- Rat
- Worm
- Fly
- Zebrafish
- Arabidopsis

More organisms will be supported in the future.

<h3>Custom file</h3>

The user can provide the app with any custom text file containing a score for each biological element.
The first line of the text file must contain the common name of one of the supported organisms (e.g. 'Human', 'Mouse') to which biological elements belong. 
The following lines indicate for each biological entity (gene, miRNAs, compound) a perturbation score.
Nodes are referred to by the id of the corresponding biological element, which must be: 
- The Entrez ID for genes
- The miRBase entry name for miRNAs
- The ChEBI id, the KEGG COMPOUND id or the KEGG GLYCAN id for chemical compounds
- The KEGG DRUG id for drugs
Perturbation score can be any real number (positive, negative or zero).

Example:

<pre>
Human
6197	-3.74
6196	0
hsa-miR-199b-3p	0
hsa-miR-128-3p	2.5
hsa-miR-214-3p	0
chebi:43474	0
</pre>

<h3>MITHrIL perturbation file</h3>

MITHrIL evaluate the de-regulation of pathways due to expression changes of one or more pathway elements. MITHrIL requires a list of biological entities with their log-fold expression changes as input. De-regulation of pathway nodes is measured by a perturbation score, which can be positive or negative. MITHril is available as a command line tool (<a href="https://github.com/alaimos/mithril-standalone">https://github.com/alaimos/mithril-standalone</a>), and its output is a perturbation file that can be directly uploaded as input to PACO. The scores associated with nodes and used by PACO are the values of the "Perturbation" column.

Example:

<pre>
# Pathway Id	Pathway Name	Gene Id	Gene Name	Perturbation	Accumulator	pValue
path:hsa00190	Oxidative phosphorylation - Enriched	64077	LHPP	0.45	1.36	1.0
path:hsa00190	Oxidative phosphorylation - Enriched	5464	PPA1	2.10	3.56	1.0
path:hsa00190	Oxidative phosphorylation - Enriched	hsa-miR-101-3p	hsa-miR-101-3p	-0.79	-2.45	1.0
</pre>

<h3>PHENSIM simulation file</h3>

PHENSIM (PHENotype SIMulator) (<a href="https://phensim.tech/">https://phensim.tech/</a>) is a tool developed to simulate the de-regulation of pathways biological elements, as a result of the over- or under-expression of user-specified molecules (e.g. genes or miRNAs). De-regulation of pathway nodes is measured by an activity score, which can be positive or negative, denoting a biological element which is more or less active than normal condition. Scores associated to nodes are the values of "Activity Score" column.

Example:

<pre>
# Pathway Id	Pathway Name	Node Id	Node Name	Is Endpoint	Is Direct Target	Activity Score	P-Value	Adjusted P-Value	Log-Probabilities (Activation, Inhibition, Others)	Pathway Activity Score	Pathway p-value	Pathway Adjusted p-value	Pathway Log-Probabilities (Activation, Inhibition, Others)	Direct Targets	Average Node Perturbation	Average Pathway Perturbation
R-HSA-198753	ERK/MAPK targets	6197	RPS6KA3	Yes	No	0.0	0.9980000000000008	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9	0.0	0.5760000000000004	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9		0.0	-2.4006391023997053E-5
R-HSA-198753	ERK/MAPK targets	6196	RPS6KA2	No	No	0.0	0.9960000000000008	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9	0.0	0.5760000000000004	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9		0.0	-2.4006391023997053E-5
R-HSA-198753	ERK/MAPK targets	hsa-miR-199b-3p	hsa-miR-199b-3p	No	No	0.0	0.9970000000000008	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9	0.0	0.5760000000000004	1.0	-20.72326583994641,-20.72326583994641,-1.999999945436137E-9		-6.314636151692764E-13	-2.4006391023997053E-5
</pre>

<h2>Examples of usage</h2>

We want to compare perturbed immune system pathways (specifically the "Interferon Signaling" pathway) in mice and humans, following the upregulation of Interferon (IFN)-stimulated gene 15 (ISG15), a ubiquitin-like protein that functions both as an extracellular cytokine and an intracellular protein modifier.

As a preliminary step, we first run two simulations using PHENSIM by upregulating ISG15 in humans and mice, respectively.

Connect to the <a href="https://phensim.tech/">PHENSIM web portal</a>[^1]. Click on "Simulations" on the left sidebar. From the simulation page that appears, clik on "New simple simulation". 

![PHENSIM_Portal](https://github.com/user-attachments/assets/0eca77c6-0251-44cd-bee1-700c17196415)

Following the guided procedure, indicate a name for the new simulation and select "Homo sapiens" as organism. Next, write "ISG15" in the filter of "NAME" column and click on the red up-arrow to include ISG15 as up-regulated gene in the simulation. Finally, select the option "Add REACTOME pathways" in the "SELECT OPTIONAL PARAMETERS" tab to include REACTOME pathways in the simulation and click on "Create simulation" to run PHENSIM.

![PHENSIM_Simulation](https://github.com/user-attachments/assets/7a7bc7ee-3cf2-4209-b0b0-42f87fcd5041)

The launched simulation will then appear on the list of all simulations launched so far by the user using the PHENSIM portal, accessible by clicking on "Simulations" on the left sidebar. When the simulation ends, i.e. the value of the "STATUS" column in the table is "Completed", download the PHENSIM simulation file as follows:
- Go to the "Simulations" panel
- Click on "Show simulation" ('eye' icon on the right)
- Go to the "Download results" box at the bottom of the page and click "Download raw results."

![PHENSIM_Results](https://github.com/user-attachments/assets/f52c9d9b-45ce-40bf-99bc-394b99759c53)

To perform a PHENSIM simulation by upregulating ISG15 in mice, repeat the same steps on the PHENSIM portal, selecting now "Mus Musculus" as organism.

Next, upload the two files into PACO and click on the "Compare" button. 

![PACO_Upload](https://github.com/user-attachments/assets/9a1ad865-2e21-433d-948e-904986d64cde)

<h3>Search by pathway</h3>

In the first case study, we want to compare the interferon signaling pathway in Human and in Mouse and, in particular, the role of human ISG15 gene and its ortholog in mouse in this pathway. 

First, in the "2. Search by" Section of the visualization panel, select "Pathway". Then, in the Section "3. Compare & filter", choose "Interferon Signaling" as "Pathway" and "ISG15" as "Gene".

![PACO_ByPathways](https://github.com/user-attachments/assets/3bb47fe2-7d5f-4ff7-87b8-8c718ee0a21b)

A multilayer network will be visualized, where layers correspond to the perturbation score files we are comparing. Each layer shows the ISG15 gene, its direct neighbors and the links between all these nodes in the selected "Interferon Signaling" pathway of the organism to which the score file refers. Dashed lines connect homologous genes, equal chemical compounds, equal drugs or equal miRNAs in the compared layers. Perturbation scores are represented by using colors (red for positive scores, blue for negative scores) and square nodes represent pathway endpoints (nodes with no incoming edges in the pathway).

An alternative visualization consists in showing all the paths in the interferon signaling pathway going backward from ISG15 gene to its upstream predecessors up to a specified distance.

For example, we want to visualize all backward paths from ISG15 to its predecessors up to distance 3. To do this, select the visualization option 
"Show all paths to selected nodes" and set "Max hops" to 3.

![PACO_ByPathwaysHops](https://github.com/user-attachments/assets/3f26b92f-55e7-45f7-a1b9-f3acf26d01b9)

In the visualized network of all possible paths, the size of the nodes is proportional to their distance to ISG15: the higher is the distance, the smaller is the node.

<h3>Search by node</h3>

As a second case study, let's investigate all pathways in which ISG15 gene is involved and how ISG15 gene is connected to the other nodes in these pathways.

First, in the "2. Search by" Section of the visualization panel, select "Node". Then, in the Section "3. Compare & filter", choose "ISG15" as "Node". In the "Pathway" dropdown-menu click on "Select All" to select all pathways in which ISG15 gene is present in all compared layers. Finally, select (if not already done) "Show selected nodes in selected pathways". The result is the same type of multilayer network descripted in the first case study, depicting ISG15 gene, its direct neighbors and the links between all these genes in the network given by the union of all selected pathways.

![PACO_ByNodes](https://github.com/user-attachments/assets/1d21af5c-9dfb-4821-a3af-6343373b8228)

Likewise, it is possible to visualize all backward paths from ISG15 to upstream nodes up to a certain distance in the network given by the union of all selected pathways, by clicking on "Show all paths to selected nodes".

<h2>References:</h2>

- Micale G, Alaimo S, Pulvirenti A (2025). <i>PACO: a Shiny app for comparing perturbed pathways associated with different phenotypes.</i> Bioinformatics Advances 5(1). <a href="https://doi.org/10.1093/bioadv/vbaf212">https://doi.org/10.1093/bioadv/vbaf212</a>

- Alaimo S, Rapicavoli RV, Marceca GP, La Ferlita A, Serebrennikova OB, et al. (2021). <i>PHENSIM: Phenotype Simulator.</i> PLOS Computational Biology 17(6): e1009069. <a href="https://doi.org/10.1371/journal.pcbi.1009069">https://doi.org/10.1371/journal.pcbi.1009069</a>

- Alaimo S, Giugno R, Acunzo M, Veneziano D, Ferro A, Pulvirenti A (2016). <i>Post-transcriptional knowledge in pathway analysis increases the accuracy of phenotypes classification.</i> Oncotarget 7(34):54572-54582. <a href="https://doi.org/10.18632/oncotarget.9788">https://doi.org/10.18632/oncotarget.9788</a>


[^1]: Before using PHENSIM web portal, a registration is required. After registration, log in to start new simulations.
