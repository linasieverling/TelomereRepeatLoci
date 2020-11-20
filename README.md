# TelomereRepeatLoci
## Snakemake workflow for detection of telomere repeat loci from WGS data

<img src="resources/images/telomere_repeat_locus_example.png" alt="Example of telomere repeat locus" width="500" />

This snakemake workflow detects telomere repeat loci within cancer genomes from WGS data. The input are BAM files from a tumor and a control sample (if available). In the first step, telomeric reads are extracted using the tool [TelomereHunter](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0). From the extracted telomeric reads, discordant reads are retrieved, where one mate is intratelomeric and the other mate is mapped to the chromosome. In regions with discordant reads, it then searches for clipped reads to find the precise position of the inserted telomere sequence.

<img src="resources/images/telomere_repeat_locus_schematic.png" alt="Detection of telomere repeat loci" width="700" />


If you are using the workflow, please cite:

**TelomereHunter – in silico estimation of telomere content and composition from cancer genomes**

Lars Feuerbach, Lina Sieverling, Katharina I. Deeg, Philip Ginsbach, Barbara Hutter, Ivo Buchhalter, Paul A. Northcott, Sadaf S. Mughal, Priya Chudasama, Hanno Glimm, Claudia Scholl, Peter Lichter, Stefan Fröhling, Stefan M. Pfister, David T. W. Jones, Karsten Rippe & Benedikt Brors

*BMC Bioinformaticsvolume 20, Article number: 272 (2019)*


**Alternative lengthening of telomeres in childhood neuroblastoma from genome to proteome**

Sabine A. Hartlieb, Lina Sieverling, Michal Nadler-Holly, Matthias Ziehm, Umut H. Toprak, Carl Herrmann, Naveed Ishaque, Konstantin Okonechnikov, Moritz Gartlgruber, Young-Gyu Park, Elisa Maria Wecht, Kai-Oliver Henrich, Larissa Savelyeva, Carolina Rosswog, Matthias Fischer, Barbara Hero, David T.W. Jones, Elke Pfaff, Olaf Witt, Stefan M. Pfister, Katharina Kiesel, Karsten Rippe, Sabine Taschner-Mandl, Peter Ambros , Benedikt Brors , Matthias Selbach, Lars Feuerbach, Frank Westermann

*under revision*

The workflow was also used in the following publication:

**Genomic footprints of activated telomere maintenance mechanisms in cancer**

Lina Sieverling, Chen Hong, Sandra D. Koser, Philip Ginsbach, Kortine Kleinheinz, Barbara Hutter, Delia M. Braun, Isidro Cortés-Ciriano, Ruibin Xi, Rolf Kabbe, Peter J. Park, Roland Eils, Matthias Schlesner, PCAWG-Structural Variation Working Group, Benedikt Brors, Karsten Rippe, David T. W. Jones, Lars Feuerbach & PCAWG Consortium

*Nature Communications volume 11, Article number: 733 (2020)*






### Detailed description of individual steps in the workflow

<img src="resources/images/TelomereRepeatLoci_workflow.png" alt="TelomereRepeatLoci snakemake workflow" width="700" />

#### 1. Run TelomereHunter

  Information on TelomereHunter can be found in the [publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0).
  

#### 2. Find candidate regions with discordant reads
#### 3. Find precise insertion sites with clipped reads
#### 4. Construct telomeric sequences at the telomere insertion sites
#### 5. Make IGV-like plot

