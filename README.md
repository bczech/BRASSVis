About
=====
[Arriba](https://github.com/suhrig/arriba) software provides a script for visualization of gene fusions detected in RNA-seq data. It accurately visualizes all detected fusions by presenting their breakpoints, coverage, cytobands, protein domains, and additional alignment informations. However, this tool is not suitable for WGS data.

Hence we created a BRASSViz that is a wrapper of modified Arriba visualization tool. Using this tool one is able to visualize BEDPE files, generated by BRASS structural variant caller.

Repository structure
=====
    .
    ├── annotations             # Annotation files (for human genome GRCh37 and GRCh38)
    ├── ref                     # Reference files with cytobands and protein domains
    ├── src                     # Source files
    ├── tst                     # Example files (BEDPE file used in the example and reports)
    ├── LICENSE
    └── README.md
    

Installation
=====
First of all, clone this repo by: `git clone https://github.com/bczech/BRASSViz`

In order to use a BRASSViz, installatin of R 4.0 is required. Check [CRAN](https://cran.r-project.org/) website to download the latest R version.

BRASSViz requires installing additional R packages that can be downloaded by running R script stored in [src/install_packages.R](src/install_packages.R).


Usage
=====

BRASSViz provides multiple options to visualize your BEDPE file:

```
usage: BRASSViz.R [options]

options:
--inputFile=<file> BRASS file path.

--exonsFile=<file> annnotation file [default: ./annotation/Homo_sapiens.GRCh38.101.chr.gtf].

--outputFile=<file> output file name [default: ./output.pdf].

--confidence=<character> confidence of fusions [default: high].

--fusionFlagThreshold=<numeric> fusion flag threshold [default: 720].

--transcriptSelection=<character> transcript selection [default: provided].

--pdfWidth=<numeric> pdf width [default: 11.692].

--pdfHeight=<numeric> pdf height [default: 8.267].

--fusionsFile=<file> fusions file title [default: Title].

--cytobandsFile=<file> cytobans file.

--proteinDomainsFile=<character> protein domains file.

--squishIntrons=<logical> squish introns [default: TRUE].

--render3dEffect=<logical> render 3D effects [default: FALSE].

--alignmentsFile=<file> alignment file.

--color1=<character> color 1 used in the visualization [default: #e5a5a5].

--color2=<character> color 2 used in the visualization [default: #a7c4e5].

--showVicinity=<numeric> show vicinity [default: 0].

--printExonLabels=<logical> print exon labels [default: TRUE].

--minConfidenceForCircosPlot=<character> minimal confidence for circos plot [default: medium].

--mergeDomainsOverlappingBy=<numeric> merge domains overlapping by [default: 0.9].

--optimizeDomainColors=<logical> optimize domain colors [default: FALSE]

--fontSize=<numeric> font size [default: 1]
```
#### Generating visualizations

##### Generating simple visualization
```
src/BRASSViz.R -i tst/T1_runs1-2_60X_vs_N1_SRR7890942_SRR7890943_30X.brass.annot.bedpe \\
-o tst/report.pdf \\
-e annotation/Homo_sapiens.GRCh37.75.gtf.gz \\
--proteinDomainsFile ref/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3 \\
--cytobandsFile ref/cytobands_hg19_hs37d5_GRCh37_v2.1.0.tsv
```
This command generates a visualization. `-i` indicates the input bedpe file, `-o` path of output file (visualization), `-e` path to genome annotation in GTF format, `--proteinDomainsFile` path to gff3 file with protein domains (included in [ref][ref/] directory), `--cytobandsFile` indicates a path to cytobands (included in [ref][ref/] directory as well).
Note that some fusions are excluded from the analysis based on the fusion flag threshold criterion (720 set by default).

Exemplary output can be seen [here.](tst/report.pdf)

Authors
=====
* **Bartosz Czech**
* **Paweł Sztromwasser**
* **Marzena Wojtaszewska**

