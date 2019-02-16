# essentials

ESSENTIALS: Software for rapid analysis of transposon insertion sequencing data
Aldert Zomer1,2, Peter Burghout1, Hester Bootsma1, Peter WM Hermans1 and Sacha van Hijum2,3,4,5

1 Laboratory of Pediatric Infectious Diseases, Radboud University Nijmegen Medical Centre, P.O. Box 9101, 6500 HB, Nijmegen, The Netherlands.
2 Centre for Molecular and Biomolecular Informatics, Nijmegen Centre for Molecular Life Sciences, Radboud University Nijmegen Medical Centre, P.O. Box 9101, 6500 HB, Nijmegen, The Netherlands.
3 NIZO food research, Kluyver Centre for Genomics of Industrial Fermentation, P.O. Box 20, 6710 BA Ede, The Netherlands
4 TI Food and Nutrition, P.O. Box 557, 6700 AN Wageningen, The Netherlands
5 Netherlands Bioinformatics Centre, 260 NBIC, P.O. Box 9101, 6500 HB Nijmegen, the Netherlands

High-throughput analysis of genome-wide transposon mutant libraries is a powerful tool for (conditional) essential gene discovery. Recently, several next generation sequencing approaches, e.g. Tn-seq, INseq and TraDIS, have been developed that accurately map the site of transposon insertions by mutant-specific amplification and sequence readout of DNA flanking the transposon insertions site, assigning a measure of essentiality based on the number of reads per gene or per mutant. However, analysis of these large and complex datasets is hampered by the lack of an easy to use and automated tool for transposon insertion sequencing data.

To fill this gap, we developed ESSENTIALS, an open source, web-based software tool for researchers in the genomics field utilizing transposon insertion sequencing analysis. It accurately predicts (conditionally) essential genes and offers the flexibility of using different sample normalization methods, genomic location bias correction, data preprocessing steps, appropriate statistical tests and various visualizations to examine the results, while requiring only a minimum of input and hands-on work from the researcher.

We successfully applied ESSENTIALS to in-house and published Tn-seq, TraDIS and HITS datasets and we show that the various pre- and post-processing steps on the sequence reads and count data with ESSENTIALS markedly improve the sensitivity and specificity of predicted gene essentiality.

Essentials at CMBI

A working version of the ESSENTIALS pipeline can be accessed at http://bamics2.cmbi.ru.nl/websoftware/essentials/essentials_start.php
