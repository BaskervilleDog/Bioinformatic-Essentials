## Bioinformatics and Predictions

Use of experimentation its too expensive and take too long to characterize organisms. Therefore, using bioinformatics software's and algorithms helps speed up the process. The statistical analysis of rates of homologous recombination of different genes could determine their order on certain chromosome, helping in creating rough localization of know genes relative to each other.

#### **Ab Initio Prediction**

Ab initio gene prediction is a computational approach that predicts gene structures based solely on the DNA sequence itself, without relying on existing annotations. This method uses models of gene structure and signals to identify potential genes within genomic sequences.

1. **Principles:**
    
    - **Signal Sensors:** These are short sequence motifs that are indicative of specific gene features, including:
        - **Splice Sites:** The locations where RNA splicing occurs, separating exons from introns.
        - **Branch Points and Polypyrimidine Tracts:** Key sequences involved in splicing regulation.
        - **Start and Stop Codons:** Indicators of the beginning and end of a protein-coding region.
   
    - **Content Sensors:** These analyze the broader sequence context, such as:
        - **Codon Usage Patterns:** Specific to species, these patterns help distinguish coding regions from non-coding regions using statistical models.
        - **Nucleotide Composition:** Variations in GC content, repetitive elements, and sequence complexity that differ between coding and non-coding regions.

2. **Algorithms and Models:**
    
    - **Dynamic Programming:** A computational approach that breaks down complex problems into simpler sub-problems, used for optimizing gene structure predictions.
    - **Linear Discriminant Analysis:** A statistical method that classifies sequences as coding or non-coding based on multiple features.
    - **Linguistic Methods:** These apply principles of language processing to predict gene structures by recognizing patterns similar to grammatical rules in sequences.
    - **Hidden Markov Models (HMMs):** Probabilistic models that capture sequence features and their dependencies, widely used for gene prediction.
    - **Neural Networks:** Machine learning models that learn complex patterns in data, allowing for the prediction of gene structures with high accuracy.

3. **Ab Initio Gene Prediction Software:**

- **Eukaryotes:**
     
| Year | Tool Name | Publication | Method | Organism |
| ---- | --------- | ----------- | ------ | -------- |
| 1994 | GeneLang | Dong, S., & Searls, D. B. (1994). Gene structure prediction by linguistic methods. Genomics, 23(3), 540–551. https://doi.org/10.1006/geno.1994.1541.<br> | Linguistic method, Hidden Markov Model, Dynamic Programming, Weight Array Matrix | Eukaryote |
| 1995 | Fgeneh / GeneFinder | Solovyev, V. V., Salamov, A. A., & Lawrence, C. B. (1995). Identification of human gene structure using linear discriminant functions and dynamic programming. Proceedings. International Conference on Intelligent Systems for Molecular Biology, 3, 367–375.<br> | Hidden Markov Model, Dynamic Programming, Linear Discriminant Analysis | Human |
| 1997 | HMMGene | Krogh, A. (1997). Two methods for improving performance of an HMM and their application for gene finding. Proceedings. International Conference on Intelligent Systems for Molecular Biology, 5, 179–186.<br> | Hidden Markov Model Conditional Random Field | Vertebrate and C. elegans |
| 1997 | GenScan | Burge, C., & Karlin, S. (1997). Prediction of complete gene structures in human genomic DNA. Journal of Molecular Biology, 268(1), 78–94. https://doi.org/10.1006/jmbi.1997.0951.<br> | Generalized Hidden Markov Model ||
| 2003 | AUGUSTUS | Stanke, M., & Waack, S. (2003). Gene prediction with a hidden Markov model and a new intron submodel. Bioinformatics (Oxford, England), 19 Suppl 2, ii215-225. https://doi.org/10.1093/bioinformatics/btg1080<br> | Hidden Markov Model | Eukaryote |
| 2004 | GlimmerHMM | Majoros, W. H., Pertea, M., & Salzberg, S. L. (2004). TigrScan and GlimmerHMM: Two open source ab initio eukaryotic gene-finders. Bioinformatics (Oxford, England), 20(16), 2878–2879. https://doi.org/10.1093/bioinformatics/bth315<br> | Generalized Hidden Markov Model | Eukaryote |
| 2005 | GeneMark-ES | Lomsadze, A., Ter-Hovhannisyan, V., Chernoff, Y. O., & Borodovsky, M. (2005). Gene identification in novel eukaryotic genomes by self-training algorithm. Nucleic Acids Research, 33(20), 6494–6506. https://doi.org/10.1093/nar/gki937.<br> | | Eukaryote |
| 2005 | BGF (Beijing Gene Finder) | Li, H., Liu, J.-S., Xu, Z., Jin, J., Fang, L., Gao, L., Li, Y.-D., Xing, Z.-X., Gao, S.-G., Liu, T., Li, H.-H., Li, Y., Fang, L.-J., Xie, H.-M., Zheng, W.-M., & Hao, B.-L. (2005). Test Data Sets and Evaluation of Gene Prediction Programs on the Rice Genome. Journal of Computer Science and Technology, 20(4), 446–453. https://doi.org/10.1007/s11390-005-0446-x.<br> | Semi Hidden Markov Model | Eukaryote |
| 2010 | Gnomon | Souvorov, A., Kapustin, Y., Kiryutin, B., Chetvernin, V., Tatusova, T., & Lipman, D. (2010). Gnomon – NCBI eukaryotic gene prediction tool. Available in: <https://www.ncbi.nlm.nih.gov/core/assets/genome/files/Gnomon-description.pdf> <br>| Hidden Markov Model, Weight Array Matrix| |
| 2014 | GeneMark-ET | Lomsadze, A., Burns, P. D., & Borodovsky, M. (2014). Integration of mapped RNA-Seq reads into automatic training of eukaryotic gene finding algorithm. Nucleic Acids Research, 42(15), e119. https://doi.org/10.1093/nar/gku557.<br> | Hidden Markov Model | Eukaryote |

- **Prokaryotes:**

| Year | Tool Name | Publication | Method | Organism |
| ---- | --------- | ----------- | ------ | -------- |
| 1998 | GeneMark.hmm | Lukashin, A. V., & Borodovsky, M. (1998). GeneMark.hmm: New solutions for gene finding. Nucleic Acids Research, 26(4), 1107–1115. https://doi.org/10.1093/nar/26.4.1107.<br> | Hidden Markov Model | Prokaryote / Archaea |
| 2001 | GeneHacker Plus | Yada, T., Totoki, Y., Takagi, T., & Nakai, K. (2001). A novel bacterial gene-finding system with improved accuracy in locating start codons. DNA Research: An International Journal for Rapid Publication of Reports on Genes and Genomes, 8(3), 97–106. https://doi.org/10.1093/dnares/8.3.97.<br> | Hidden Markov Model | Prokaryote |
| 2001 | GeneMarkS | Besemer, J., Lomsadze, A., & Borodovsky, M. (2001). GeneMarkS: A self-training method for prediction of gene starts in microbial genomes. Implications for finding sequence motifs in regulatory regions. Nucleic Acids Research, 29(12), 2607–2618. https://doi.org/10.1093/nar/29.12.2607.<br> | Hidden Markov Model | Prokaryote |
| 2003 | AMIGene | Bocs, S., Cruveiller, S., Vallenet, D., Nuel, G., & Médigue, C. (2003). AMIGene: Annotation of MIcrobial Genes. Nucleic Acids Research, 31(13), 3723–3726.<br> | Hidden Markov Model | Prokaryote |
| 2003 | EasyGene | Larsen, T. S., & Krogh, A. (2003). EasyGene – a prokaryotic gene finder that ranks ORFs by statistical significance. BMC Bioinformatics, 4(1), 21. https://doi.org/10.1186/1471-2105-4-21.<br> | Hidden Markov Model | Prokaryote / Archaea |
| 2003 | ZCurve | Guo, F.-B., Ou, H.-Y., & Zhang, C.-T. (2003). ZCURVE: A new system for recognizing protein-coding genes in bacterial and archaeal genomes. Nucleic Acids Research, 31(6), 1780–1789. https://doi.org/10.1093/nar/gkg254.<br> | Z-curve | Prokaryote / Archaea |
| 2006 | MetaGeneAnnotator (MGA) | Noguchi, H., Taniguchi, T., & Itoh, T. (2008). MetaGeneAnnotator: Detecting Species-Specific Patterns of Ribosomal Binding Site for Precise Gene Prediction in Anonymous Prokaryotic and Phage Genomes. DNA Research, 15(6), 387–396. https://doi.org/10.1093/dnares/dsn027.<br> | | Prokaryote |
| 2007 | GISMO | Krause, L., McHardy, A. C., Nattkemper, T. W., Pühler, A., Stoye, J., & Meyer, F. (2007). GISMO--gene identification using a support vector machine for ORF classification. Nucleic Acids Research, 35(2), 540–549. https://doi.org/10.1093/nar/gkl1083.<br> | Support Vector Machine | Prokaryote |
| 2010 | Prodigal (PROkaryotic DYnamic programming Gene-finding ALgorithm) | Hyatt, D., Chen, G.-L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: Prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics, 11(1), 119. https://doi.org/10.1186/1471-2105-11-119. <br>| Dynamic Programming + Hidden Markov Model | Prokaryote |
| 2014 | ZUPLS | Song, K., Tong, T., & Wu, F. (2014). Predicting essential genes in prokaryotic genomes using a linear method: ZUPLS. Integrative Biology: Quantitative Biosciences from Nano to Macro, 6(4), 460–469. https://doi.org/10.1039/c3ib40241j.<br> | Z-curve | Prokaryote |

- **Virus:**

| Year | Tool Name | Publication | Method | Organism |
| ---- | --------- | ----------- | ------ | -------- |
| 2003 | GeneMarkS (virus version) | Mills, R., Rozanov, M., Lomsadze, A., Tatusova, T., & Borodovsky, M. (2003). Improving gene annotation of complete viral genomes. Nucleic Acids Research, 31(23), 7041–7055. https://doi.org/10.1093/nar/gkg878.<br>  | Hidden Markov Model | Virus |
| 2006 | ZCURVE_V | Guo, F.-B., & Zhang, C.-T. (2006). ZCURVE_V: A new self-training system for recognizing protein-coding genes in viral and phage genomes. BMC Bioinformatics, 7(1), 9. https://doi.org/10.1186/1471-2105-7-9. <br>  | Z-curve | Virus |
| 2007 | GLIMMER3 | Delcher, A. L., Bratke, K. A., Powers, E. C., & Salzberg, S. L. (2007). Identifying bacterial genes and endosymbiont DNA with Glimmer. Bioinformatics, 23(6), 673–679. https://doi.org/10.1093/bioinformatics/btm009. <br>  | Interpolated Markov Model | Bacteria, Archaea, and Viruses |
| 2019 | Vgas (Viral Genome Annotation System) | Zhang, K.-Y., Gao, Y.-Z., Du, M.-Z., Liu, S., Dong, C., & Guo, F.-B. (2019). Vgas: A Viral Genome Annotation System. Frontiers in Microbiology, 10. https://doi.org/10.3389/fmicb.2019.00184.<br>  | ZCURVE_V + BLASTp | Virus |


- **Metagenome:**

| Year | Tool Name | Publication | Method | Organism |
| ---- | --------- | ----------- | ------ | -------- |
| 2006 | Metagene | Noguchi, H., Park, J., & Takagi, T. (2006). MetaGene: Prokaryotic gene finding from environmental genome shotgun sequences. Nucleic Acids Research, 34(19), 5623–5630. https://doi.org/10.1093/nar/gkl723.<br>  | | Metagenome |
| 2010 | FragGeneScan | Rho, M., Tang, H., & Ye, Y. (2010). FragGeneScan: Predicting genes in short and error-prone reads. Nucleic Acids Research, 38(20), e191. https://doi.org/10.1093/nar/gkq747. <br>  | Hidden Markov Model | Metagenome |
| 2010 | MetaGeneMark | Zhu, W., Lomsadze, A., & Borodovsky, M. (2010). Ab initio gene identification in metagenomic sequences. Nucleic Acids Research, 38(12), e132. https://doi.org/10.1093/nar/gkq275.<br>  | Hidden Markov Model | Metagenome |
| 2013 | MetaGUN | Liu, Y., Guo, J., Hu, G., & Zhu, H. (2013). Gene prediction in metagenomic fragments based on the SVM algorithm. BMC Bioinformatics, 14(5), S12. https://doi.org/10.1186/1471-2105-14-S5-S12.<br>  | Support Vector Machine | Metagenome |


3. **Advantages:**
    
    - **Independence from Known Sequences:** Ab initio methods can predict novel genes that have no homologs in existing databases, making them important tools for annotating newly sequenced genomes.
    - **Versatility:** Applicable to a wide range of organisms, including those with few or no closely related reference genomes.

4. **Limitations:**
    
    - **Lower Accuracy for Complex Genomes:** Ab initio predictions can be less accurate for genomes with complex structures, such as those with extensive alternative splicing or overlapping genes.
    - **High False Positive Rates:** Without homology evidence, ab initio methods may generate predictions that do not correspond to functional genes, necessitating further validation.

---

### Homology-Based Approach

The homology-based approach utilizes evolutionary conservation to predict genes by comparing unknown sequences with those from other organisms with known annotations. This method is grounded in the principle that functional regions, like exons, are more conserved across species than non-functional regions, such as intergenic or intronic sequences. By identifying similarities between sequences, this approach infers the structure and potential function of unknown genomic regions.

#### Principles:

1. **Sequence Conservation:** Functional regions (e.g., exons) are typically conserved across species, making sequence similarity a robust indicator of gene presence and structure.
2. **Homology Inference:** By aligning sequences from expressed sequence tags (ESTs), proteins, or well-annotated genomes, researchers can extrapolate gene structures in the target genome, effectively transferring annotations from model organisms.

#### Homology Based Gene Prediction Software:

| Year | Tool Name | Publication | Method | Organism |
| ---- | --------- | ----------- | ------ | -------- |
| 1999 | CRITICA (Coding Region Identification Tool Invoking Comparative Analysis) | Badger, J. H., & Olsen, G. J. (1999). CRITICA: Coding region identification tool invoking comparative analysis. Molecular Biology and Evolution, 16(4), 512–524. https://doi.org/10.1093/oxfordjournals.molbev.a026133.<br> | Comparative| Prokaryote / Archaea |
| 2000 | CEM | Bafna, V., & Huson, D. H. (2000). The conserved exon method for gene finding. Proceedings. International Conference on Intelligent Systems for Molecular Biology, 8, 3–12.<br>| Comparative genomics | |
| 2000 | Rosetta | Batzoglou, S., Pachter, L., Mesirov, J. P., Berger, B., & Lander, E. S. ([s.d.]). Human and Mouse Gene Structure: Comparative Analysis and Application to Exon Prediction.<br>| Comparative genomics | |
| 2001 | SGP-1 (Syntenic Gene Prediction) | Wiehe, T., Gebauer-Jung, S., Mitchell-Olds, T., & Guigó, R. (2001). SGP-1: Prediction and validation of homologous genes based on sequence alignments. Genome Research, 11(9), 1574–1583. https://doi.org/10.1101/gr.177401.<br> | Comparative | Vertebrates and plants |
| 2001 | GenomeScan | Yeh, R. F., Lim, L. P., & Burge, C. B. (2001). Computational inference of homologous gene structures in the human genome. Genome Research, 11(5), 803–816. https://doi.org/10.1101/gr.175701.<br>| Comparative | |
| 2001 | Twinscan | Korf, I., Flicek, P., Duan, D., & Brent, M. R. (2001). Integrating genomic homology into gene structure prediction. Bioinformatics, 17(suppl_1), S140–S148. https://doi.org/10.1093/bioinformatics/17.suppl_1.S140.<br>| Comparative-genomics-based | |
| 2002 | AGenDA (Alignment-based Gene-Detection Algorithm) | Rinner, O., & Morgenstern, B. (2002). AGenDA: Gene prediction by comparative sequence analysis. In Silico Biology, 2(3), 195–205.<br>| Comparative | Eukaryote |
| 2002 | DOUBLESCAN | Meyer, I. M., & Durbin, R. (2002). Comparative ab initio prediction of gene structures using pair HMMs. Bioinformatics (Oxford, England), 18(10), 1309–1318. https://doi.org/10.1093/bioinformatics/18.10.1309.<br> | Comparative | |
| 2002 | GAZE | Howe, K. L., Chothia, T., & Durbin, R. (2002). GAZE: A generic framework for the integration of gene-prediction data by dynamic programming. Genome Research, 12(9), 1418–1427. https://doi.org/10.1101/gr.149502.<br>| Comparative / combiner | |
| 2003 | SGP2 | Parra, G., Agarwal, P., Abril, J. F., Wiehe, T., Fickett, J. W., & Guigó, R. (2003). Comparative Gene Prediction in Human and Mouse. Genome Research, 13(1), 108–117. https://doi.org/10.1101/gr.871403.<br>| Comparative | Eukaryote |
| 2003 | SLAM | Alexandersson, M., Cawley, S., & Pachter, L. (2003). SLAM: Cross-species gene finding and alignment with a generalized pair hidden Markov model. Genome Research, 13(3), 496–502. https://doi.org/10.1101/gr.424203.<br> | Comparative | Eukaryote |
| 2003 | ETOPE | Nekrutenko, A., Chung, W.-Y., & Li, W.-H. (2003). ETOPE: Evolutionary test of predicted exons. Nucleic Acids Research, 31(13), 3564–3567. https://doi.org/10.1093/nar/gkg597.<br>| Comparative / evolutionary | Eukaryote |
| 2003 | EvoGene | Pedersen, J. S., & Hein, J. (2003). Gene finding with a hidden Markov model of genome structure and evolution. Bioinformatics (Oxford, England), 19(2), 219–227. https://doi.org/10.1093/bioinformatics/19.2.219.<br>| Comparative / evolutionary | |
| 2004 | Projector | Meyer, I. M., & Durbin, R. (2004). Gene structure conservation aids similarity based gene prediction. Nucleic Acids Research, 32(2), 776–783. https://doi.org/10.1093/nar/gkh211.<br>| Comparative | |
| 2005 | TWAIN | Majoros, W. H., Pertea, M., & Salzberg, S. L. (2005). Efficient implementation of a generalized pair hidden Markov model for comparative gene finding. Bioinformatics (Oxford, England), 21(9), 1782–1788. https://doi.org/10.1093/bioinformatics/bti297.<br> | Comparative | |
| 2006 | DOGFISH (for ‘detection of genomic features in sequence homologies’) | Carter, D., & Durbin, R. (2006). Vertebrate gene finding from multiple-species alignments using a two-level strategy. Genome Biology, 7(Suppl 1), S6. https://doi.org/10.1186/gb-2006-7-s1-s6.<br>| Comparative | Vertebrate |
| 2006 | N_Scan_EST | Wei, C., & Brent, M. R. (2006). Using ESTs to improve the accuracy of de novo gene prediction. BMC Bioinformatics, 7, 327. https://doi.org/10.1186/1471-2105-7-327.<br>| Comparative + Evidence | |
| 2006 | TWINSCAN_EST | Wei, C., & Brent, M. R. (2006). Using ESTs to improve the accuracy of de novo gene prediction. BMC Bioinformatics, 7, 327. https://doi.org/10.1186/1471-2105-7-327.<br> | Comparative + Evidence | |
| 2007 | Contrast | Gross, S. S., Do, C. B., Sirota, M., & Batzoglou, S. (2007). CONTRAST: A discriminative, phylogeny-free approach to multiple informant de novo gene prediction. Genome Biology, 8(12), R269. https://doi.org/10.1186/gb-2007-8-12-r269.<br>| Comparative | |
| 2007 | Conrad | DeCaprio, D., Vinson, J. P., Pearson, M. D., Montgomery, P., Doherty, M., & Galagan, J. E. (2007). Conrad: Gene prediction using conditional random fields. Genome Research, 17(9), 1389–1398. https://doi.org/10.1101/gr.6558107.<br>| Comparative | |
| 2015 | GASS ( Genome Annotation based on Species Similarity) | Wang, Y., Chen, L., Song, N., & Lei, X. (2015). GASS: Genome structural annotation for Eukaryotes based on species similarity. BMC Genomics, 16(1), 150. https://doi.org/10.1186/s12864-015-1353-3.<br>| Comparative | |
| 2016 | AugustusCGP | König, S., Romoth, L. W., Gerischer, L., & Stanke, M. (2016). Simultaneous gene finding in multiple genomes. Bioinformatics (Oxford, England), 32(22), 3388–3395. https://doi.org/10.1093/bioinformatics/btw494.<br>| Comparative | Eukaryote |
| 2016 | CESAR | Sharma, V., & Hiller, M. (2019). Coding Exon-Structure Aware Realigner (CESAR): Utilizing Genome Alignments for Comparative Gene Annotation. Methods in Molecular Biology (Clifton, N.J.), 1962, 179–191. https://doi.org/10.1007/978-1-4939-9173-0_10.<br>| Comparative | |



#### Advantages:

- **High Accuracy:** This approach achieves high accuracy when homologous sequences are available, enabling precise gene models based on existing annotations.
- **Functional Insights:** By aligning unknown sequences to known references, similarity-based methods not only predict gene structures but also provide functional annotations, such as identifying protein-coding genes, motifs, and other regulatory elements.

#### Limitations:

- **Dependency on Known Data:** The success of similarity-based methods is limited by the availability and quality of homologous sequences in databases. Novel or highly divergent genes may be overlooked if no similar sequences are present.
- **Low Resolution in Non-Conserved Regions:** Non-conserved regions, such as species-specific genes or non-coding RNAs, are often missed due to insufficient sequence similarity.

---

#### **Integrated Approaches**

To overcome the limitations of each individual method, integrated approaches combine similarity-based and ab initio predictions, along with experimental data (e.g., RNA-Seq, ESTs), to enhance the accuracy and reliability of gene annotations.

- **Hybrid Models:** Tools like MAKER, Ensembl, and NCBI’s RefSeq pipelines utilize integrated models that combine ab initio predictions with sequence homology and evidence-based data.
- **Evidence Integration:** Incorporating transcriptome data (e.g., RNA-Seq) provides direct evidence of gene expression, improving the prediction of exon-intron structures and alternative splicing events.
- **Iterative Refinement:** Combining predictions from multiple methods allows for iterative refinement, reducing false positives and enhancing the resolution of gene models.

Gene prediction remains a dynamic and evolving field, with ongoing improvements driven by advances in computational methods, machine learning, and the increasing availability of genomic and transcriptomic data. These methods collectively enable researchers to construct detailed and accurate maps of genomes, supporting the study of gene function, evolution, and the complex regulatory networks that underpin biological systems.

#### Integrated Based Gene Prediction Software:

| Year | Tool Name | Publication | Method | Organism |
| ---- | --------- | ----------- | ------ | -------- |
| 2003 | GeneComber | Shah, S. P., McVicker, G. P., Mackworth, A. K., Rogic, S., & Ouellette, B. F. F. (2003). GeneComber: Combining outputs of gene prediction programs for improved results. Bioinformatics (Oxford, England), 19(10), 1296–1297. https://doi.org/10.1093/bioinformatics/btg139.<br>| EUI, GI and EUI frame algorithms | |
| 2004 | Combiner | Allen, J. E., Pertea, M., & Salzberg, S. L. (2004). Computational Gene Prediction Using Multiple Sources of Evidence. Genome Research, 14(1), 142–148. https://doi.org/10.1101/gr.1562804.<br>| Linear Combiner with voting function | |
| 2004 | Reganor | McHardy, A. C., Goesmann, A., Pühler, A., & Meyer, F. (2004). Development of joint application strategies for two microbial gene finders. Bioinformatics, 20(10), 1622–1631. https://doi.org/10.1093/bioinformatics/bth137.<br>| Combiner: abinitio + evidence           | Prokaryote / Archaea |
| 2005 | JIGSAW | Allen, J. E., & Salzberg, S. L. (2005). JIGSAW: Integration of multiple sources of evidence for gene prediction. Bioinformatics, 21(18), 3596–3603. https://doi.org/10.1093/bioinformatics/bti609.<br>| GHMM-like algorithm | |
| 2007 | GLEAN | Elsik, C. G., Mackey, A. J., Reese, J. T., Milshina, N. V., Roos, D. S., & Weinstock, G. M. (2007). Creating a honey bee consensus gene set. Genome Biology, 8(1), R13. https://doi.org/10.1186/gb-2007-8-1-r13.<br>| HMM | Eukaryote |
| 2007 | Genomix | Coghlan, A., & Durbin, R. (2007). Genomix: A method for combining gene-finders’ predictions, which uses evolutionary conservation of sequence and intron-exon structure. Bioinformatics (Oxford, England), 23(12), 1468–1475. https://doi.org/10.1093/bioinformatics/btm133.<br> | DP | Eukaryote |
| 2008 | SCGPred | Li, X., Ren, Q., Weng, Y., Cai, H., Zhu, Y., & Zhang, Y. (2008). SCGPred: A score-based method for gene structure prediction by combining multiple sources of evidence. Genomics, Proteomics & Bioinformatics, 6(3–4), 175–185. https://doi.org/10.1016/S1672-0229(09)60005-X.<br> | | Eukaryote |
| 2008 | Evidence Modeler (EVM) | Haas, B. J., Salzberg, S. L., Zhu, W., Pertea, M., Allen, J. E., Orvis, J., White, O., Buell, C. R., & Wortman, J. R. (2008). Automated eukaryotic gene structure annotation using EVidenceModeler and the Program to Assemble Spliced Alignments. Genome Biology, 9(1), R7. https://doi.org/10.1186/gb-2008-9-1-r7.<br> | | |
| 2008 | Evigan | Liu, Q., Mackey, A. J., Roos, D. S., & Pereira, F. C. N. (2008). Evigan: A hidden variable model for integrating gene evidence for eukaryotic gene prediction. Bioinformatics (Oxford, England), 24(5), 597–605. https://doi.org/10.1093/bioinformatics/btn004.<br>| Dynamic Bayes networks (DBNs) | Eukaryote |
| 2008 | Maker | Cantarel, B. L., Korf, I., Robb, S. M. C., Parra, G., Ross, E., Moore, B., Holt, C., Sánchez Alvarado, A., & Yandell, M. (2008). MAKER: An easy-to-use annotation pipeline designed for emerging model organism genomes. Genome Research, 18(1), 188–196. https://doi.org/10.1101/gr.6743907.<br> | | |
| 2012 | eCRAIG (ensemble CRAIG) | Bernal, A., Crammer, K., & Pereira, F. (2012). Automated gene-model curation using global discriminative learning. Bioinformatics (Oxford, England), 28(12), 1571–1578. https://doi.org/10.1093/bioinformatics/bts176.<br>| CRF-based | |
| 2015 | Ipred | Zickmann, F., & Renard, B. Y. (2015). IPred—Integrating ab initio and evidence based gene predictions to improve prediction accuracy. BMC Genomics, 16(1), 134. https://doi.org/10.1186/s12864-015-1315-9.<br> | Evidence-based combiner | |



#### Pipeline Based Gene Prediction Software:
| Year | Tool Name | Publication | Method | Organism |
| ---- | --------- | ----------- | ------ | -------- |
| 2003 | PASA (Program to Assemble Spliced Alignments) | Haas, B. J., Delcher, A. L., Mount, S. M., Wortman, J. R., Smith, R. K., Hannick, L. I., Maiti, R., Ronning, C. M., Rusch, D. B., Town, C. D., Salzberg, S. L., & White, O. (2003). Improving the Arabidopsis genome annotation using maximal transcript alignment assemblies. Nucleic Acids Research, 31(19), 5654–5666. https://doi.org/10.1093/nar/gkg770<br> | Pipeline - combiner - evidence-based | |
| 2006 | MaGe | Vallenet, D., Labarre, L., Rouy, Z., Barbe, V., Bocs, S., Cruveiller, S., Lajus, A., Pascal, G., Scarpelli, C., & Médigue, C. (2006). MaGe: A microbial genome annotation system supported by synteny results. Nucleic Acids Research, 34(1), 53–65. https://doi.org/10.1093/nar/gkj406.<br> | Pipeline | Bacteria |
| 2006 | FGENESH++ | Solovyev, V., Kosarev, P., Seledsov, I., & Vorobyev, D. (2006). Automatic annotation of eukaryotic genes, pseudogenes and promoters. Genome Biology, 7(1), S10. https://doi.org/10.1186/gb-2006-7-s1-s10.<br> | Pipeline hybrid: HMM + extrinsic | |
| 2008 | RAST | Aziz, R. K., Bartels, D., Best, A. A., DeJongh, M., Disz, T., Edwards, R. A., Formsma, K., Gerdes, S., Glass, E. M., Kubal, M., Meyer, F., Olsen, G. J., Olson, R., Osterman, A. L., Overbeek, R. A., McNeil, L. K., Paarmann, D., Paczian, T., Parrello, B., … Zagnitko, O. (2008). The RAST Server: Rapid annotations using subsystems technology. BMC Genomics, 9, 75. https://doi.org/10.1186/1471-2164-9-75.<br> | Pipeline | Bacterial and archaeal |
| 2011 | MAKER2 | Holt, C., & Yandell, M. (2011). MAKER2: An annotation pipeline and genome-database management tool for second-generation genome projects. BMC Bioinformatics, 12, 491. https://doi.org/10.1186/1471-2105-12-491.<br>| Pipeline / combiner | Eukaryote / Prokaryote |
| 2011 | VMGAP (TheViral MetaGenome Annotation Pipeline) | Lorenzi, H. A., Hoover, J., Inman, J., Safford, T., Murphy, S., Kagan, L., & Williamson, S. J. (2011). TheViral MetaGenome Annotation Pipeline(VMGAP):an automated tool for the functional annotation of viral Metagenomic shotgun sequencing data. Standards in Genomic Sciences, 4(3), 418–429. https://doi.org/10.4056/sigs.1694706.<br>| Pipeline | Viruses |
| 2012 | MOCAT | Kultima, J. R., Sunagawa, S., Li, J., Chen, W., Chen, H., Mende, D. R., Arumugam, M., Pan, Q., Liu, B., Qin, J., Wang, J., & Bork, P. (2012). MOCAT: A metagenomics assembly and gene prediction toolkit. PloS One, 7(10), e47656. https://doi.org/10.1371/journal.pone.0047656.<br>| Pipeline Use Prodigal or MetaGeneMark | Metagenome |
| 2014 | OMIGA (Optimized Maker-Based Insect Genome Annotation) | Liu, J., Xiao, H., Huang, S., & Li, F. (2014). OMIGA: Optimized Maker-Based Insect Genome Annotation. Molecular Genetics and Genomics: MGG, 289(4), 567–573. https://doi.org/10.1007/s00438-014-0831-7.<br>| Pipeline (MAKER) Augustus, Snap, GeneMark | Insect |
| 2014 | Prokka | Seemann, T. (2014). Prokka: Rapid prokaryotic genome annotation. Bioinformatics (Oxford, England), 30(14), 2068–2069. https://doi.org/10.1093/bioinformatics/btu153.<br>| Pipeline Ab initio + evidence-based for functional annotation | Prokaryote |
| 2014 | DFAST (DDBJ Fast Annotation and Submission Tool) | Tanizawa, Y., Fujisawa, T., & Nakamura, Y. (2018). DFAST: A flexible prokaryotic genome annotation pipeline for faster genome publication. Bioinformatics (Oxford, England), 34(6), 1037–1039. https://doi.org/10.1093/bioinformatics/btx713.<br>| Pipeline | Prokaryote |
| 2016 | PGAP | Tatusova T. et al. (2016) NCBI prokaryotic genome annotation pipeline. Nucleic Acids Res ., 44, 6614–6624. | Pipeline GenemarkS+ Glimmer + extrinsic data | Prokaryote |
| 2016 | CAT Comparative Analysis toolkit | Fiddes, I. T., Armstrong, J., Diekhans, M., Nachtweide, S., Kronenberg, Z. N., Underwood, J. G., Gordon, D., Earl, D., Keane, T., Eichler, E. E., Haussler, D., Stanke, M., & Paten, B. (2018). Comparative Annotation Toolkit (CAT)—Simultaneous clade and personal genome annotation. Genome Research, 28(7), 1029–1038. https://doi.org/10.1101/gr.233460.117.<br> | Pipeline Evidence-based, comparative-abinitio (AugustusCGP) | |
| 2017 | funannotate | Jon Love, Jon Palmer, Jason Stajich, Tyler Esser, Erik Kastman, Daniel Bogema, & David Winter. (2019). nextgenusfs/funannotate: funannotate v1.7.0 (1.7.0). Zenodo. https://doi.org/10.5281/zenodo.3534297.<br> | Pipeline Evidence Modeler + Augustus + GeneMark-ES/ET + evidence + PASA | Built specifically for fungi, but will also work with higher eukaryotes |
| 2018 | FunGAP | Min, B., Grigoriev, I. V., & Choi, I.-G. (2017). FunGAP: Fungal Genome Annotation Pipeline using evidence-based gene model evaluation. Bioinformatics (Oxford, England), 33(18), 2936–2937. https://doi.org/10.1093/bioinformatics/btx353.<br>| Pipeline | |
| 2019 | VAPiD (Viral Annotation Pipeline and iDentification) | Shean, R. C., Makhsous, N., Stoddard, G. D., Lin, M. J., & Greninger, A. L. (2019). VAPiD: A lightweight cross-platform viral annotation pipeline and identification tool to facilitate virus genome submissions to NCBI GenBank. BMC Bioinformatics, 20(1), 48. https://doi.org/10.1186/s12859-019-2606-y.<br>| Pipeline | Virus |
| 2019 | GAAP | Kong, J., Huh, S., Won, J.-I., Yoon, J., Kim, B., & Kim, K. (2019). GAAP: A Genome Assembly + Annotation Pipeline. BioMed Research International, 2019, 4767354. https://doi.org/10.1155/2019/4767354.<br>| Pipeline Augustus, EVM, MAKER, PASA | |
| 2021 | Bakta | Schwengers, O., Jelonek, L., Dieckmann, M. A., Beyvers, S., Blom, J., & Goesmann, A. (2021). Bakta: Rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11), 000685. https://doi.org/10.1099/mgen.0.000685.<br>| Pipeline | Bacteria + plasmid |
| 2021 | FINDER | Banerjee, S., Bhandary, P., Woodhouse, M., Sen, T. Z., Wise, R. P., & Andorf, C. M. (2021). FINDER: An automated software package to annotate eukaryotic genes from RNA-Seq data and associated protein sequences. BMC Bioinformatics, 22(1), 205. https://doi.org/10.1186/s12859-021-04120-9.<br>| Pipeline Protein + RNAseq alignment + Abinitio (BRAKER) | Eukaryote |
| 2022 | Pharokka | Bouras, G., Nepal, R., Houtak, G., Psaltis, A. J., Wormald, P.-J., & Vreugde, S. (2023). Pharokka: A fast scalable bacteriophage annotation tool. Bioinformatics (Oxford, England), 39(1), btac776. https://doi.org/10.1093/bioinformatics/btac776.<br>| Pipeline Prodigal is used for CDS | Bacteriophage |
| 2023 | FLAG | Troy, W., Damas, J., Titus, A. J., & Cantarel, B. L. (2023). FLAG: Find, Label, Annotate Genomes, a fully automated tool for genome gene structural and functional annotation of highly fragmented non-model species (p. 2023.07.14.548907). bioRxiv. https://doi.org/10.1101/2023.07.14.548907.<br>| Pipeline Protein + RNAseq alignment + Abinitio | Eukaryote |
| 2023 | Galba | Brůna, T., Li, H., Guhlin, J., Honsel, D., Herbold, S., Stanke, M., Nenasheva, N., Ebel, M., Gabriel, L., & Hoff, K. J. (2023). Galba: Genome annotation with miniprot and AUGUSTUS. BMC Bioinformatics, 24(1), 327. https://doi.org/10.1186/s12859-023-05449-z.<br>| Pipeline Hybrid = Ab-initio evidence-driven | Eukaryote |


---
