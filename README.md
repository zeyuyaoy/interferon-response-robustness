### Context Dependence and Measurement Sensitivity of Interferon-β Induction Contrasts in Human Monocytes and NK Cells

Interferon (IFN) responses are frequently compared between immune cell types by scoring a set of IFN-stimulated genes in each population. Quantitative comparisons between cell types depend on normalization, log transformation and sequencing depth, as well as induction of the scored genes.

We examined the extent to which one such comparison depends on analytical choices. We analysed single-cell RNA-sequencing data from peripheral blood cells of eight donors with systemic lupus erythematosus (SLE). Matched aliquots from each donor were either left untreated or stimulated with IFN-β for 6 h (GEO GSE96583). Counts from CD14$`^+`$ monocytes and natural killer (NK) cells were aggregated by donor, cell type and condition into pseudobulk profiles, and the induction of a published ten-gene IFN panel was quantified.

We found that mean induction was 1.834 $`\log_2`$ units greater in monocytes than in NK cells, and greater in monocytes in every donor. Every alternative analysis preserved the positive mean. Its magnitude ranged from 1.23 to 1.87, and one donor’s sign changed across specifications. The most influential choices were how low counts were regularized before log transformation (the prior count, or a fixed offset in counts per million) and the number of cells aggregated per pseudobulk. In an independent cohort of five healthy donors whose monocytes and NK cells were enriched and cultured separately for 21 h (HIRISA, GEO GSE306664), the contrast was reversed, coinciding with high panel expression in untreated monocytes of four of the five donors; the regularization of low counts changed the contrast by only 0.155, as expected for pseudobulk libraries that were 94-378 times larger. The depositors' own differential-expression results for this cohort also show greater induction in NK cells for all ten genes.

Cell-type-specific IFN responses in this dataset were reported previously ([Kang et al., 2018](https://doi.org/10.1038/nbt.4042)); our study further evaluated three properties of the comparison:

1) its average direction,
2) its magnitude, and
3) its consistency across individual donors.

#### Study design

We used `muscData::Kang18_8vs8()` (ExperimentHub EH2259; GSE96583), comprising 35,635 genes and 29,065 droplets. Cells from each of eight SLE donors were cultured for 6 h with or without 100 U/mL IFN-β. Cells from all donors were then pooled by condition and sequenced ([Kang et al., 2018](https://doi.org/10.1038/nbt.4042)). "Control" denotes the untreated aliquot cultured in parallel with the stimulated aliquot.

We retained cells annotated as singlets and as CD14$`^+`$ monocytes or NK cells, with at least 200 detected genes. The final dataset comprised 7,346 cells in 32 groups (8 donors × 2 cell types × 2 conditions). The smallest group, the control NK cells of donor 1039, contained 18 cells.

Counts were summed within each group to generate 32 pseudobulk profiles. Monocyte and NK profiles were normalized separately using the trimmed mean of M-values (TMM) method in edgeR, which adjusts for differences in library size and composition between samples ([Robinson & Oshlack, 2010](https://doi.org/10.1186/gb-2010-11-3-r25); [Chen et al., 2025](https://doi.org/10.1093/nar/gkaf018)). Lowly expressed genes were excluded only during estimation of the normalization factors; all genes were retained for scoring. Expression was expressed as $`\log_2`$ counts per million (CPM), with a prior count $`p`$ added before transformation to avoid undefined values for zero counts ($`p = 2`$ in the primary analysis).

The panel comprised IFI44L, ISG15, IFIT3, XAF1, MX1, IFI6, IFIT1, TRIM22, MX2 and RSAD2. These genes were defined in an independent study by their consistent induction across type I IFNs and cell types ([Rigby et al., 2026](https://doi.org/10.1038/s44319-026-00750-3)). No gene was selected or excluded on the basis of the present data.

All analyses are specified in `analysis/primary_contrast/scripts/config.R` and the files in `analysis/inputs/`:
- a lower prior count (0.5);
- sequential exclusion of each donor;
- 100 random subsamples of 18 cells per group, equalizing cell numbers across groups;
- exclusion of cells above the 99.5th percentile of UMI counts within each group, as a safeguard against doublets;
- comparison with random, expression-matched panels of non-IFN genes;
- estimation of TMM factors without IFN-related genes;
- a revised identifier mapping of the Hallmark IFN-α gene set;
- the Hallmark set excluding the genes shared with the panel;
- ordering of donors by baseline expression;
- an alternative expression measurement, $`\log_2(\mathrm{CPM}+1)`$ with a fixed offset of 1 CPM (`fixed_offset_cpm`; `analysis/inputs/fixed_offset_specification.tsv`).

#### Results

For donor $`d`$ and cell type $`\ell`$ (monocyte or NK), **induction** $`R_{d\ell}`$ was defined as the mean difference between stimulated and control $`\log_2`$ CPM across the ten panel genes. The **contrast** $`\delta_d`$ is the difference between monocyte and NK induction within a donor; $`\theta`$ is the mean contrast across the eight donors:

$$
R_{d\ell} = \frac{1}{10} \sum_{g=1}^{10} \left( y_{dg\ell}^{\text{stim}} - y_{dg\ell}^{\text{ctrl}} \right), \qquad
\delta_d = R_{d,\text{mono}} - R_{d,\text{NK}}, \qquad
\theta = \frac{1}{8} \sum_{d=1}^{8} \delta_d .
$$

Here, $`y_{jg} = \log_2\!\left[10^6 (C_{jg} + p_j)/(E_j + 2p_j)\right]`$, where $`C_{jg}`$ is the count for gene $`g`$ in pseudobulk $`j`$. $`E_j`$ is the **effective library size**, defined as the total count multiplied by the TMM factor. $`p_j = p\,E_j/\bar{E}_\ell`$ scales the prior count by relative library size, as implemented in edgeR. Each gene contributes equally to $`R`$; each donor contributes equally to $`\theta`$, irrespective of cell number. A value of $`\theta = 1.834`$ therefore denotes a mean difference in induction of 1.834 $`\log_2`$ units on this normalized relative scale.

> [!NOTE]
> $`\theta`$ measures a contrast on a normalized relative transcript scale; absolute transcript abundance and signalling activity require additional measurements.

| Analysis                                                  | $`\hat{\theta}`$  | Interval                       | Donors with $`\delta_d > 0`$ |
|-----------------------------------------------------------|-------------------|--------------------------------|------------------------------|
| Ten-gene panel, prior count 2 (primary)                   | 1.834             | $`[1.163, 2.506]`$             | 8/8                          |
| Ten-gene panel, prior count 0.5                           | 1.595             | $`[0.905, 2.284]`$             | 7/8                          |
| Ten-gene panel, TMM factors estimated without IFN genes   | 1.853             | $`[1.171, 2.534]`$             | 8/8                          |
| Ten-gene panel, top 0.5% UMI cells per group removed (38) | 1.870             | -                              | -                            |
| Ten-gene panel, $`\log_2(\mathrm{CPM}+1)`$, fixed offset  | 1.233             | $`[0.692, 1.773]`$             | 7/8                          |
| Hallmark IFN-α set, 96 genes                              | 1.452             | $`[1.196, 1.709]`$             | 8/8                          |
| Hallmark IFN-α set excluding 5 panel genes, 91 genes      | 1.402             | $`[1.175, 1.629]`$             | 8/8                          |
| Subsamples of 18 cells per group, prior count 2           | 1.450<sup>‡</sup> | $`[1.196, 1.724]`$<sup>‡</sup> | 8/8 in 89 of 100 draws       |
| Subsamples of 18 cells per group, prior count 0.5         | 1.473<sup>‡</sup> | $`[1.036, 1.859]`$<sup>‡</sup> | 8/8 in 54 of 100 draws       |

Values are in $`\log_2`$ units. Unless otherwise indicated, intervals are 95% $`t`$-intervals across the eight donor contrasts (7 degrees of freedom).

<sup>‡ Median and 2.5th-97.5th percentiles across 100 random cell subsamples; hallmark analyses score a larger and different gene set, and their values are therefore not directly comparable with those of the ten-gene panel.</sup>

##### Monocytes exhibited greater induction than NK cells in all eight donors

All ten genes were induced in both cell types in every donor (monocyte $`R`$, 4.18-7.66; NK $`R`$, 3.50-5.05), consistent with the original definition of the panel ([Rigby et al., 2026](https://doi.org/10.1038/s44319-026-00750-3)). Induction was greater in monocytes in every donor. Donor contrasts ranged from 0.232 (donor 1016) to 2.912 (donor 107), with a mean of $`\hat\theta = 1.834`$ (SD, 0.804; 95% $`t`$-interval, 1.163-2.506).

The interval assumes that the eight donor contrasts are independent and approximately normally distributed; it does not incorporate uncertainty arising from the joint estimation of normalization factors or from the choice among alternative analytical specifications; it is therefore described as conditional. An exact one-sided sign-flip test, which compares the observed mean with the means obtained under all $`2^8`$ sign assignments of the eight contrasts, yielded $`p = 2^{-8} \approx 0.0039`$. This is the minimum attainable value with eight donors; it arises whenever all contrasts are positive. Donors were treated as the units of replication, because cells from the same donor are not statistically independent ([Squair et al., 2021](https://doi.org/10.1038/s41467-021-25960-2); [Zimmerman et al., 2021](https://doi.org/10.1038/s41467-021-21038-1)).

Studies of purified human blood cells have shown that both IFN-induced transcriptional programs ([Waddell et al., 2010](https://doi.org/10.1371/journal.pone.0009753)) and proximal type I IFN signalling ([van Boxel-Dezaire et al., 2006](https://doi.org/10.1016/j.immuni.2006.08.014)) differ between cell types. Such differences provide possible explanations for the observed contrast; mechanistic attribution requires additional measurements.

##### No single donor determined the mean contrast

To assess the influence of individual donors, each donor was excluded in turn. Exclusion of a single donor's contrast yielded means of 1.680-2.063. Re-running the complete pipeline, including normalization, without that donor yielded 1.694-2.095. The largest change upon exclusion of a contrast was 0.229, below the threshold of 0.25 specified in the analysis plan for identifying an influential donor. Donor 1016 exerted the greatest influence: its exclusion *increased* the mean to 2.063 (2.095 after re-running), because its contrast was the smallest. As all eight contrasts were positive, exclusion of a single donor could not reverse the sign of the mean. This analysis therefore assesses the stability of the magnitude.

##### The prior count altered the magnitude of the contrast and the sign for one donor

The prior count is a constant added to all counts before log transformation. Reducing it from 2 to 0.5 decreased $`\hat\theta`$ to 1.595 (0.905-2.284) and changed the contrast for donor 1016 from 0.232 to −0.034. Induction *increased* in both cell types (monocytes, 6.158 to 6.394; NK cells, 4.323 to 4.799). The contrast decreased because the increase was larger in NK cells.

This asymmetry arises because edgeR scales the prior count by library size. NK pseudobulks contained substantially fewer counts than monocyte pseudobulks (mean effective library size, 117,580 versus 793,081), principally because NK groups comprised fewer cells. A prior count of 2 consequently added 17.01 CPM to each NK gene and 2.52 CPM to each monocyte gene. A larger additive constant raises near-zero control values disproportionately and thereby reduces the estimated induction. In NK cells, this effect was most pronounced for genes with minimal expression before stimulation. Three such genes, RSAD2, IFIT1 and IFI44L, were detected in 1.2-4.5% of control cells; they accounted for 0.156 of the 0.239 reduction in $`\hat\theta`$ (65.1%). The sensitivity of log-transformed single-cell data to the prior count is well documented ([Ahlmann-Eltze & Huber, 2023](https://doi.org/10.1038/s41592-023-01814-1)). In the present data, an identical setting acted unequally on two cell types that differed in library size.

##### A fixed offset of 1 CPM reduced the magnitude of the contrast and again reversed the sign for donor 1016

To test whether the sensitivity of the contrast is associated with the library-size-dependent regularization in edgeR log-CPM, we replaced $`\log_2`$ CPM with $`\log_2(\mathrm{CPM}+1)`$, which adds a fixed offset of 1 CPM to every pseudobulk profile. The primary TMM effective library sizes, pseudobulk profiles, genes and donors were reused without refitting. Because the edgeR denominator term cancels in stimulated-minus-control differences, the two measurements differ only in the CPM added to each value: 2.52 (monocytes) and 17.01 (NK cells) in the primary analysis, and 1 in both cell types here.

With the fixed offset, $`\hat\theta`$ was 1.233 (0.692-1.773), 0.601 lower than in the primary analysis. The contrast decreased in every donor, by 0.202-0.953 (mean paired difference, −0.601; 95% $`t`$-interval, −0.862 to −0.341). NK induction increased by 0.778 and monocyte induction by 0.176. Seven donor contrasts remained positive (1.131-2.152); that of donor 1016 changed from 0.232 to −0.142, as it had at prior count 0.5. RSAD2, IFIT1 and IFI44L accounted for 0.442 of the 0.601 reduction (73.5%), almost entirely through greater NK induction; these are the same sparsely detected control NK genes that accounted for the effect of the prior count. Under the criteria in `analysis/inputs/fixed_offset_specification.tsv`, the result shows sensitivity of magnitude (change ≥ 0.25) and a donor-level change of direction confined to donor 1016; the positive mean was unchanged.

These results indicate that approximately one third of the primary magnitude depends on the larger CPM offset that edgeR applies to the smaller NK libraries. The appropriate value remains dependent on the measurement specification. An equal offset on the CPM scale corresponds to unequal regularization on the count scale (approximately 0.79 counts per monocyte and 0.12 per NK pseudobulk at the mean library size); the fixed-offset estimate therefore retains dependence on library size through count-scale regularization.

##### Equalizing cell numbers preserved the direction of the contrast and reduced its magnitude

Because monocyte groups contained more cells than NK groups, 18 cells (the size of the smallest group) were drawn from every group, and the analysis was repeated across 100 draws. Identical cells were used at both prior counts. The mean contrast was positive in all 100 draws at both prior counts, with medians of 1.450 (prior count 2) and 1.473 (prior count 0.5). All eight donor contrasts were simultaneously positive in 89 draws at prior count 2 and in 54 draws at prior count 0.5. The contrast for donor 1016 was positive in 90 and 59 draws, respectively; that of every other donor was positive in at least 95.

Subsampling reduced the estimate: the full-data value (1.834) exceeded the 97.5th percentile of the prior-count-2 draws (1.724), and induction was lower in both cell types (monocytes, 4.59; NK cells, 3.15). In subsamples, reducing the prior count had little effect (change in $`\hat\theta`$, −0.023, compared with 0.239 using all cells). Subsampling simultaneously alters the cells included, library sizes, the genes that pass filtering, the TMM factors and the CPM added by the prior count. The present analysis cannot attribute the reduction to any one of these factors. Both estimates are conditional on the number of cells aggregated.

> [!NOTE]
> All draws derive from the same donors and sequenced cells; they quantify sensitivity to cell selection within the existing biological sample.

##### Between-donor variation reflected differences in control monocyte expression

To identify the source of between-donor variation, panel expression in control and stimulated monocytes was compared across donors. All monocyte profiles were normalized jointly, so these values are directly comparable between donors. Stimulated monocyte expression was similar across donors (10.74-11.24 $`\log_2`$ CPM); control expression ranged more widely, from 3.58 to 6.66. The ranking of donors by control monocyte expression was the exact inverse of their ranking by contrast (Spearman ρ = −1). This relationship is largely a consequence of the definition: when stimulated levels are nearly constant, higher control expression necessarily corresponds to lower induction. The correlation therefore identifies the source of between-donor variation and is descriptive.

Donor 1016 had the highest control monocyte expression (6.66), the lowest monocyte induction (4.18) and the smallest contrast. It was also the only donor in which control monocytes expressed the panel at a higher level than control NK cells (+0.28). Peripheral blood cells from many patients with SLE express IFN-inducible genes at elevated levels as a result of ongoing IFN exposure in vivo ([Baechler et al., 2003](https://doi.org/10.1073/pnas.0337679100); [Bennett et al., 2003](https://doi.org/10.1084/jem.20021553)). Elevated control expression reduces estimated induction at a given stimulated expression level and may contribute to the profile of donor 1016. Evaluation of this explanation would require clinical, serum IFN or treatment data, which are not available for this dataset.

An analogous decomposition can be expressed across cell types. Averaged across donors, monocyte expression exceeded NK expression by 0.487 $`\log_2`$ units after stimulation and was 1.347 units lower in controls, giving 0.487 − (−1.347) = 1.834. This partition does not have a biological interpretation; because the two cell types were normalized separately, a constant shift in all monocyte values would alter both terms equally without changing $`\theta`$.

Individual genes differed in induction amplitude. Nine of the ten genes were induced more strongly in monocytes, with the largest differences for RSAD2 (+3.41), ISG15 (+2.91) and IFIT1 (+2.50). IFI6 was induced more strongly in NK cells (−0.37). IFI6 also exhibited the largest difference in control detection between cell types: it was detected in 23% of control monocytes and 9% of control NK cells. This observation parallels the donor-level pattern, in which higher control expression accompanied lower induction; a result obtained at the panel level should therefore not be attributed to each constituent gene.

##### A larger IFN-α gene set reproduced the direction of the contrast with reduced between-donor variability

To determine whether the result depended on the choice of ten genes, the analysis was repeated with the Hallmark IFN-α response gene set ([Liberzon et al., 2015](https://doi.org/10.1016/j.cels.2015.12.004)). The version specified in the analysis plan could not be evaluated: only 86 of the 97 gene symbols (88.7%) mapped uniquely to the dataset, below the specified minimum of 90%. In a revised analysis, symbols were mapped through HGNC and Ensembl identifiers from GRCh37, the genome build of the dataset ([Seal et al., 2023](https://doi.org/10.1093/nar/gkac888); [Yates et al., 2015](https://doi.org/10.1093/bioinformatics/btu613)). This procedure resolved 96 genes and yielded $`\hat\theta = 1.452`$ (1.196-1.709; SD, 0.307). Five of these genes are also members of the ten-gene panel; their exclusion left 91 genes and yielded 1.402 (1.175-1.629; SD, 0.271), with a minimum donor contrast of 0.790. All donor contrasts were positive in both analyses.

Between-donor standard deviation was 62-66% lower than for the ten-gene panel (0.804); the contrast for donor 1016 increased from 0.23 to 0.76. Averaging across a larger number of genes is likely to reduce the influence of a small number of sparsely detected genes. Because the 91-gene analysis excludes all panel genes, the greater induction in monocytes is not restricted to those ten genes. These analyses assess gene-set sensitivity within the same donors and experiment.

> [!CAUTION]
> Scores derived from gene sets that differ in size and composition are not on a common scale ([Andreatta & Carmona, 2021](https://doi.org/10.1016/j.csbj.2021.06.043); [Ciernik et al., 2026](https://doi.org/10.1101/gr.280926.125)), so the values 1.45 and 1.83 should not be interpreted as estimates of a single quantity.

##### The panel contrast exceeded that of all expression-matched random panels

Genes with low control expression may exhibit large apparent induction for technical reasons. To assess this possibility, the panel was compared with random non-IFN genes of similar control expression, an approach established in single-cell gene-set scoring ([Tirosh et al., 2016](https://doi.org/10.1126/science.aad0501)). For each panel gene, a substitute was drawn from the 50 non-IFN genes nearest in control expression in both cell types. This procedure was repeated to generate 1,000 random ten-gene panels. The contrasts of these panels had a median of 0.233 (2.5th-97.5th percentiles, −0.312 to 0.837), and the IFN panel contrast of 1.834 exceeded all of them.

Two considerations limit this comparison. First, the median of the random panels exceeded zero, indicating that genes outside the IFN sets also showed a small monocyte-greater difference on this scale. Second, close matches were unavailable for the three most sparsely detected panel genes (RSAD2, IFI44L and IFIT1): 85-91% of their substitutes differed by more than 1 $`\log_2`$ unit. In addition, the sampling procedure does not select uniformly among all possible matched panels. The comparison is therefore reported as a descriptive rank.

##### Estimating TMM factors without IFN-related genes changed the contrast by approximately 0.02

TMM normalization assumes that most genes are not differentially expressed between samples. Stimulation violated this assumption for a subset of genes; IFN-related genes (the panel and the Hallmark IFN-α and IFN-γ genes) accounted for 23-27% of counts in stimulated monocyte pseudobulks, compared with 6-7.5% in controls. The corresponding proportions in NK pseudobulks were 16-19% and 8.5-10%; monocyte TMM factors shifted accordingly, from 1.11-1.24 in controls to 0.80-0.96 after stimulation. TMM is designed to accommodate such shifts by trimming genes with extreme fold changes ([Robinson & Oshlack, 2010](https://doi.org/10.1186/gb-2010-11-3-r25)).

To verify that this occurred, TMM factors were re-estimated without IFN-related genes. The prior count was adjusted to add the same CPM as in the primary analysis, so that the comparison isolated the effect of the normalization factors. $`\hat\theta`$ changed by 0.018 at prior count 2 and by 0.019 at prior count 0.5. Median induction of non-IFN genes remained between −0.10 and 0.02 in every donor and cell type, consistent with the null change expected for unaffected genes. IFN-related genes therefore did not materially distort normalization. Exclusion of cells above the 99.5th percentile of UMI counts within each group (38 cells) yielded $`\hat\theta = 1.870`$.

##### In an independent cohort the contrast was reversed, and the regularization of low counts had little effect

To assess whether the direction of the contrast and its sensitivity to low-count regularization transfer beyond GSE96583, the same panel, estimand and transformations were applied to the Human IFN Response in Immune Subsets Atlas (HIRISA; GEO GSE306664; [Kuan et al., 2025](https://doi.org/10.21203/rs.3.rs-8197575/v1)). CD14$`^+`$ monocytes and NK cells were enriched separately from cryopreserved PBMC of five healthy donors, cultured for 21 h with or without 100 U/mL IFN-β, and profiled by 10x Genomics Flex. Control and stimulated samples of each donor and cell type were multiplexed into a single processing pool, allowing within-pool comparisons of stimulation. The analysis was specified in `analysis/independent_cohort/specification.tsv`, which was locked after a one-donor pilot and before stimulated or NK data of the remaining four donors were downloaded. Cells labelled by the depositors as CD14 monocytes or NK cells with at least 200 detected genes were retained (217,231 cells in 20 groups; smallest group, the control monocytes of donor 3283BW, 977 cells). Each cell type was normalized separately by TMM across the ten profiles.

| Quantity                                                        | GSE96583 (8 SLE donors) | HIRISA (5 healthy donors)   |
|-----------------------------------------------------------------|-------------------------|-----------------------------|
| Stimulation; preparation                                        | 6 h; PBMC               | 21 h; enriched cell types   |
| $`\hat\theta`$, prior count 2                                   | 1.834 (1.163, 2.506)    | −4.609 (−7.158, −2.060)     |
| $`\hat\theta`$, prior count 0.5                                 | 1.595                   | −4.616                      |
| $`\hat\theta`$, $`\log_2(\mathrm{CPM}+1)`$                      | 1.233                   | −4.454                      |
| $`D = \hat\theta_{\text{prior }2} - \hat\theta_{\text{CPM}+1}`$ | 0.601 (0.341, 0.862)    | −0.155 (−0.239, −0.070)     |
| Donors with $`\delta_d > 0`$ at prior count 2                   | 8/8                     | 0/5                         |
| Mean effective library size, monocytes / NK cells               | 793,081 / 117,580       | 74.9 million / 44.5 million |
| CPM added by prior count 2, monocytes / NK cells                | 2.52 / 17.01            | 0.027 / 0.045               |

<sup>Intervals are 95% $`t`$-intervals across donor contrasts (7 and 4 degrees of freedom); they are descriptive and conditional on the specification.</sup>

The contrast was negative in all five donors ($`\hat\theta = -4.609`$). NK induction (5.99-7.17) exceeded that in GSE96583 (3.50-5.05); monocyte induction was small in four donors (0.16-1.90). In these four donors the panel was already highly expressed in untreated monocytes: mean control expression, on the scale used for GSE96583 (TMM-normalized log-CPM, prior count 2), was 9.32-11.45, compared with 5.48 in donor 3283BW, and ISG15 reached 1,409-11,406 CPM, compared with 109. Such levels exceed the range of control monocytes in GSE96583 on the same scale (3.58-6.66), with comparability still limited by assay differences despite the common transformation (`results/cross_study_synthesis/baseline_harmonization.csv`). The specified rule for anomalous controls used a related quantity, $`\log_2(\mathrm{CPM}+1)`$ on raw library sizes (9.57-11.32 and 5.80). The rule specified to identify anomalous controls flagged no donor, because elevated control expression occurred in most donors. In donor 3283BW, the only donor with low control expression, monocyte induction was 5.86 and the contrast −0.963. Elevated control expression after 21 h of culture of isolated monocytes, accompanied by CXCL10 and IFNB1 transcripts, is compatible with IFN production during culture, donor state, enrichment and assay effects. Direct IFN measurements would be required to distinguish these explanations.

On average, HIRISA pseudobulks contained 94 (monocytes) and 378 (NK cells) times more counts than those of GSE96583, and the prior count of 2 consequently added only 0.027 CPM to monocyte and 0.045 CPM to NK profiles. Reducing the prior count to 0.5 changed $`\hat\theta`$ by 0.007. Under $`\log_2(\mathrm{CPM}+1)`$, $`\hat\theta`$ was −4.454; the difference $`D = -0.155`$ (donor values −0.27 to −0.09) was below the threshold of 0.25, and no donor changed sign, which meets the prespecified criterion for agreement. $`D`$ is negative here because the fixed offset of 1 CPM now exceeds the offset applied by edgeR. As in GSE96583, its effect arose principally from NK genes that were sparsely detected in controls: IFIT1 and IFI44L (detected in 0.2-3.3% of control NK cells) accounted for 0.097 of the 0.155. In 100 subsamples of 977 cells per group, the median $`\hat\theta`$ was −4.571 (2.5th-97.5th percentiles, −4.663 to −4.474) and the median $`D`$ −0.081.

HIRISA produced a contrast with a different direction and magnitude. Differences in donor health, stimulation duration, cell preparation and sequencing platform make this a transferability assessment with unresolved causal attribution. The regularization behaved as its definition implies: the CPM added by a prior count is inversely proportional to library size, so its effect was small when both cell types had large libraries. The empirical observation common to both cohorts is narrower: the effect was concentrated in NK genes with sparse control expression.

##### Final triangulation with the published HIRISA differential-expression results

The HIRISA investigators released per-gene summaries of donor-level MAST comparisons (median and least favourable log$`_2`$ fold change and adjusted p-value across the five donors). Their monocyte population is identical to the CD14 monocytes analysed here, and their NK population is the same. Both methods found all ten genes induced in both cell types, and greater induction in NK cells for 10 of 10 genes. Both also identified weak or absent monocyte induction in at least one donor for the same five genes (IFI44L, MX1, MX2, TRIM22 and XAF1). Pseudobulk magnitudes were larger, by 1.17-2.68-fold in NK cells, and the excess was greatest for genes sparsely detected in control NK cells. Agreement between the two methods supports analytical robustness of the reversal within the same cohort. Shared data and design limitations still constrain interpretation. A median of donor-level MAST fold changes has a different scale from $`\hat\theta`$. Published values were used solely for comparison of per-gene summaries. Per-gene tables and the comparison figure are in `results/independent_cohort/triangulation/`.

##### Exploratory time course in mixed PBMCs (GSE226572)

To examine whether exposure duration alone moves the contrast, the same panel, estimand and transformations were applied to GSE226572 ([de Cevins et al., 2023](https://doi.org/10.1016/j.xcrm.2023.101333)), in which PBMCs of three healthy donors were exposed to 1,000 U/mL IFN-β for 1-36 h and profiled with 10x Single Cell 3' v3 chemistry. This is a new exploratory analysis, separate from the study above: cells were labelled with a fixed marker rule that excludes IFN genes, and the protocol (`analysis/time_course/specification.tsv`) was locked after a feasibility gate and before any panel-gene value was computed. Its primary quantity is the change in the contrast between 4 h and 24 h of exposure, $`\Gamma`$. The contrast was −0.173 at 4 h and 0.382 at 24 h, so that Γ = 0.555; the change was positive in all three donors (0.44-0.67) and under every prespecified sensitivity analysis, and the result was inconsistent with the prespecified hypothesis of a decrease. The increase reflected waning NK induction (from 4.24 to 3.73); monocyte induction was stable, consistent with the kinetics described by the original authors. Longer exposure increased the contrast in this design. Estimates are descriptive given the three-donor sample. Used alone as the reference, each of the two untreated libraries per donor gave Γ of 0.555 and 0.553 (pooled, 0.555), and widening the monocyte definition to CD14 and CD16 monocytes gave 0.595. The marker rule lost CD14 with exposure: of the cells that the depositors labelled classical monocytes at 24 h, the rule labelled 453 CD16 monocytes and 442 CD14 monocytes. In a prespecified check (`analysis/time_course/annotation_validation/specification.tsv`, locked before any expression value of these cells was read), the depositors' labels were read from the metadata of their deposited Seurat object and matched to all 115,503 of their cells. With the depositors' classical monocytes and NK cells, the contrast was −0.139 at 4 h and 0.458 at 24 h, so that Γ = 0.597, positive in all three donors (0.39-0.77); including non-classical monocytes gave 0.553. The increase was preserved across all four tested annotation definitions (`results/time_course/annotation_validation/gamma_by_annotation.csv`). Two limitations apply. Only three donors were available, so no donor-level test can reach conventional significance (the smallest possible two-sided sign-flip $`p`$ is 0.25). The nature of the second untreated library of each donor is undocumented. Per-step outputs are in `results/time_course/`, and the feasibility assessment that preceded the analysis (20 feasibility gates and 2 annotation-rule checks passed) is recorded in `results/time_course/feasibility_checks.csv`.

##### One fixed Hallmark IFN-α set across the three cohorts

To determine whether the HIRISA reversal and the rise of the contrast in GSE226572 depended on the ten-gene panel, a single Hallmark IFN-α gene set was applied to all three cohorts. The 93 genes that resolve in GSE96583 and occur as exactly one feature in both HIRISA and GSE226572 were used, and a second set excluded the five panel genes that are Hallmark members (88 genes). The specification (`analysis/hallmark_gene_set/specification.tsv`) was locked before any Hallmark value was computed in HIRISA or GSE226572, and the analysis was run once. The previously reported 96-gene GSE96583 result was known when this secondary analysis was specified.

| Cohort         | Estimand               | Ten-gene panel          | Hallmark, 93 genes | Hallmark without panel genes, 88 |
|----------------|------------------------|-------------------------|--------------------|----------------------------------|
| GSE96583 (6 h) | $`\hat\theta`$         | 1.834 (8/8 donors > 0)  | 1.478 (8/8)        | 1.427 (8/8)                      |
| HIRISA (21 h)  | $`\hat\theta`$         | −4.609 (5/5 donors < 0) | −0.990 (4/5)       | −0.758 (4/5)                     |
| GSE226572      | $`\Gamma`$             | 0.555 (3/3 donors > 0)  | 0.322 (3/3)        | 0.312 (3/3)                      |
| GSE226572      | $`\hat\theta`$ at 24 h | 0.382                   | 0.057              | 0.039                            |

Under the prespecified direction-only rule, all three directions persisted. Two qualifications are important. The HIRISA reversal was much weaker on the broader set: its 95% interval (−2.35 to 0.37) included zero, and donor 3283BW, the only donor without elevated control monocyte expression, changed sign (−0.96 to +0.87); the four donors with elevated controls stayed negative (−1.01 to −2.00). The strongly negative contrast therefore appears concentrated in the ten-gene panel and in donors with high baseline expression, suggesting a possible contribution from baseline expression that remains causally unresolved. In GSE226572, the rise of the contrast with exposure time generalized (Γ positive in every donor on both sets), and the positive contrast at 24 h varied across specifications. Magnitudes are not compared across gene sets or cohorts. Outputs are in `results/hallmark_gene_set/`.

##### Comparison of the three cohorts

The three cohorts differ in every design factor at once, so no cross-cohort difference can be assigned to one factor. No donor contrasts were pooled, and the cohorts were not treated as exchangeable.

|                                                                     | GSE96583                               | HIRISA (GSE306664)                                      | GSE226572                                  |
|---------------------------------------------------------------------|----------------------------------------|---------------------------------------------------------|--------------------------------------------|
| Donors; health                                                      | 8; SLE                                 | 5; healthy                                              | 3; healthy (all male)                      |
| Preparation                                                         | Mixed PBMC; donors pooled by condition | Monocytes and NK cells enriched and cultured separately | Mixed PBMC; one library per donor and time |
| IFN-β dose; exposure                                                | 100 U/mL; 6 h                          | 100 U/mL; 21 h                                          | 1,000 U/mL; 1-36 h                         |
| Assay                                                               | 10x droplet 3'                         | 10x Flex                                                | 10x 3' v3                                  |
| Cell labels                                                         | Depositor                              | Depositor                                               | Marker rule; depositor labels as a check  |
| $`\hat\theta`$ (prior count 2)                                      | 1.834                                  | −4.609                                                  | −0.173 (4 h); 0.382 (24 h)                 |
| Change under $`\log_2(\mathrm{CPM}+1)`$                             | −0.601                                 | +0.155                                                  | −0.009                                     |
| Mean control monocyte panel expression, TMM log-CPM (prior count 2) | 3.58-6.66                              | 5.48 (one donor); 9.32-11.45 (four donors)              | 5.71-7.27                                  |

Control expression was placed on one scale for all cohorts in `results/cross_study_synthesis/baseline_harmonization.csv`; cross-assay comparisons of absolute levels remain approximate. Four points follow:

1) The same estimand, computed identically, took clearly different values in three designs and changed sign with time within one of them, supporting context dependence with the responsible factors remaining unresolved.
2) The rise in GSE226572 challenges a simple duration-only explanation for the HIRISA reversal. Interactions between stimulation duration, culture conditions and other factors remain possible. The causal contributions of enrichment, baseline activation, disease status, dose and assay remain unresolved.
3) GSE226572 controls cultured for 36 h in mixed PBMC stayed near the GSE96583 range, indicating that prolonged culture can occur with lower baseline panel expression.
4) The effect of low-count regularization followed library size in all three cohorts (0.60 with libraries of 0.1-0.8 million counts, at most 0.16 with 44-75 million, at most 0.01 with 4-7 million), which supports the library-size account and adds no new empirical finding. Establishing the cause of the reversal and of the elevated HIRISA control expression would require an experiment that varies preparation within donors at a fixed dose and duration, with IFN measured in the culture.

#### Conclusions

The results support five conclusions of differing strength:

1) In GSE96583, monocytes exhibited greater mean induction than NK cells in every analysis, in every subsample and with both larger IFN gene sets.
2) The magnitude of the contrast ranged from 1.23 to 1.87 in analyses using all cells and was approximately 1.45-1.47 in subsamples. It depended on the regularization of low counts (the prior count or a fixed CPM offset), the number of cells aggregated and the gene set. Any single estimate should therefore be reported together with these analytical choices.
3) The contrast was positive in seven donors across all analyses. Its sign for donor 1016 depended on the regularization of low counts and on cell selection; this donor also had the highest control monocyte expression.
4) In an independent cohort of healthy donors whose monocytes and NK cells were cultured separately for 21 h, the contrast was reversed (−4.609). The reversal coincided with high panel expression in untreated monocytes of four of five donors, whose cause could not be determined; the depositors' independent analysis of the same cells by MAST agreed in direction for all ten genes. The regularization of low counts changed this contrast by only 0.155, consistent with the much larger and more similar library sizes. The direction of the GSE96583 contrast should therefore not be generalized to other experimental designs.
5) In GSE226572, the monocyte–NK contrast increased between 4 h and 24 h in all three donors under four annotation definitions (Γ = 0.555, 0.595, 0.597 and 0.553). Independent annotations support this sensitivity result within the same three donors (`results/time_course/annotation_validation/gamma_by_annotation.csv`). The original publication examined IFN-response dynamics ([de Cevins et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10772457/)); the present contribution concerns the fixed donor-level contrast and its analytical robustness.

> [!TIP]
> These findings suggest four reporting practices for studies that compare IFN scores between cell types:
>
> - Report the direction, the magnitude and the donor-level results separately.
> - Report the library size of each cell type, and the CPM added by the prior count or offset, with any log-CPM comparison, because an identical prior count acts differently on cell types with different library sizes, and an identical CPM offset corresponds to different numbers of counts.
> - Report control expression together with induction, because higher control expression can reduce estimated induction; baseline signalling activity requires direct measurement.
> - Treat scores from a larger gene set as a distinct measurement; their scale depends on gene-set composition, which also influenced between-donor variability in this analysis.

#### Limitations

In the original experiment, cells from all donors were pooled by condition, so that all control cells and all stimulated cells were each processed as a single pool. Effects of stimulation therefore cannot be fully separated from differences in the processing of the two pools. This confounding persists after donor pairing and subsampling.

Secondly, the primary analysis comprises eight SLE donors at one IFN-β dose and time point and two broadly defined cell types. Clinical information was not available. The independent cohort examined (HIRISA) differs from GSE96583 in donor health, stimulation duration, cell preparation and platform, comprises five donors, and therefore assesses transferability across experimental contexts. GSE226572 adds only three healthy donors in a third experimental setting; its annotation definitions reuse those donors. Its second untreated library per donor has an undocumented role, and its smallest attainable two-sided donor sign-flip p-value is 0.25. The cohorts are not exchangeable and are not pooled.

Furthermore, the control NK group of donor 1039 contained only 18 cells. All 13 mitochondrial genes have zero counts in the supplied matrix, which precluded filtering on mitochondrial content. The deposited annotation lists 24,679 singlets; the original article reports 24,305; the source of this discrepancy could not be determined.

Finally, all results concern normalized transcript levels. Protein, signalling and functional measurements were not available.

A more definitive test would specify the gene set and all analytical parameters and apply them to an independent donor cohort, with control and stimulated cells distributed across multiple processing batches. HIRISA met the batch condition. Its untreated monocytes, enriched and cultured alone for 21 h, expressed the panel highly in four of five donors, with the cause unresolved; a replication would also require mixed-cell culture and a stimulation duration comparable to that of GSE96583. A time course with matched signalling measurements would be required to establish whether the transcriptional difference reflects a difference in cellular responsiveness to IFN.

#### Figures

##### Donor-level robustness

![Donor robustness](results/primary_contrast/figures/collection1_donor_robustness.png)

**(a) Donor contrasts.** Each point represents the monocyte-minus-NK induction for one donor; the diamond and bar indicate the mean and its 95% interval. All eight contrasts exceed zero, which constitutes the principal result. Donor 1016 (0.23) is distinct from the remaining seven donors (1.42-2.91) and is the donor whose sign subsequently varies. **(b) Mean contrast after exclusion of each donor.** Circles denote exclusion of the donor contrast; triangles denote re-running of the complete pipeline without the donor. The mean remains between 1.68 and 2.10 (dashed line, all donors), indicating that no single donor determines the result. Agreement between the two procedures indicates that re-estimating normalization without a donor has little effect. **(c) Mean contrast across 100 subsamples of 18 cells per group.** All draws exceed zero at both prior counts. The full-data estimate (dashed line) lies at the upper extreme of the prior-count-2 distribution, indicating that the number of aggregated cells affects the magnitude of the estimate. The distribution at prior count 0.5 is wider. **(d) Donor contrasts across analyses.** The ordering of donors is broadly preserved across ten-gene analyses. Only donor 1016 falls below zero, and only at prior count 0.5 (−0.034; −0.043 with TMM factors estimated without IFN genes). The Hallmark gene sets (squares) reduce the dispersion between donors.

##### Sensitivity to analytical choices

![Specification sensitivity](results/primary_contrast/figures/collection2_specification_sensitivity.png)

**(e) Mean contrast across analyses.** All estimates and intervals exceed zero, with point estimates between 1.40 and 1.87; this panel summarizes the robustness of the direction. Solid bars denote 95% intervals across donors, and dashed bars denote the range across cell subsamples; these quantify different sources of variability. **(f) Quality control and prior count.** Exclusion of high-UMI cells has minimal effect on the mean (1.870); reduction of the prior count produces the largest change. **(g) Expression-matched random panels.** The histogram shows contrasts for 1,000 random non-IFN panels; the blue line indicates the IFN panel. The IFN panel exceeds all random panels, indicating a substantially larger difference than for non-IFN genes of comparable control expression. The random panels are centred slightly above zero (0.233), indicating a small general monocyte-greater difference on this scale. **(h) Induction in each cell type.** Reducing the prior count increases induction in both cell types, more so in NK cells, which accounts for the reduced contrast. Subsampling reduces induction in both cell types, by 1.2-1.6 $`\log_2`$ units at prior count 2 and by 0.5-0.7 at prior count 0.5. Induction within each cell type is therefore less stable than the difference between cell types.

##### Expression and normalization

![Expression and normalization](results/primary_contrast/figures/collection3_expression_normalization.png)

**(i) Control and stimulated expression of each panel gene.** Circles indicate control expression and arrowheads stimulated expression; arrow length corresponds to induction. All genes are induced in both cell types, to differing extents. IFI6 is the exception, with higher control expression and slightly lower induction in monocytes. Reducing the prior count principally alters the control values of sparsely detected genes such as IFI44L and IFIT1. Because the two cell types were normalized separately, horizontal positions are not directly comparable between them. **(j) Decomposition of each donor contrast.** For each donor, the contrast (×) equals the monocyte-minus-NK difference after stimulation (▲) plus the NK-minus-monocyte difference in controls (●). Because the cell types were normalized separately, the partition between these terms is arithmetic and lacks biological interpretation. The comparison across donors is informative: donor 1016 alone has a negative control term, consistent with its near-zero contrast. **(k) Gene contributions to the effect of the prior count.** Using all cells, the reduction in the contrast at prior count 0.5 arises principally from increased NK induction of RSAD2, IFIT1 and IFI44L. In subsamples, the monocyte and NK changes approximately cancel (net change, −0.023). The effect of the prior count therefore depends on the number of aggregated cells. **(l) Induction of all genes used for normalization.** Non-IFN genes are centred near zero in every donor and cell type, as expected under appropriate normalization. Hallmark IFN genes are induced more strongly in monocytes than in NK cells in every donor, indicating that the pattern extends beyond the ten panel genes. Panel genes are among the most strongly induced genes in every sample.

##### Fixed-offset expression measurement

![Fixed-offset comparison](results/primary_contrast/figures/14_fixed_offset_donor_contrasts.png)

**(m) Donor contrasts under edgeR log-CPM (prior count 2) and $`\log_2(\mathrm{CPM}+1)`$.** Circles show the primary contrasts and triangles the fixed-offset contrasts; grey segments connect the two values for each donor, and the bars at the top show each mean with its 95% interval. Every contrast decreases under the fixed offset, and only donor 1016 falls below zero (−0.142). Both measurements use identical counts, TMM factors, genes and donors; they differ only in the CPM added before log transformation.

##### Independent cohort (HIRISA)

![HIRISA donor contrasts](results/independent_cohort/figures/h2_donor_contrasts.png)

**(n) HIRISA donor contrasts under the three measurements.** Filled circles show edgeR log-CPM at prior count 2, open circles prior count 0.5 (largely hidden, as the values nearly coincide), and triangles $`\log_2(\mathrm{CPM}+1)`$; dashed lines show the means. All five contrasts are negative. Donor 3283BW, the only donor with low control monocyte expression, has the contrast closest to zero (−0.963). The three measurements differ by at most 0.29 in any donor.

#### Data and code availability

The input data is publicly available from [GEO GSE96583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583) and was accessed through `muscData::Kang18_8vs8()` (ExperimentHub EH2259, muscData 1.26.0). The local count cache has the MD5 checksum `b890d4b6a805f6f7be00e89927d4b34e`. The independent cohort comprises 20 per-sample count files from [GEO GSE306664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306664), listed with their MD5 checksums in `analysis/independent_cohort/samples.tsv` and `analysis/independent_cohort/downloads.tsv`, and read with rhdf5 2.56.1. GSE226572 provides 24 count files and a deposited Seurat object for the three-donor time course; URLs, sizes and checksums are recorded in `analysis/time_course/samples.tsv`, `analysis/time_course/downloads.tsv` and `analysis/time_course/annotation_validation/downloads.tsv`. No new experimental data were generated.

This project is implemented in R 4.6.2 with Bioconductor 3.23. The following packages were used:

- edgeR 4.10.7,
- SingleCellExperiment 1.34.0,
- Matrix 1.7-6, and
- ggplot2 4.0.3.

The CRAN packages `msigdbr` and `httr2` were also used to retrieve gene annotations.

Verify the archived release first:

```bash
git clone --branch v1.2.0 https://github.com/zeyuyaoy/interferon-response-robustness.git
cd interferon-response-robustness
Rscript release/verify_release.R
```

For a full rerun, use a separate checkout and install R 4.6.2, Bioconductor 3.23 and the packages recorded in `release/renv.lock`. With `renv` installed, restore with `Rscript -e 'renv::restore(lockfile="release/renv.lock", prompt=FALSE); renv::activate()'`. The activation creates a local `.Rprofile` so subsequent `Rscript` commands use the restored project library. The saved `sessionInfo.txt` files describe each original run, and `release/software_versions.tsv` records the finalization environment. Restoration requires network access and appropriate platform libraries (including HDF5); R and operating-system libraries require separate installation.

Run the commands below from the repository root. The dataset is loaded from `cache/Kang18_8vs8.Rda` (MD5 above; not tracked by Git) if present, and otherwise downloaded through `muscData::Kang18_8vs8()`. The analysis specification comprises `analysis/primary_contrast/scripts/config.R` and the files in `analysis/inputs/` (listed in `analysis/inputs/manifest.txt`). Each output table is labelled with the identifier of the run that produced it; the identifier ends in `-dirty` if the analysis sources differed from the recorded commit.

```bash
Rscript analysis/primary_contrast/tests/test_utils.R
Rscript analysis/primary_contrast/scripts/run_all.R
```

> [!IMPORTANT]
> Each run clears `results/primary_contrast/` and `results/primary_contrast/figures/`; existing outputs should therefore be preserved by running the analysis in a separate working copy. Other experiment scripts also overwrite their output files. Reruns change timestamps and run identifiers, so the archived checksum/run-ID verifier applies to the untouched release. Compare regenerated numerical tables with the archived tables and inspect each new run’s checks and status records.

The independent-cohort analysis is separate from the pipeline above and writes to `results/independent_cohort/` and `results/independent_cohort/figures/`. It downloads about 727 MB to `cache/hirisa/` (not tracked by Git) and requires about 5 GB of memory:

```bash
Rscript analysis/independent_cohort/scripts/00_download.R
Rscript analysis/independent_cohort/scripts/01_pilot.R
Rscript analysis/independent_cohort/scripts/02_analysis.R
```

The comparison with the published HIRISA differential-expression results downloads two tables to `cache/hirisa_deg/` and writes to `results/independent_cohort/triangulation/`. The exploratory time-course analysis downloads 24 count files (about 1.2 GB) from [GEO GSE226572](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226572) to `cache/gse226572/` and writes to `results/time_course/`. Its annotation check also downloads the depositors' Seurat object (10.4 GB) to the same folder and reads only its cell metadata:

```bash
Rscript analysis/independent_cohort/triangulation/01_triangulation.R
Rscript analysis/time_course/scripts/00_manifest.R
Rscript analysis/time_course/scripts/01_download.R
Rscript analysis/time_course/scripts/02_annotation.R
Rscript analysis/time_course/tests/test_gse226572.R
Rscript analysis/time_course/scripts/03_analysis.R
Rscript analysis/time_course/annotation_validation/00_download_author_object.R
Rscript analysis/time_course/annotation_validation/01_extract_author_metadata.R
Rscript analysis/time_course/annotation_validation/02_map_author_labels.R
Rscript analysis/time_course/annotation_validation/test_rds_stream.R
Rscript analysis/time_course/annotation_validation/03_annotation_validation.R
Rscript analysis/hallmark_gene_set/scripts/00_eligible_genes.R
Rscript analysis/hallmark_gene_set/scripts/01_hallmark.R
Rscript analysis/cross_study_synthesis/baseline_harmonization.R
```

Raw inputs and derived cell caches belong under Git-ignored `cache/`, including the 10,435,582,500-byte Seurat object. They are not included in the release. The frozen gene sets and annotation snapshot are versioned; do not refresh annotations to reproduce this study. Download scripts validate existing manifests. To reproduce the primary historical table comparison after a full run in a separate checkout:

```bash
Rscript release/compare_runs.R /path/to/archived/results/primary_contrast results/primary_contrast /tmp/reproduction_comparison.csv
```

The comparison ignores run IDs and reports numerical and text differences; inspect its CSV, since it does not exit nonzero on a mismatch. Historical paths in locked specifications and original run logs are intentional; see [the path map](release/FINALIZATION.md). Automated checks do not resolve the study’s small donor counts, design confounding or unmeasured mechanisms.

#### Final research artifacts

The scientific investigation is closed after author-annotation validation. Release v1.2.0 preserves the [concise report](release/research_report.md), [seven-slide presentation](release/presentation.pptx), [release notes](release/RELEASE_NOTES.md), and [claim audit and provenance record](release/FINALIZATION.md). The report states the study’s contribution and limitations, including the original publication’s prior analysis of IFN-response dynamics.

#### References

BibTeX entries for all references are in `references/references.bib`.

#### Layout

Each experiment has a folder under `analysis/` (code and locked inputs) and a folder of the same name under `results/` (outputs and figures).

```
README.md
references/        references.bib
release/           checksums.sha256, verification and reproduction records
analysis/
  inputs/          frozen inputs shared by all experiments (gene sets, annotation snapshot)
  shared/          utils_*.R used by all experiments
  primary_contrast/        GSE96583, eight SLE donors
    scripts/       config.R, run_all.R, 00_fetch_annotation.R ... 14_fixed_offset_measurement.R
    tests/         test_utils.R
  independent_cohort/      HIRISA, GSE306664
    specification.tsv, spec_history.tsv, samples.tsv, downloads.tsv
    scripts/       00_download.R, 01_pilot.R, 02_analysis.R, utils_hirisa.R
    triangulation/ comparison with published HIRISA MAST results
  time_course/             GSE226572, exploratory
    specification.tsv, spec_history.tsv, samples.tsv, downloads.tsv
    scripts/, tests/
    annotation_validation/ author cell labels from deposited Seurat object
  hallmark_gene_set/       one fixed Hallmark IFN-alpha set across three cohorts
    specification.tsv, spec_history.tsv, eligible_genes.tsv, scripts/
  cross_study_synthesis/   baseline_harmonization.R
results/
  primary_contrast/        <step>/ ..., figures/
  independent_cohort/      tables, figures/, triangulation/ (tables, figures/)
  time_course/             tables, figures/, annotation_validation/
  hallmark_gene_set/
  cross_study_synthesis/
```
