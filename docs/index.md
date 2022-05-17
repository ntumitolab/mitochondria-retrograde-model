# Mitochondria retrograde ODE model

![GitHub repo size](https://img.shields.io/github/repo-size/ntumitolab/mitochondria-retrograde-model)

> Codes, data and notebooks for the manuscript titled "Mathematical Modeling and Analysis of Mitochondrial Retrograde Signaling Dynamics"

## Graphical abstract

![Graphical abstract](https://user-images.githubusercontent.com/29009898/130342513-081f4592-3cc6-4468-ba3b-868416f3be6b.png)

## Mitochondrial Retrograde Signalling

Mitochondrial retrograde signaling reports mitochondrial status to the nucleus. However, there is a lack of understanding of how the nucleus capture mitochondrial status in dynamics and information processing. It is a complicated biochemical reaction that occurs in most eukaryotic organisms. In this repository, we focus on the RTG pathway in yeast. This pathway is the simplest retrograde signaling pathway that has been investigated thoroughly. Data are collected from [1] and [2] (See [Dataset](#dataset)). This repository aims to compose known protein interactions and nucleus relocation that fulfills all known responses of the yeast RTG pathway. Monte-Carlo approach is used to solve this Boolean satisfiability problem, and the parameter searching/ simulation/ threading is facilitated by DifferentialEquations.jl [3].

## Installation

Julia v1.7+ is needed (https://julialang.org/).

## RTG protein levels

<https://github.com/ntumitolab/rnaseq_rtg_expression>

## Parameters of retrograde signaling model

> Table of parameters: [`solution_rtgM4.csv`](https://github.com/ntumitolab/mitochondria-retrograde-model/blob/main/src/data/solution_rtgM4.csv)

This file is generated and modified from [a script](https://github.com/stevengogogo/RetroSignalModel.jl/blob/82cc8fd91b1fa03983c08a5eed1144bac3e93dcc/scripts/find_valid_solutions.jl). All solutions are corresponding to the knockout experiments of [1] and [2] with the conditions in [boolean_table_RTG13.csv](https://github.com/ntumitolab/mitochondria-retrograde-model/blob/main/src/data/boolean_table_RTG13.csv). The solutions are stored in [`solution_rtgM4.csv`](https://github.com/ntumitolab/mitochondria-retrograde-model/blob/main/src/data/solution_rtgM4.csv).

## Responses of Yeast RTG proteins to mitochondrial damage

> Data: [boolean_table_RTG13.csv](https://github.com/ntumitolab/mitochondria-retrograde-model/blob/main/src/data/boolean_table_RTG13.csv)

This folder contains summarized responses of mitochondrial retrograde signaling in yeast.

### Components

| Standard Name | Variable Name | Details                                      |
| ------------- | ------------- | -------------------------------------------- |
| RTG1          | `rtg1`        | https://www.yeastgenome.org/locus/S000005428 |
| RTG2          | `rtg2`        | https://www.yeastgenome.org/locus/S000005428 |
| RTG3          | `rtg3`        | https://www.yeastgenome.org/locus/S000000199 |
| Mks1          | `mks1`        | https://www.yeastgenome.org/locus/S000005020 |

### Definition of Response

In [1] and [2], RTG response is observed via GFP tags on either RTG1 or RTG3. In wild-type, mitochondrial damage can cause these proteins to accumulate in the nucleus, resulting in the intensified brightness of the nucleus region observed by fluorescent microscopy. As shown in [boolean_table_RTG13.csv](src/data/boolean_table_RTG13.csv), the responses are categorized in binary results: whether GFP is accumulated in the nucleus in a given condition. Based on [1] and [2], there are 20 reactions listed in the table.

For example, the following is one of the conditions mentioned in [1]:

| Rtg1 | Rtg2 | Rtg3 | s   | Mks | gfp  | Trans2Nuc |
| ---- | ---- | ---- | --- | --- | ---- | --------- |
| 0    | 0    | 1    | 1   | 1   | rtg3 | 1         |

Under the columns of `Rtg1`, `Rtg2`, `Rtg3` and `Mks`, `0` means that the given protein is suppressed by knockout. On the other hand, `1` represent an expression of wild type. Also, `1` in `s` represent mitochondrial dysfunction, and `0` means the absence of mitochondrial damage. The `gfp` column describes the location of GFP tag. In this example, GFP tag is on `Rtg3`. As known in [1], `Rtg3-GFP` translocates to the nucleus under this condition. Therefore, `Trans2Nuc` is marked as `1`, which means the GFP tags nucleus translocation happens.

### Reactions

There are 20 reactions summarized in the table. Some conditions are yet to be explored; some are from [1] (Sekito et al. 2000) or [2] (Sekito et al. 2002). Missing conditions are labeled with `NA`.

| Line Number | Reference |
| ----------- | --------- |
| 2           | NA        |
| 3           | [1]       |
| 4           | [1]       |
| 5           | [1]       |
| 6           | [1]       |
| 7           | [1]       |
| 8           | [1]       |
| 9           | [1]       |
| 10          | NA        |
| 11          | NA        |
| 12          | [1]       |
| 13          | [1]       |
| 14          | [1]       |
| 15          | [1]       |
| 16          | [1]       |
| 17          | [1]       |
| 18          | [2]       |
| 19          | [2]       |
| 20          | [2]       |
| 21          | [2]       |

## Notebook execution status

```{nb-exec-table}
```

## References

1. Sekito, Takayuki, Janet Thornton, and Ronald A. Butow. "Mitochondria-to-nuclear signaling is regulated by the subcellular localization of the transcription factors Rtg1p and Rtg3p." Molecular biology of the cell 11.6 (2000): 2103-2115. URL: https://doi.org/10.1091/mbc.11.6.2103
2. Sekito, Takayuki, Zhengchang Liu, Janet Thornton, and Ronald A. Butow. “RTG-Dependent Mitochondria-to-Nucleus Signaling Is Regulated by MKS1 and Is Linked to Formation of Yeast Prion [URE3].” Molecular Biology of the Cell 13, no. 3 (March 2002): 795–804. https://doi.org/10.1091/mbc.01-09-0473.
3. Rackauckas, Christopher, and Qing Nie. “DifferentialEquations.Jl – A Performant and Feature-Rich Ecosystem for Solving Differential Equations in Julia.” Journal of Open Research Software 5, no. 1 (May 25, 2017): 15. https://doi.org/10.5334/jors.151.
4. Gasch, Audrey P., et al. "Single-cell RNA sequencing reveals intrinsic and extrinsic regulatory heterogeneity in yeast responding to stress." PLoS biology 15.12 (2017): e2004050. URL: https://doi.org/10.1371/journal.pbio.2004050
5. Delmans, Mihails, and Martin Hemberg. "Discrete distributional differential expression (D3E)-a tool for gene expression analysis of single-cell RNA-seq data." BMC bioinformatics 17.1 (2016): 1-13. URL: https://doi.org/10.1186/s12859-016-0944-6
