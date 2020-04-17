# M3S: A user-friendly tool for comprehensive evaluation of multi-modality and statistical model selection for single-cell RNA sequencing data.

A large number of single-cell RNA sequencing (scRNA-seq) data set was
recently generated to characterize the heterogeneous in cell types or
states in a complex tissue or biological process. Noting that the gene
expression in a single cell is purely determined by the transcriptional
regulatory signals, which are largely varied through different cells,
single gene's expression profile naturally form into a multi-modal
distribution, with each component corresponds to a possible regulatory
state. Multiple statistical models have been developed to capture such
multi-modalities (MM), for the data generated of different conditions or
generated by different experimental platforms. Such models include
Quasi-Poisson (QP), Negative Binomial (NB), Zero Inflated Negative
Binomial (ZINB), Zero Inflated Gaussian (ZIG), Mixture Gaussian (MG),
Beta Poisson (BP), Zero Inflated Mixture Gaussian (ZIMG), and Left
Truncated Mixture Gaussian (LTMG), majorly differ by their assumptions
of "drop-out", errors and MM. We have recently developed a system
biological model to decompose possible contributors of the
multi-modality and errors of scRNA-seq data. Our analysis and other
recent works clearly suggested varied sources of the errors and MM that
are highly affected by experimental condition and platform. However,
there is lack of a computational capability to select the most proper
statistical model for each scRNA-seq data.

Motivated by this idea, we developed a user-friendly R package,
M<sup>3</sup>S, to (1) select the most proper statistical models and
differential gene expression analysis method, (2) characterize varied
transcriptional regulatory state, and (3) identify differentially varied
gene expressions, for a given scRNA-seq. In addition, the package can be
also applied to bulk tissue transcriptomics or other omics data for
characterization of the MM of features.

#### [Download Paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3243-1)

Zhang, Yu, Changlin Wan, Pengcheng Wang, Wennan Chang, Yan Huo, Jian Chen, Qin Ma, Sha Cao, and Chi Zhang. "M3S: A comprehensive model selection for multi-modal single-cell RNA sequencing data." BMC bioinformatics 20, no. 24 (2019): 1-5.

## Citations
If you find the code helpful in your resarch or work, please cite the following papers.
```BibTex
@article{zhang2019m3s,
  title={M3S: A comprehensive model selection for multi-modal single-cell RNA sequencing data},
  author={Zhang, Yu and Wan, Changlin and Wang, Pengcheng and Chang, Wennan and Huo, Yan and Chen, Jian and Ma, Qin and Cao, Sha and Zhang, Chi},
  journal={BMC bioinformatics},
  volume={20},
  number={24},
  pages={1--5},
  year={2019},
  publisher={BioMed Central}
}
```
