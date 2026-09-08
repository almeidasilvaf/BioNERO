# BioNERO: Biological Network Reconstruction Omnibus

BioNERO aims to integrate all aspects of biological network inference in
a single package, including data preprocessing, exploratory analyses,
network inference, and analyses for biological interpretations. BioNERO
can be used to infer gene coexpression networks (GCNs) and gene
regulatory networks (GRNs) from gene expression data. Additionally, it
can be used to explore topological properties of protein-protein
interaction (PPI) networks. GCN inference relies on the popular WGCNA
algorithm. GRN inference is based on the "wisdom of the crowds"
principle, which consists in inferring GRNs with multiple algorithms
(here, CLR, GENIE3 and ARACNE) and calculating the average rank for each
interaction pair. As all steps of network analyses are included in this
package, BioNERO makes users avoid having to learn the syntaxes of
several packages and how to communicate between them. Finally, users can
also identify consensus modules across independent expression sets and
calculate intra and interspecies module preservation statistics between
different networks.

## See also

Useful links:

- <https://github.com/almeidasilvaf/BioNERO>

- Report bugs at <https://github.com/almeidasilvaf/BioNERO/issues>

## Author

**Maintainer**: Fabricio Almeida-Silva
<fabricio_almeidasilva@hotmail.com>
([ORCID](https://orcid.org/0000-0002-5314-2964))

Authors:

- Thiago Venancio ([ORCID](https://orcid.org/0000-0002-2215-8082))
