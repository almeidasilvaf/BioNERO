# Package index

## All functions

- [`PC_correction()`](PC_correction.md) : Apply Principal Component
  (PC)-based correction for confounding artifacts
- [`SFT_fit()`](SFT_fit.md) : Pick power to fit network to a scale-free
  topology
- [`ZKfiltering()`](ZKfiltering.md) : Filter outlying samples based on
  the standardized connectivity (Zk) method
- [`check_SFT()`](check_SFT.md) : Check scale-free topology fit for a
  given network
- [`consensus_SFT_fit()`](consensus_SFT_fit.md) : Pick power to fit
  networks to scale-free topology
- [`consensus_modules()`](consensus_modules.md) : Identify consensus
  modules across independent data sets
- [`consensus_trait_cor()`](consensus_trait_cor.md) : Correlate
  set-specific modules and consensus modules to sample information
- [`cor2adj()`](cor2adj.md) : Calculate an adjacency matrix from a
  correlation matrix
- [`cormat_to_edgelist()`](cormat_to_edgelist.md) : Transform a
  correlation matrix to an edge list
- [`detect_communities()`](detect_communities.md) : Detect communities
  in a network
- [`dfs2one()`](dfs2one.md) : Combine multiple expression tables (.tsv)
  into a single data frame
- [`enrichment_analysis()`](enrichment_analysis.md) : Perform
  overrepresentation analysis for a set of genes
- [`exp2cor()`](exp2cor.md) : Calculate pairwise correlations between
  genes in a matrix
- [`exp2gcn()`](exp2gcn.md) : Infer gene coexpression network from gene
  expression
- [`exp2gcn_blockwise()`](exp2gcn_blockwise.md) : Infer gene
  coexpression network from gene expression in a blockwise manner
- [`exp2grn()`](exp2grn.md) : Infer gene regulatory network from
  expression data
- [`exp_genes2orthogroups()`](exp_genes2orthogroups.md) : Collapse
  gene-level expression data to orthogroup level
- [`exp_preprocess()`](exp_preprocess.md) : Preprocess expression data
  for network reconstruction
- [`filt.se`](filt.se.md) : Filtered maize gene expression data from
  Shin et al., 2021.
- [`filter_by_variance()`](filter_by_variance.md) : Keep only genes with
  the highest variances
- [`gene_significance()`](gene_significance.md) : Calculate gene
  significance for a given group of genes
- [`get_HK()`](get_HK.md) : Get housekeeping genes from global
  expression profile
- [`get_edge_list()`](get_edge_list.md) : Get edge list from an
  adjacency matrix for a group of genes
- [`get_hubs_gcn()`](get_hubs_gcn.md) : Get GCN hubs
- [`get_hubs_grn()`](get_hubs_grn.md)
  [`get_hubs_ppi()`](get_hubs_grn.md) : Get hubs for gene regulatory
  network
- [`get_neighbors()`](get_neighbors.md) : Get 1st-order neighbors of a
  given gene or group of genes
- [`grn_average_rank()`](grn_average_rank.md) : Rank edge weights for
  GRNs and calculate average across different methods
- [`grn_combined()`](grn_combined.md) : Infer gene regulatory network
  with multiple algorithms and combine results in a list
- [`grn_filter()`](grn_filter.md) : Filter a gene regulatory network
  based on optimal scale-free topology fit
- [`grn_infer()`](grn_infer.md) : Infer gene regulatory network with one
  of three algorithms
- [`is_singleton()`](is_singleton.md) : Logical expression to check if
  gene or gene set is singleton or not
- [`modPres_WGCNA()`](modPres_WGCNA.md) : Calculate module preservation
  between two expression data sets using WGCNA's algorithm
- [`modPres_netrep()`](modPres_netrep.md) : Calculate module
  preservation between two expression data sets using NetRep's algorithm
- [`module_enrichment()`](module_enrichment.md) : Perform enrichment
  analysis for coexpression network modules
- [`module_preservation()`](module_preservation.md) : Calculate network
  preservation between two expression data sets
- [`module_stability()`](module_stability.md) : Perform module stability
  analysis
- [`module_trait_cor()`](module_trait_cor.md) : Correlate module
  eigengenes to trait
- [`net_stats()`](net_stats.md) : Calculate network statistics
- [`og.zma.osa`](og.zma.osa.md) : Orthogroups between maize and rice
- [`osa.se`](osa.se.md) : Rice gene expression data from Shin et al.,
  2021.
- [`parse_orthofinder()`](parse_orthofinder.md) : Parse orthogroups
  identified by OrthoFinder
- [`plot_PCA()`](plot_PCA.md) : Plot Principal Component Analysis (PCA)
  of samples
- [`plot_dendro_and_colors()`](plot_dendro_and_colors.md) : Plot
  dendrogram of genes and modules
- [`plot_eigengene_network()`](plot_eigengene_network.md) : Plot
  eigengene network
- [`plot_expression_profile()`](plot_expression_profile.md) : Plot
  expression profile of given genes across samples
- [`plot_gcn()`](plot_gcn.md) : Plot gene coexpression network from edge
  list
- [`plot_gene_significance()`](plot_gene_significance.md) : Plot a
  heatmap of gene significance
- [`plot_grn()`](plot_grn.md) : Plot gene regulatory network from edge
  list
- [`plot_heatmap()`](plot_heatmap.md) : Plot heatmap of hierarchically
  clustered sample correlations or gene expression
- [`plot_module_trait_cor()`](plot_module_trait_cor.md) : Plot a heatmap
  of module-trait correlations
- [`plot_ngenes_per_module()`](plot_ngenes_per_module.md) : Plot number
  of genes per module
- [`plot_ppi()`](plot_ppi.md) : Plot protein-protein interaction network
  from edge list
- [`q_normalize()`](q_normalize.md) : Quantile normalize the expression
  data
- [`remove_nonexp()`](remove_nonexp.md) : Remove genes that are not
  expressed based on a user-defined threshold
- [`replace_na()`](replace_na.md) : Remove missing values in a gene
  expression data frame
- [`zma.interpro`](zma.interpro.md) : Maize Interpro annotation
- [`zma.se`](zma.se.md) : Maize gene expression data from Shin et al.,
  2021.
- [`zma.tfs`](zma.tfs.md) : Maize transcription factors
