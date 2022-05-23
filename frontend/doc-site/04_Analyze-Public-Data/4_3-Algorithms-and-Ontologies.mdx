# Algorithms and Ontologies


## Algorithms used by Cellxgene

Differential expression

We're . The  assumes normally distributed values on a linear scale.

Cellxgene's differential expression uses a , which assumes that the two populations are normally distributed, but may have unequal variance. We use a two-sided t-test against the null hypothesis that the two populations have equal means. P-values are adjusted with the .

To help avoid spurious results, we only return genes with log fold change greater than 0.01:

|log2( mean(set1) / mean(set2) )| > 0.01

The log fold change threshold can be configured on cellxgene Desktop with the  command. We then sort genes by their associated |t value| and return the top 15 genes.


## Cell type ontology ordering

Cell types in a dot plot can be ordered using lineage information to maximize the utility of scExpression for pattern discovery and readability. The cell ontology (CL) directed acyclic graph (DAG) provides a good approximation to biological lineages. Thus we have a developed a method to create an ordered list of cell types based on the information contained in the CL DAG.

The method traverses the DAG starting from the root until reaching the leftmost bottom node adding nodes present in the data along the way. Once it reaches the bottom, it recursively repeats the process for any siblings. After, it traverses back and repeats the sibling process until exhaustion.


## Algorithm

Load the CL ontology, using either the cl.obo or cl.owl from the source file .

Subset the ontology to only contain the cell types of interest along with their ancestors.

Create a 2-dimensional representation of the DAG using the dot method of graphiz. See  for a full description of this method.

Traverse the tree starting at the root:

If node is present in cell types of interest, add to the final ordered list.

Remove the node from nodes of interest.

Get children nodes

Perform recursion for each child node, going left-to-right through the children nodes.

Return the ordered list


## Example

The following is an example of how a subset of the CL DAG looks like when rendering it with the dot algorithm from graphiz. Orange nodes represent cell types that would be present in the data corpus and are the target labels to be ordered.

The cell type ordering algorithm described above would return the following ordered list of cell types.

Native Cell

Secretory Cell

Club Cell

Chondrocyte

Goblet Cell

Mast Cell

T Cell

CD4 positive, Alpha beta T Cell

CD8 positive, Alpha beta T Cell

B Cell

Natural Killer

Dendritic Cell

Monocyte

Macrophage

Alveolar Macrophage

Basal Cell

Fibroblast

Megakaryocyte

Enucleate Erythrocyte

Myofibroblast Cell

Mesothelial Cell

Endothelial Cell

Blood Vessel Endothelial Cell

Smooth Muscle Cell

Vascular Associated Smooth Muscle Cell

Epithelial Cell

Squamous Epithelial Cell

Ciliated Cell
