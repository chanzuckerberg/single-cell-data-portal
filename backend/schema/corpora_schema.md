# Corpora Schema

Authors: mkinsella@chanzuckerberg.com

Document Status: _Approved_

Version: 1.0.0

Date Last Modified: 2020-07-21


## Background

Corpora is for "publishing, exploring, integrating, and re-using" single-cell datasets. There is a tension in those goals. On one hand, focusing on
publication and exploration means Corpora should let submitters be in full control of how their datasets appear on the platform. This includes everything 
from metadata fields to the colors in their plots. On the other hand, focusing on integration and re-use means Corpora should require datasets to follow a
standard schema even if that's at odds with submitters' preferences. Complicating this further, there is no current consensus on exactly how to annotate single
cell data, and any schema that Corpora adopts will evolve significantly over time.

Corpora will resolve this tension by storing and presenting two forms of submitted matrices. The first form, the "original" matrix, will follow the preferences
of submitters except where they conflict with a few fundamental requirements. The second form, the Corpora "remix", will follow a more extensive Corpora schema
that will aim to enable data integration. Having both matrices will allow Corpora to pursue both its publication and integration goals rather than having to
choose one over the other.

Maintaining multiple dataset forms does introduce some complexity into dataset validation and curation as well as the Corpora schema definition. Below, the
two-tiered approach to both is described.


## Validation and Curation

When a dataset is submitted to Corpora, Corpora's validation and curation processes should produce two forms of the dataset and verify that each meets
Corpora's requirements.


1. **Original** - This is the dataset with the metadata, gene identifiers, etc. exactly as it was received from the submitters. The original dataset is
most useful for publication, since it will match related journal articles or preprints as well as the submitters' opinions about how their results are
best presented.
2. **Remix** - This is the dataset intended for integration. It has been adapted to follow the schema for the Corpora remix. That schema will change
over time, but this adaptation involves changing metadata field names and values, inserting new metadata where necessary, mapping gene ids from one
annotation to another, introducing ontology values where appropriate, etc

Note that over time there may be more than one "remix" prepared for each dataset as integration use cases evolve. Also note that the "dataset form"
concept is not related to mere changes to the _format_ of the dataset, like converting from loom to an h5ad. 

The multiple forms of the dataset naturally mean that validation and curation split into multiple levels.

First, there is "basic" validation and curation that is applied to the original form of the dataset. This checks that Corpora's
[fundamental requirements](#Original-Dataset-Requirements) are met. These are things that are just absolutely needed for the system to work. For example,
Corpora will not host PII even if a submitter really wants it to. Second, there is the curation of the integration-focused remix dataset. This adapts
the dataset to all the requirements of the [Corpora remix schema](#Remix-Metadata).

Corpora's curation policies will be very strict regarding the fundamental requirements. Datasets that do not satisfy them will be rejected. But the
integration-focused curation will be more flexible. While curation is still driven by human interaction, we will strongly urge submitters to supply all
the required information. And we will also urge them to change their own datasets to match the Corpora schema so that the difference between the "original"
and the "remix" is very small. If submitters decline to do so, we'll talk with them about their reasons and record their feedback, but we will still add
the dataset to Corpora, and we will still create the remix variant of the submitted data.


### Submitter Incentives

We recognize that Corpora cannot _make_ submitters do anything; following our schema is voluntary. But, Corpora can create incentives. For example, if a
submitter doesn't want to include "tissue" then their dataset would not show up in tissue-related search results in the Portal. Or if a submitter uses a
gene annotation that is incompatible with the Corpora schema annotation, then gene card information will not be available for their dataset in cellxgene.

Over time we expect these incentives will lead to broader adoption of the Corpora schema and less friction around creating the remix datasets.


### Software vs Policy

Finally, we should draw a distinction between how we approach the schema and validation as a _policy_ versus how we implement it in _software_. The schema,
especially anything outside of fundamental requirements, will change over time. It might change frequently and significantly as we get feedback from submitters
and users. So while we might be strict about our validation policies at any given time, our software needs to be very flexible. We need to avoid a situation
where we can't evolve schema and curation policies because it requires expensive software updates. There will be situations where software features or
capabilities will be hard or impossible because we prioritize this flexibility.

## Schema

### Original Dataset Requirements


These are the fundamental requirements that all dataset must meet, the handful of things our software depends on and mandatory legal requirements. All
submitted data must comply with these requirements.


*   **Unique observation identifiers**. Each observation (usually a cell) must have an id that is unique within the dataset. This is usually already
    present in every submission, for example as a barcode.
*   **Unique feature identifiers**. Every feature (usually a gene or transcript) also needs a unique identifier. This is occasionally not present because
    of one-to-many mappings between gene symbols and other gene ids. In cases where there are duplicated feature identifiers, they will be summed during
    curation.
*   **Type annotations or inferrable types**. For both expression values and metadata, we don't want to have to guess at the type. Typing is already explicit
    in anndata, loom, Seurat and Bioconductor formats.
*   **No PII**. No metadata can be personally identifiable: no names, dates of birth, specific locations, etc. There's a
    [list](https://docs.google.com/document/d/1nlUuRiQ8Awo_QsiVFi8OWdhGONjeaXE32z5GblVCfcI/edit?usp=sharing).


#### Matrix Layers

The count matrix itself can exist in several forms. Generally, there is a "raw" count matrix without scaling, filtering, normalization, etc. And there is
a "final" matrix that is used in the publication with QC and various corrections. There can also be several intermediate matrices.

These different transformations of the count matrix are called "layers" in several file formats. Submissions to Corpora must contain at least one layer,
the "raw" layer. They may optionally also have a "final" layer used for presentation, and each layer must be identified.


### Remix Metadata

The remaining part of the schema defines what's needed for the Corpora remix. The words "must" and "requirement" are used throughout, but note the
discussion above regarding curation policies as well the [caution against](#Software-vs-policy) introducing any
software dependencies on these fields.


#### Integration Metadata

To support data integration, each cell must have the following metadata:

**Field name**|**Constraints**
:--|:--
tissue|string (see note below regarding cell culture and organoids)
assay|string
disease|string
cell\_type|string
sex|"male", "female", "mixed", "unknown"
ethnicity|string, "na" if non-human, "unknown" if not available
development\_stage|string, "unknown" if not available

The `tissue` field must be appended with "(cell culture)" or "(organoid)" if appropriate.
In addition to these free text fields (except sex), each cell must also have ontology values:

**Field name**|**Constraints**
:--|:--
tissue\_ontology\_term\_id|UBERON term
assay\_ontology\_term\_id|EFO term
disease\_ontology\_term\_id|MONDO term or [PATO:0000461](http://bioportal.bioontology.org/ontologies/PATO?p=classes&conceptid=PATO%3A0000461)
cell\_type\_ontology\_term\_id|CL term
ethnicity\_ontology\_term\_id|HANCESTRO term, "na" if non-human
development\_stage\_ontology\_term\_id|HsapDv term if human, child of EFO:0000399 otherwise

The `tissue_ontology_term_id` field must be appended with "(cell culture)" or "(organoid)" if appropriate.

The `ontology_term_id` fields may be empty strings when no appropriate ontology value is available, for example if the
`cell_type` field describes a cell type for which not ontology term exists or if the `ethnicity` is unknown. If the
field is not empty, then it must be an [OBO-format ID](http://www.obofoundry.org/id-policy.html), meaning it is a CURIE
where the prefix identifies the ontology. For example `EFO:0000001` is a term in the `EFO` ontology.

If ontology terms are missing in the submitted dataset, then as part of preparing the remix dataset, Corpora will insert
appropriate ontology terms based on text from the submitters.

If the features of the dataset are human genes, then the feature ids must be [HGNC](https://www.genenames.org/about/guidelines/#!/#tocAnchor-1-7) approved
gene symbols. If the original dataset does not use HGNC symbols, Corpora will perform a conversion from the submitted features to HGNC symbols. In
cases where multiple gene symbols are merged into a single symbol, Corpora will sum the counts in the raw data layer.

Similarly if the features are mouse genes, then the feature ids must be [MGI](http://www.informatics.jax.org/mgihome/nomen/gene.shtml) gene symbols.

Finally, the whole dataset must be annotated with fields that indicate the organism and describe the meanings of the [matrix layers](#Matrix-Layers):

**Field name**|**Constraints**
:--|:--
organism|String
organism\_ontology\_term\_id|NCBITaxon term
layer\_descriptions|Mapping from {layer\_name: layer\_description, ...} Each description is free text, though one layer must be described as "raw".

As with the ontology fields above, the ontology is added as part of Corpora curation.


#### Presentation Metadata

There are also fields that are required so that Corpora (both Portal and cellxgene) can present datasets appropriately.

Each dataset must have at least one **embedding**, a mapping from each cell to a tuple of floats of length at least 2. These are usually generated by algorithms
like umap or tsne and are used to display the dataset in cellxgene. If a submitted dataset does not have an embedding, the umap embedding will be calculated
and inserted during curation.

Datasets also must have a few other metadata fields for presentation:

**Field name**|**Description**
:--|:--
title|String that identifies the dataset
contributors|A list of {"name": ..., "institution": ..., "email": ...}, where each value is a string. Only "name" is required.


If the dataset has an associated preprint or journal publication, these fields are required:

**Field name**|**Description**
:--|:--
preprint\_doi|DOI of the associated preprint.
publication\_doi|DOI of the associated journal publication.


#### Presentation Hints

The metadata fields below are fully optional. They aren't needed for integration, and Corpora can present the data fine without them, but if they are
present Corpora will do something with them. This allows submitters to fine-tune how their datasets are presented, which is a common request.


<table>
  <tr>
   <td><strong>Field name</strong>
   </td>
   <td><strong>Description</strong>
   </td>
  </tr>
  <tr>
   <td>color_map
   </td>
   <td>Submitters can include a field called "{field}_colors" for any other categorical integer metadata field. The value must be an array of one of
       fourcolor specifications:
<ul>

<li>CSS4 color name, as supported by matplotlib
    <a href="https://matplotlib.org/3.1.0/gallery/color/named_colors.html">https://matplotlib.org/3.1.0/gallery/color/named_colors.html</a>

<li>RGB tuple/list with values ranging from 0.0 to 1.0, as in [0.5, 0.75, 1.0]

<li>RFB tuple/list with values ranging from 0 to 255, as in [128, 192, 255]

<li>Hex triplet string, as in "#08c0ff"and each string must be a hex color code.

</li>
</ul>
The color code at the nth position in the array corresponds to category n in the metadata field.

   </td>
  </tr>
  <tr>
   <td>tags
   </td>
   <td>A list of strings that the Portal or other tools may use for searching, sorting, or filtering.
   </td>
  </tr>
  <tr>
   <td>default_field
   </td>
   <td>Name of another field that should be the "default" for display. So for example, when cellxgene starts, this is the field that the cells should be
       colored by.
   </td>
  </tr>
  <tr>
   <td>default_embedding
   </td>
   <td>Name of the embedding that should be the default for display.
   </td>
  </tr>
</table>


Additionally, submitters may include project-level metadata inside their dataset. This information is used to populate project metadata in Corpora when
the submitted dataset is used to create a new project. 


**Field name**|**Description**
:--|:--
project\_name|String
project\_description|Longer description of the project.
project\_links|List of links associaed with the project.

Much of the project metadata is links. Each link is a dictionary with keys "link_url", "link_name", and "link_type".
"link_name" gives the text that should be displayed for the link. "link_type" is one of `PROTOCOL`, `RAW_DATA`,
`SUMMARY`, or `OTHER`. Exactly one link should be of type `SUMMARY`.

#### System Fields

Finally, Corpora inserts a **corpora_dataset_revision** field for internal bookkeeping. This must not be present in the original submitted matrix, and Corpora
will insert it during curation. Corpora also reserves any metadata fields that start with "corpora_" for future system use, so submitted datasets must not
have any such fields.


## **Implementations**

Corpora requires submitted count matrices and associated metadata to be in one of three formats: AnnData, Loom, or a Seurat v3 RDS. Other formats are rejected.
Each of these formats has a way to include metadata along with the count data, so a submission can be entirely contained within a single file. 


### **Scanpy/AnnData**

[Anndata](https://anndata.readthedocs.io/en/stable/) is the hdf5-based file format used by scanpy. There is a python library for interacting with it. The
count data is stored in an attribute `X` of shape (# of cells, # of genes). `X` can either be a numpy.ndarray or a scipy.sparse.spmatrix.

Information about cells and genes are stored in `obs` and `var` dataframes, respectively. Each of those has an index which can serve as cell and gene ids
if all the values are unique.

The `obs` dataframe must store all the cell-level metadata except for embeddings. So, there must be columns named "tissue", "assay", "disease", and
"cell_type".

Embeddings are stored in `obsm` with a key name of "X_{description}" where description  provides some information about how the embedding was generated,
X_umap for example.

Finally, the dataset-level metadata is stored in `uns`, which is just key-value pairs.

Note that anndata supports "layers" and "raw" values for counts. Those are permitted, but Corpora will treat `X` as the "final" matrix for further analysis
and visualization. Once anndata [unifies its treatment of layers](https://github.com/theislab/anndata/issues/244), Corpora will use the "default" as the
final matrix, however that ends up being specified. In anndata files, the layer_descriptions dictionary should have a key "X" and optionally "raw.X" to
describe those layers.


### **Loom**

[Loom](http://linnarssonlab.org/loompy/format/index.html) is also hdf5-based. The count data is stored in the "main matrix" of shape
(# of genes, # of cells), which is the `/matrix` hdf5 dataset in the file.

Information about cells and genes are stored in the hdf5 `/col_attrs` and `/row_attrs` groups, respectively, and are available via the `ca` and `ra`
python attributes. The unique id for cells must be stored in the CellID column attribute, while the unique id for genes must be stored in the "Accession"
row attribute.

The column attributes must also contain values for "tissue", "assay", "disease", and "cell_type".

Embeddings are not handled in a special way, they are just more column attributes. So, embeddings should be stored as an ndarray under a key in column
attributes that ends with "_embeddings".

actuallyDataset-level metadata is stored in the loom file's `/attrs` hdf5 group, and are accessible as a mapping via the `attrs` attributie in the python
interface.

Like anndata, loom supports layers, but Corpora will focus on the "main matrix" as the "final" matrix for visualization.


### **Seurat**

[Seurat](https://github.com/satijalab/seurat/wiki) is an R library with a bunch of methods for single-cell data analysis. It defines a "Seurat object"
which contains all the data and metadata for a dataset. This object is commonly serialized using R's saveRDS and readRDS functions. There are two major
versions of Seurat objects: 2.0 and 3.0. Corpora requires 3.0 objects.

Seurat objects can store multiple "assays" with the default identified as the "`active.assay`". Corpora will only use the data from the active assay.
Expression data is stored in the data slot of an assay, and is shaped (# of genes, # of cells). Cell and gene ids are available as `colnames` and `rownames`
of the object, respectively.

Cell-level metadata is stored in the `meta.data` slot of the Seurat object, and that should have columns named "tissue", "assay", "disease", and "cell_type".

Embeddings are stored in the `reductions` slot of the object, and each embedding is a `DimReduc` object. The embedding values themselves are stored in
`cell.embeddings` slot.

Finally, Seurat objects have a misc slot where the dataset-level metadata must be stored. The misc slot can store values in its own slots and functions
as a key-value store.
