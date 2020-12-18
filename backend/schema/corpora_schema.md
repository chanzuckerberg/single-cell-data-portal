# cellxgene Data Integration Schema

Authors: mkinsella@chanzuckerberg.com

Document Status: _Approved_

Version: 1.1.0

Date Last Modified: 2020-12-11


## Background

cellxgene aims to support the publication and sharing of single-cell datasets as well as the construction of a data corpus that facilitates data integration
across multiple tissues and experiments. Achieving the latter goal requires some harmonization of metadata and features in the cellxgene Data Portal. But if
that harmonization is too onerous, it will burden the goal of rapid data sharing.

We balance these two goals by requiring datasets hosted in the cellxgene Data Portal to follow a small schema with only a few required fields. These fields
are expected to be very useful for data integration and also simply and readily available from data submitters.

Note that the requirements in the schema are just the minimum required information. We expect that datasets will have additional metadata, and this
will be preserved in datasets submitted to the Data Portal.

## Curation and Validation

When a submitter is preparing a dataset for the cellxgene Data Portal, they will be able to use a command line tool to apply changes to their dataset so that it
follows the schema. When that tool successfully completes, it will give the submitter a message that the dataset is ready to be submitted to the Data Portal,
and it will write the schema version number into the dataset's metadata.

Then, when the submitter uploads the dataset to the Data Portal, the Data Portal will verify that the dataset does indeed follow an accepted version of the
schema. If it does not, it will reject the dataset with an appropriate error message.

## Schema

### Basic Requirements

There are a few requirements that are needed for the cellxgene Explorer to work correctly or are mandatory legal requirements:

*   **Unique observation identifiers**. Each observation (usually a cell) must have an id that is unique within the dataset. This is usually already
    present in every submission, for example as a barcode.
*   **Unique feature identifiers**. Every feature (usually a gene or transcript) also needs a unique identifier. This is occasionally not present because
    of one-to-many mappings between gene symbols and other gene ids. In cases where there are duplicated feature identifiers, they will need to be appropriately
    combined before submission. For example, raw counts will be summed and logged counts will be exponentiated, summed, and logged.
*   **No PII**. No metadata can be personally identifiable: no names, dates of birth, specific locations, etc. There's a
    [list](https://docs.google.com/document/d/1nlUuRiQ8Awo_QsiVFi8OWdhGONjeaXE32z5GblVCfcI/edit?usp=sharing).

### Matrix Layers

The count matrix itself can exist in several forms. Generally, there is a "raw" count matrix without scaling, filtering, normalization, etc. And there is
a "final" matrix that is used in the publication with QC and various corrections. There can also be several intermediate matrices. We require that these matrices
have the same dimensions, so the raw count matrix should include the same cells and genes as the final.

These different transformations of the count matrix are called "layers" in several file formats. Submissions to the Data Portal must contain the "raw"
layer. Note that AnnData objects have an attribute called [raw](https://anndata.readthedocs.io/en/latest/anndata.AnnData.raw.html#anndata.AnnData.raw). This
could be a natural place to store the raw count matrix, but as long as the locations of the layers are specified, it can go anywhere.


### Schema Version

Datasets in the Data Portal must store the version of the schema they follow (that is, the version of this document) as
well as the version of the particular encoding used. The encoding is documented
[elsewhere](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/schema/corpora_schema_h5ad_implementation.md) and describes techincal details
of how the schema should be serialized in a particular file format.

**Field name**|**Constraints**
:--|:--
corpora_schema_version|Follows [SemVer](https://semver.org/) conventions
corpora_encoding_version|Follows [SemVer](https://semver.org/) conventions

### Integration Metadata

To support data integration, each cell must have the following metadata:

**Field name**|**Constraints**
:--|:--
tissue|string (see note below regarding cell culture and organoids)
assay|string
disease|string
cell\_type|string
sex|"male", "female", "mixed", "unknown", or "other"
ethnicity|string, "na" if non-human, "unknown" if not available
development\_stage|string, "unknown" if not available

The `tissue` field must be appended with " (cell culture)" or " (organoid)" if appropriate. Also, if the source of cells is cell culture or organiod, the
`cell_type` field can be left empty.
In addition to these free text fields (except sex), each cell must also have ontology values:

**Field name**|**Constraints**
:--|:--
tissue\_ontology\_term\_id|UBERON term
assay\_ontology\_term\_id|EFO term
disease\_ontology\_term\_id|MONDO term or [PATO:0000461](http://bioportal.bioontology.org/ontologies/PATO?p=classes&conceptid=PATO%3A0000461)
cell\_type\_ontology\_term\_id|CL term
ethnicity\_ontology\_term\_id|HANCESTRO term, "na" if non-human
development\_stage\_ontology\_term\_id|HsapDv term if human, child of EFO:0000399 otherwise

The `tissue_ontology_term_id` field must be appended with " (cell culture)" or " (organoid)" if appropriate.

The `ontology_term_id` fields may be empty strings when no appropriate ontology value is available, for example if the
`cell_type` field describes a cell type for which no ontology term exists or if the `ethnicity` is unknown. If the
field is not empty, then it must be an [OBO-format ID](http://www.obofoundry.org/id-policy.html), meaning it is a CURIE
where the prefix identifies the ontology. For example `EFO:0000001` is a term in the `EFO` ontology.

If the features of the dataset are human genes, then the feature ids must be [HGNC](https://www.genenames.org/about/guidelines/#!/#tocAnchor-1-7) approved
gene symbols.

Similarly if the features are mouse genes, then the feature ids must be [MGI](http://www.informatics.jax.org/mgihome/nomen/gene.shtml) gene symbols.

Finally, the whole dataset must be annotated with fields that indicate the organism and describe the meanings of the [matrix layers](#Matrix-Layers):

**Field name**|**Constraints**
:--|:--
organism|String
organism\_ontology\_term\_id|NCBITaxon term
layer\_descriptions|Mapping from {layer\_name: layer\_description, ...} Each description is free text, though one layer must be described as "raw".


### Presentation Metadata

There are also fields that are required so that the cellxgene Data Portal and Explorer can present datasets appropriately.

Each dataset must have at least one **embedding**, a mapping from each cell to a tuple of floats of length at least 2. These are usually generated by algorithms
like umap or tsne and are used to display the dataset in the Explorer. cellxgene provides a command line tool to add embeddings if they are missing.

Datasets also must have a few other metadata fields for presentation:

**Field name**|**Description**
:--|:--
title|String that identifies the dataset


#### Presentation Hints

The metadata fields below are optional. They aren't needed for integration, and cellxgene can display the data fine without them, but if they are
included cellxgene will do something with them. This allows submitters to fine-tune how their datasets are presented, which is a common request.


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


## **Implementations**

The Data Portal requires submitted count matrices and associated metadata to be in one of three formats: AnnData, Loom, or a Seurat v3 RDS. Other formats are rejected.
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

Note that anndata supports "layers" and "raw" values for counts. Those are permitted, but cellxgene will treat `X` as the "final" matrix for further
visualization. Once anndata [unifies its treatment of layers](https://github.com/theislab/anndata/issues/244), cellxgene will use the "default" as the
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

Dataset-level metadata is stored in the loom file's `/attrs` hdf5 group, and are accessible as a mapping via the `attrs` attributie in the python
interface.

Like anndata, loom supports layers, but cellxgene will focus on the "main matrix" as the "final" matrix for visualization.


### **Seurat**

[Seurat](https://github.com/satijalab/seurat/wiki) is an R library with a bunch of methods for single-cell data analysis. It defines a "Seurat object"
which contains all the data and metadata for a dataset. This object is commonly serialized using R's saveRDS and readRDS functions. There are two major
versions of Seurat objects: 2.0 and 3.0. cellxgene requires 3.0 objects.

Seurat objects can store multiple "assays" with the default identified as the "`active.assay`". cellxgene will only use the data from the active assay.
Expression data is stored in the data slot of an assay, and is shaped (# of genes, # of cells). Cell and gene ids are available as `colnames` and `rownames`
of the object, respectively.

Cell-level metadata is stored in the `meta.data` slot of the Seurat object, and that should have columns named "tissue", "assay", "disease", and "cell_type".

Embeddings are stored in the `reductions` slot of the object, and each embedding is a `DimReduc` object. The embedding values themselves are stored in
`cell.embeddings` slot.

Finally, Seurat objects have a misc slot where the dataset-level metadata must be stored. The misc slot can store values in its own slots and functions
as a key-value store.
