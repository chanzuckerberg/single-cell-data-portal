# Corpora Schema AnnData Implementation

Corpora datasets are stored in the HDF5-backed AnnData format as written by version 0.7 of the anndata library. One
goal of using this format is to have all metadata stored in the same file as the data, so cellxgene will have access to
all Corpora metadata within a single file.

The Corpora schema requirements and definitions for the AnnData `X`, `uns`, `obs`, and `obsm` attributes are below.


### Note on types
The types below are python3 types. Note that a python3 `str` is a sequence of Unicode code points, which is stored
null-terminated and UTF-8-encoded by anndata.

For succinctness, there is a `link` type used below that is a `dict` {"link_name": `str`, "link_url": `str`}.

## [`X`](#X)

Corpora does not impose any additional constraints on the `X` data matrix. So it may be sparse or dense and any
numeric `numpy.dtype`.

## [`uns`](#uns)

Recall that `uns` is a mapping with `str`s as keys. Corpora requires the following keys and values in `uns`:

**Key**|**Value Type**|**Notes**
-----|-----|-----
title|`str`|Title of the dataset
contributors|`list` of `dict`s. each `dict` must a "name": `str` item and may optionally have "email": `str` and "institution": `str`|This is the person who gave us the dataset, not necessarily the author or actual data generator.
layer\_descriptions|dict with string keys and values|One key must be "X" which describes the transformations (if any) performed to produce the X matrix cellxgene displays.
organism|`str`| 
organism\_ontology|`str`| 
project\_name|`str`|Name of the project that contains this dataset
project\_description|`str`| 
project\_protocol\_links|`list` of `link`s|May be empty list
project\_raw\_data\_links|`list` of `link`s|May be empty list
project\_other\_links|`list` of `link`s|May be empty list

The `project\_protocol\_links`, `project\_raw\_data\_links`, and `project\_other\_links` keys are always present in `uns` as
they are used in the creation of "projects during Corpora ingestion. But if the values for these keys are empty lists, then
wrt creating a user interface, that is equivalent to the fields being missing.

Corpora defines the following optional keys and values in `uns`. If they key is present, then the value must not be empty.

**Key**|**Value Type**|**Notes**
-----|-----|-----
preprint\_doi|`str`|Suitable for submission to the CrossRef API
publication\_doi|`str`|Suitable for submission to the CrossRef API
default\_embedding|`str`|Must match a key to an embedding in `obsm`.
default\_field|`str`|Must match a column name in `obs`
tags|`list` of `str`s| 
<obs\_column>\_colors|`list` of color specifications (see anndata/cellxgene documentation)|<obs\_column> must be a column name in `obs`. There may be multiple keys like this.


## [`obsm`](#obsm)

`obsm` is a mapping from `str`s to matrices of shape (n, m) where `n` is the number of rows in `X` and `m >= 1`.
Corpora requires one value in `obms` to be a two-dimensional embedding, meaning it is of shape (n, 2) and is a
numeric type.

## [`obs`](#obs)

Corpora requires a number of columns be present in the `obs` dataframe:

**Key**|**Value Type**|**Notes**
-----|-----|-----
tissue|`str`| 
assay|`str`| 
disease|`str`| 
cell\_type|`str`| 
sex|`str`| 
ethnicity|`str`| 
development\_stage|`str`| 
tissue\_ontology|`str`|Not necessarily human-readable
assay\_ontology|`str`|Not necessarily human-readable
disease\_ontology|`str`|Not necessarily human-readable
cell\_type\_ontology|`str`|Not necessarily human-readable
ethnicity\_ontology|`str`|Not necessarily human-readable
development\_stage\_ontology|`str`|Not necessarily human-readable
