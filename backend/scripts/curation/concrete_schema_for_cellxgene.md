# Corpora Schema AnnData Implementation

Corpora datasets are stored in the HDF5-backed AnnData format as written by version 0.7 of the anndata library. One
goal of using this format is to have all metadata stored in the same file as the data, so cellxgene will have access to
all Corpora metadata within a single file.

The Corpora schema requirements and definitions for the AnnData `X`, `uns`, `obs`, and `obsm` attributes are below.


### Note on types
The types below are python3 types. Note that a python3 `str` is a sequence of Unicode code points, which is stored
null-terminated and UTF-8-encoded by anndata.

For succinctness, there is a `link` type used below that is a `dict` {"link_name": `str`, "link_url": `str`}.

## `X`

Corpora does not impose any additional constraints on the `X` data matrix. So it may be sparse or dense and any
numeric `numpy.dtype`.

## `uns`

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
project\_protocol\_links|`list` of `link`s|May be empty
project\_raw\_data\_links|`list` of `link`s|May be empty
project\_other\_links|`list` of `link`s|May be empty

Corpora defines the following optional keys and values in `uns`:

**Key**|**Value Type**|**Notes**
-----|-----|-----
preprint\_doi|`str`|Suitable for submission to the CrossRef API
publication\_doi|`str`|Suitable for submission to the CrossRef API
default\_embedding|`str`|Must match a key to an embedding in `obsm`.
default\_field|`str`|Must match a column name in `obs`
tags|`list` of `str`s| 
<obs\_column>\_colors|`list` of color specifications (see anndata/cellxgene documentation)|<obs\_column> must be a column name in `obs`. There may be multiple keys like this.
