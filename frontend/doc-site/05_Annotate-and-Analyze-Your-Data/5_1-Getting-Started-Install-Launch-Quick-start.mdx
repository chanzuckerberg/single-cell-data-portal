# Getting Started: Install, Launch, Quick Start


## Install

Cellxgene desktop has two parts:

``cellxgene`` is the main explorer application, which takes an already-processed h5ad file as input. This is installed by default.

``cellxgene prepare`` provides auxiliary functionality for preparing your dataset. This is not installed by default.


### Requirements

You'll need python 3.6+ and an up-to-date version of Google Chrome. The web UI is tested on OSX and Windows using Chrome, and the python CLI is tested on OSX and Ubuntu (via WSL/Windows). It should work on other platforms, but if you run into trouble let us know [INSERT LINK TO OUR ONLINE COMMUNITY HERE]

Python.org has help on installing a recent version of Python, including the pip package manager. Chrome is available at [https://www.google.com/chrome/downloads/](https://www.google.com/chrome/downloads/).


### Basic Install Using Pip ▼

To install the `cellxgene` explorer alone, run:


```
pip install cellxgene
```


To install `cellxgene` and the optional `cellxgene prepare`, run:


```
pip install cellxgene[prepare]
```


Note: if the aforementioned optional prepare package installation fails, you can also install these packages directly:


```
pip install scanpy>=1.3.7 python-igraph louvain>=0.6
```


On Linux platforms, you may also need to install build dependencies first:

```
sudo apt-get install build-essential python-dev
pip install scanpy>=1.3.7 python-igraph louvain>=0.6
```

If you already have `cellxgene` installed, you can update to the most recent version by running:


```
pip install cellxgene --upgrade
```



### Using a conda environment ▼

To install `cellxgene` alone, run:

```
conda create --yes -n cellxgene python=3.7
conda activate cellxgene
pip install cellxgene
```

To install` cellxgene` and the optional` cellxgene prepare`, run:

```
conda create --yes -n cellxgene python=3.7
conda activate cellxgene
pip install cellxgene[prepare]
```

### Using a virtual environment  ▼

To install cellxgene alone, run:

```
ENV_NAME=cellxgene
python3.7 -m venv ${ENV_NAME}
source ${ENV_NAME}/bin/activate
pip install cellxgene
```

To install cellxgene and cellxgene prepare, run:


```
ENV_NAME=cellxgene
python3.7 -m venv ${ENV_NAME}
source ${ENV_NAME}/bin/activate
pip install cellxgene[prepare]
```

### Using docker ▼

Build the image

```
docker build . -t cellxgene
```


Run the container and mount data (change data location, `--port` and `--host` parameters as needed). You will need to use `--host 0.0.0.0` to have the container listen to incoming requests from the browser.

```
docker run -v "$PWD/example-dataset/:/data/" -p 5005:5005 cellxgene launch --host 0.0.0.0 data/pbmc3k.h5ad
```

## Launch


### Launching cellxgene with your dataset ▼

Once you've prepared [link to prepare your data page] your data for cellxgene, you can launch the app using:


```
cellxgene launch mydataset.h5ad --open
```



### Launching from a URL▼

You can also launch from a URL directly like this:


```
cellxgene launch https://github.com/chanzuckerberg/cellxgene/raw/main/example-dataset/pbmc3k.h5ad
```


Support for S3 and GCS is not enabled by default. If you wish to directly access S3 or GFS, install one or both of the following packages:

-  [s3fs ](https://s3fs.readthedocs.io/en/latest/)for S3 support

- [gcsfs ](https://gcsfs.readthedocs.io/en/latest/)for GCS support

For example:

```
pip install s3fs
cellxgene launch s3://mybucket.s3-us-west-2.amazonaws.com/mydata.h5ad
```


### Options for cellxgene launch ▼

For the most up-to-date and comprehensive list of options, run cellxgene launch --help

`--open` automatically opens the web browser after launching (caveat: only works on some operating systems).

`--disable-annotations`, `--annotations-file` & `--annotations-dir` all have to do with creating new categorical annotations in the application. We have a whole separate page about their usage! :) [LINK TO ANNOTATE DATA PAGE]

`--diffexp-lfc-cutoff` as explained in the methods[LINK TO ALGORITHMS USED BY CXG PAGE] , genes are only returned in differential expression if the effect size is above the specified threshold for log fold change. Defaults to 0.01.

`--disable-diffexp` will disable and hide the Compute Differential Expression feature. For large datasets, or datasets loaded with the `--backed` option, computing differential expression may be extremely slow or use excessive resources on the host computer (e.g., memory thrashing). Disabling the feature will ensure that this computation is not initiated accidentally.

`--backed` option instructs cellxgene launch to read the H5AD file in "backed" mode

By default, cellxgene will read the entire H5AD will be into memory at startup, improving application speed and performance. Very large datasets may not fit in memory. The "--backed" mode will read the file incrementally, reducing memory use, and for large files, improving startup speed. However, this option will also significantly slow down access to gene expression histograms, and may render differential expression calculations too slow to use (see `--disable-diffexp` for an option to disable this feature).

`--embedding` restricts which embeddings will be available in the viewer. By default, all embeddings specified in `anndata.obsm['X_name']` will be loaded; if you have many embeddings, you may wish to restrict this list for a speedier launch.

`--title` adds a title to the viewer. Defaults to file name.

`--about` adds a link where users can go to find more infomation about the dataset. Requires `https` .

`--obs-names` allows you to specify which column in anndata.obs to use as anndata.obs.index .

`--var-names` allows you to specify which column in `anndata.var` to use as anndata.var.index .

`--max-category-items` omits categorical metadata fields that contain more than N _distinct_ values. Defaults to 1000.


## Quick Start


### Overview

Cellxgene Desktop Explorer is a tool in the cellxgene suite that enables collaborative exploratory analysis of single cell data and private data sharing.

This page describes how to quickly get started with exploratory analysis on your local computer. Those interested in hosting cellxgene to enable your collaborators to explore data from their browser without needing to download the data file or install cellxgene should see Self-Hosting CellxGene[LINK TO SELF-HOSTING CXG PAGE ON DOCS SITE ].


### Quick Start

To install CellxGene[LINK TO INSTALL SECTION ABOVE], you need Python 3.6+. We recommend installing cellxgene into a conda[LINK TO CONDA INSTRUCTIONS ABOVE]  or virtual environment[LINK TO VIRTUAL ENVIRONMENT INSTRUCTIONS ABOVE].

Install the package, then launch cellxgene with an example [AnnData](https://anndata.readthedocs.io/en/latest/) file.


### Example datasets

The following datasets are available on the web and can be opened directly using cellxgene launch:


#### Peripheral blood mononuclear cells

[embed table]


```
cellxgene launch https://cellxgene-example-data.czi.technology/pbmc3k.h5ad
```



#### Tabula muris

[embed table]


```
cellxgene launch https://cellxgene-example-data.czi.technology/tabula-muris.h5ad
```



#### Tabula muris senis

[embed table]


```
cellxgene launch https://cellxgene-example-data.czi.technology/tabula-muris-senis.h5ad
```
