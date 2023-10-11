import yaml

default_config = """
server:
  app:
    verbose: false
    debug: false
    host: localhost
    port : null
    open_browser: false
    force_https: false
    flask_secret_key: null
    generate_cache_control_headers: false
    server_timing_headers: false
    csp_directives: null

    # By default, cellxgene will serve api requests from the same base url as the webpage.
    # In general api_base_url and web_base_url will not need to be set.
    # Reasons to set these parameters:
    #  1. The cellxgene deploymnent is in an environment where the webpage and api have
    #     different base urls.  In this case both api_base_url and web_base_url must be set.
    #     It is up to the server admin to ensure that the networking is setup correctly for this environment.
    api_base_url: null
    web_base_url: null

  multi_dataset:
    # dataroot or dataroots must be set, then cellxgene may serve multiple datasets.
    # dataroot must be a string, representing the path to a directory or S3 prefix.  In this
    # case the datasets in that location are accessed from <server>/d/<datasetname>.
    # example:
    #     dataroot:  /path/to/datasets/
    # or
    #     dataroot:  s3://bucket/prefix/
    #
    # dataroots can be a dictionary, where a dataset key is associated with a base_url
    # and a dataroot.
    # example:
    #     dataroots:
    #         d1:
    #            base_url: set1
    #            dataroot: /path/to/set1_datasets/
    #         d2:
    #            base_url: set2/subdir
    #            dataroot: /path/to/set2_datasets/
    #
    # In this case, datasets can be accessed from <server>/set1/<datasetname> or
    # <server>/set2/subdir/<datasetname>.

    dataroot: null

    # The index page when in multi-dataset mode:
    #   false or null:  this returns a 404 code
    #   true:  loads a test index page, which links to the datasets that are available in the dataroot
    #   string/URL:  redirect to this URL:  flask.redirect(config.multi_dataset__index)
    index: false


  data_locator:
    api_base: null
    # s3 region name.
    #   if true, then the s3 location is automatically determined from the dataroot.
    #   if false/null, then do not set.
    #   if a string, then use that value (e.g. us-east-1).
    s3_region_name: true

  gene_info:
    api_base: null

  adaptor:
    cxg_adaptor:
      # The key/values under tiledb_ctx will be used to initialize the tiledb Context.
      # If 'vfs.s3.region' is not set, then it will automatically use the setting from
      # data_locator / s3 / region_name.
      tiledb_ctx:
        sm.tile_cache_size:  8589934592  # 8GiB
        py.init_buffer_bytes: 536870912  # 512MiB

  limits:
    column_request_max: 32
    diffexp_cellcount_max: null


default_dataset:
  app:
    # Scripts can be a list of either file names (string) or dicts containing keys src, integrity and crossorigin.
    # these will be injected into the index template as script tags with these attributes set.
    scripts: []
    # Inline scripts are a list of file names, where the contents of the file will be injected into the index.
    inline_scripts: []

    about_legal_tos: null
    about_legal_privacy: null

  presentation:
    max_categories: 1000
    custom_colors: true

  embeddings:
    names : []

  diffexp:
    enable: true
    lfc_cutoff: 0.01
    top_n: 10
    count: 15

  X_approximate_distribution: normal # currently fixed config
"""


def get_default_config():
    return yaml.load(default_config, Loader=yaml.Loader)
