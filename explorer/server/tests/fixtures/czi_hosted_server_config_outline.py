f"""server:
  app:
    verbose: {verbose}
    debug: {debug}
    host: {host}
    port: {port}
    open_browser: {open_browser}
    force_https: {force_https}
    flask_secret_key: {flask_secret_key}
    generate_cache_control_headers: {generate_cache_control_headers}
    server_timing_headers: {server_timing_headers}
    csp_directives: {csp_directives}
    api_base_url: {api_base_url}
    web_base_url: {web_base_url}

  multi_dataset:
    dataroot: {dataroot}
    dataroots: {dataroots}
    index: {index}

  data_locator:
    api_base: {data_locator_api_base}
    s3_region_name: {data_locator_region_name}

  gene_info:
    api_base: {gene_info_api_base}

  adaptor:
    cxg_adaptor:
      tiledb_ctx:
        sm.tile_cache_size:  {cxg_tile_cache_size}
        py.init_buffer_bytes: {cxg_tiledb_py_init_buffer_size}

  limits:
    column_request_max: {column_request_max}
    diffexp_cellcount_max: {diffexp_cellcount_max}
"""
