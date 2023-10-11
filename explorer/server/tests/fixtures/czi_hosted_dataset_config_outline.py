f"""
default_dataset:
  app:
    scripts: {scripts} #list of strs (filenames) or dicts containing keys
    inline_scripts: {inline_scripts} #list of strs (filenames)

    about_legal_tos: {about_legal_tos}
    about_legal_privacy: {about_legal_privacy}

  presentation:
    max_categories: {max_categories}
    custom_colors: {custom_colors}

  embeddings:
    names: {embedding_names}

  diffexp:
    enable: {enable_difexp}
    lfc_cutoff: {lfc_cutoff}
    top_n: {top_n}

  X_approximate_distribution: {X_approximate_distribution}
"""
