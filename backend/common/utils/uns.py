def filter_spatial_data(content, library_id):
    """
    This filters data associated with the "spatial" key in a dictionary, specifically
    retaining certain sub-items from "images" and "scalefactors" sub-dictionaries.
    https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md#uns-dataset-metadata
    """
    spatial_filtered = {}
    spatial_filtered[library_id] = {
        "images": {
            "hires": content["images"]["hires"],  # Omit hires data once deep zooming feature is implemented
            "fullres": [],  # Currently not including fullsres data, due to deep zooming feature coming soon
        },
        "scalefactors": {
            "spot_diameter_fullres": content["scalefactors"]["spot_diameter_fullres"],
            "tissue_hires_scalef": content["scalefactors"]["tissue_hires_scalef"],
        },
    }
    return spatial_filtered
