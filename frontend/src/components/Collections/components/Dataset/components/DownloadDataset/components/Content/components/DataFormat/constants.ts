/**
 * @deprecated by TOOLTIP_SLOT_PROPS constant once "DOWNLOAD_UX" feature flag is removed (#5566).
 */
export const DOWNLOAD_TOOLTIP_SLOT_PROPS = {
  arrow: {
    style: {
      color: "#cd4246",
    },
  },
  tooltip: {
    style: {
      backgroundColor: "#cd4246",
      borderRadius: 2,
      boxShadow:
        "0 0 0 1px rgba(17,20,24,.1), 0 2px 4px rgba(17,20,24,.2), 0 8px 24px rgba(17,20,24,.2)",
      fontSize: 14,
      fontWeight: 400,
      lineHeight: "inherit",
      maxWidth: "unset", // Override the max-width specification for dark sdsStyle.
      padding: "8px 9.6px",
    },
  },
};

/**
 * @deprecated by TOOLTIP_TITLE constant once "DOWNLOAD_UX" feature flag is removed (#5566).
 */
export const DOWNLOAD_TOOLTIP_TITLE =
  "A .rds (Seurat v4) download is unavailable due to limitations in the R dgCMatrix sparse matrix class.";

export const TOOLTIP_SLOT_PROPS = {
  tooltip: {
    style: {
      maxWidth: 332, // Override the max-width specification for dark sdsStyle.
    },
  },
};

export const TOOLTIP_TITLE =
  ".rds (Seurat v4) is unavailable due to limitations in the R dgCMatrix sparse matrix class.";
