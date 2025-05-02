import React, { FC, SyntheticEvent, useEffect, useMemo, useState } from "react";
import { Dataset, DATASET_ASSET_FORMAT } from "src/common/entities";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  DialogActions,
  DialogContent,
  DialogTitle,
  Tab,
  Tooltip,
} from "@czi-sds/components";
import { Button } from "src/components/common/Button";
import {
  POSSIBLE_DOWNLOAD_FORMATS,
  getNotAvailableText,
} from "./components/BrowserTab/components/DataFormat/constants";
import { downloadMultipleFiles, getDownloadLink } from "./utils";
import { BrowserTab } from "./components/BrowserTab/BrowserTab";
import { ActiveTab } from "./utils";
import { ApiTab } from "./components/ApiTab/ApiTab";
import { StyledTabs } from "./style";
interface Props {
  isError?: boolean;
  isLoading?: boolean;
  onClose: () => void;
  name: string;
  dataAssets: Dataset["dataset_assets"];
  collectionId?: string;
}

export interface DownloadLinkType {
  filetype: DATASET_ASSET_FORMAT;
  fileSize: number | undefined;
  downloadURL: string | undefined;
}

const Content: FC<Props> = ({
  isError = false,
  isLoading = false,
  onClose,
  name,
  dataAssets,
  collectionId,
}) => {
  const [activeTab, setActiveTab] = useState<ActiveTab>(ActiveTab.Browser);
  const [initialFetchDone, setInitialFetchDone] = useState(false);

  const [selectedFormats, setSelectedFormats] = useState<
    DATASET_ASSET_FORMAT[]
  >([]);
  const [downloadLinks, setDownloadLinks] = useState<DownloadLinkType[]>([]);
  const [isDownloadLinkLoading, setIsDownloadLinkLoading] =
    useState<boolean>(false);

  const noDownloadLinksAvailable = useMemo(() => {
    return (
      selectedFormats.length > 0 &&
      selectedFormats.every((format) => {
        const link = downloadLinks.find((link) => link.filetype === format);

        return (
          !link || // no matching link
          link.downloadURL === undefined ||
          link.downloadURL === getNotAvailableText(format)
        );
      })
    );
  }, [selectedFormats, downloadLinks]);

  const isDownloadDisabled =
    !selectedFormats.length ||
    isDownloadLinkLoading ||
    isError ||
    isLoading ||
    noDownloadLinksAvailable;

  const availableFormats = useMemo(
    () => new Set(dataAssets.map((dataAsset) => dataAsset.filetype)),
    [dataAssets]
  );

  // Determine formats to download, ensuring ATAC_INDEX includes ATAC_FRAGMENT
  const formatsToDownload = useMemo(() => {
    return selectedFormats.flatMap((format) =>
      format === DATASET_ASSET_FORMAT.ATAC_INDEX
        ? [DATASET_ASSET_FORMAT.ATAC_INDEX, DATASET_ASSET_FORMAT.ATAC_FRAGMENT]
        : [format]
    );
  }, [selectedFormats]);

  // Set selected formats based on available formats and ATAC index
  useEffect(() => {
    if (!initialFetchDone) return;

    const isATACIncomplete = (format: DATASET_ASSET_FORMAT) => {
      if (format !== DATASET_ASSET_FORMAT.ATAC_INDEX) return false;
      const hasIndex = availableFormats.has(DATASET_ASSET_FORMAT.ATAC_INDEX);
      const hasFragment = availableFormats.has(
        DATASET_ASSET_FORMAT.ATAC_FRAGMENT
      );
      return !(hasIndex && hasFragment);
    };

    const allFormats = POSSIBLE_DOWNLOAD_FORMATS.filter(
      ({ format }) => availableFormats.has(format) && !isATACIncomplete(format)
    ).map(({ format }) => format);

    const defaultFormats: DATASET_ASSET_FORMAT[] = [];

    // Select H5AD by default if available
    allFormats.includes(DATASET_ASSET_FORMAT.H5AD)
      ? defaultFormats.push(DATASET_ASSET_FORMAT.H5AD)
      : defaultFormats.push(allFormats[0]);

    setSelectedFormats(defaultFormats);
  }, [availableFormats, initialFetchDone]);

  useEffect(() => {
    if (initialFetchDone || dataAssets.length === 0) return;
    const fetchDownloadLinks = async () => {
      setIsDownloadLinkLoading(true);

      try {
        const links = await Promise.all(
          dataAssets.map((asset) => getDownloadLink(asset))
        );
        const cleanedLinks = links
          .filter((link): link is DownloadLinkType => link !== null)
          .map((link) => ({
            ...link,
            downloadURL: link.downloadURL ?? getNotAvailableText(link.filetype),
          }));
        setDownloadLinks(cleanedLinks);
        setInitialFetchDone(true);
      } catch (error) {
        console.error("Error fetching download links", error);
      } finally {
        setIsDownloadLinkLoading(false);
      }
    };

    fetchDownloadLinks();
  }, [dataAssets, initialFetchDone]);
  /**
   * Tracks dataset download analytics as specified by the custom analytics event.
   * @param event - Custom analytics event.
   * @param dataFormat - Data format (optionally specified with event "DOWNLOAD_DATA_FORMAT_CLICKED").
   */
  const handleAnalytics = (
    event: EVENTS,
    dataFormat?: DATASET_ASSET_FORMAT
  ) => {
    track(event, {
      dataset_name: name,
      data_format: dataFormat || selectedFormats,
    });
  };

  const handleChange = (format: DATASET_ASSET_FORMAT) => {
    if (selectedFormats.includes(format)) {
      setSelectedFormats(selectedFormats.filter((f) => f !== format));
    } else {
      setSelectedFormats([...selectedFormats, format]);
    }
    handleAnalytics(EVENTS.DOWNLOAD_DATA_FORMAT_CLICKED, format);
  };

  const handleTabsChange = (_: SyntheticEvent, activeTabValue: ActiveTab) => {
    setActiveTab(activeTabValue);
  };

  return (
    <>
      <DialogTitle title="Download Data" />
      <StyledTabs
        value={activeTab}
        underlined
        sdsSize="small"
        onChange={handleTabsChange}
      >
        <Tab label="Browser" value={ActiveTab.Browser} />
        <Tab label="Python API" value={ActiveTab.PythonApi} />
        <Tab label="R API" value={ActiveTab.RApi} />
      </StyledTabs>
      <DialogContent>
        {activeTab === ActiveTab.Browser && (
          <BrowserTab
            name={name}
            isError={isError}
            isLoading={isLoading}
            availableFormats={availableFormats}
            handleChange={handleChange}
            downloadLinks={downloadLinks}
            formatsToDownload={formatsToDownload}
            handleAnalytics={handleAnalytics}
            selectedFormats={selectedFormats}
            isDownloadLinkLoading={isDownloadLinkLoading}
          />
        )}
        {activeTab === ActiveTab.PythonApi && (
          <ApiTab
            name={name}
            censusCopyText={`import cellxgene_census

                help(cellxgene_census)
                help(cellxgene_census.get_anndata)
            `}
            discoverCopyText={
              collectionId + `/datasets/${dataAssets[0].dataset_id}`
            }
          />
        )}
        {activeTab === ActiveTab.RApi && (
          <ApiTab
            name={name}
            censusCopyText={`library("cellxgene.census")

              ?cellxgene.census::get_seurat
            `}
            discoverCopyText={
              collectionId + `/datasets/${dataAssets[0].dataset_id}`
            }
          />
        )}
      </DialogContent>
      <DialogActions>
        <Button
          isAllCaps={false}
          onClick={onClose}
          sdsStyle="minimal"
          sdsType="secondary"
        >
          Cancel
        </Button>
        <Tooltip
          placement="top"
          disableHoverListener={!isDownloadDisabled}
          disableFocusListener={!isDownloadDisabled}
          disableTouchListener={!isDownloadDisabled}
          title={"Select at least one valid format"}
        >
          <span>
            <Button
              data-testid="download-asset-download-button"
              disabled={isDownloadDisabled}
              onClick={() => {
                downloadMultipleFiles(formatsToDownload, downloadLinks);
                handleAnalytics(EVENTS.DOWNLOAD_DATA_COMPLETE);
              }}
              sdsStyle="square"
              sdsType="primary"
            >
              Download
            </Button>
          </span>
        </Tooltip>
      </DialogActions>
    </>
  );
};

export default Content;
