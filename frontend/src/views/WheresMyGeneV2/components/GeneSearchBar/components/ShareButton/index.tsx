import { Notification } from "@czi-sds/components";
import { useCallback, useContext, useEffect, useMemo, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EMPTY_OBJECT, noop } from "src/common/constants/utils";
import { isSSR } from "src/common/utils/isSSR";
import { useRouter } from "next/router";

import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../SaveExport";
import {
  StyledButtonDiv,
  StyledButtonIcon,
  StyledIcon,
  StyledLabel,
  StyledNotificationDetails,
  StyledNotificationLabel,
  StyledNotificationWrapper,
} from "./style";
import { generateAndCopyShareUrl, loadStateFromQueryParams } from "./utils";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import {
  CellTypeByTissueName,
  useCellTypesByTissueName,
  usePrimaryFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { getCompareOptionNameById } from "src/views/WheresMyGene/common/constants";
import { CellType } from "src/views/WheresMyGene/common/types";
import { useTissueMetadata } from "src/common/queries/cellGuide";

export default function ShareButton(): JSX.Element {
  const state = useContext(StateContext);
  const { data: tissues } = useTissueMetadata();
  const router = useRouter();

  const {
    selectedFilters,
    selectedGenes,
    selectedOrganismId,
    compare,
    filteredCellTypes,
  } = state;

  // (seve): this is done in multiple places, we should consolidate to a single hook

  const {
    data: rawCellTypesByTissueName,
    isLoading: isLoadingCellTypesByTissueName,
  } = useCellTypesByTissueName();

  const [cellTypesByTissueName, setCellTypesByTissueName] =
    useState<CellTypeByTissueName>(EMPTY_OBJECT);

  // This is needed to prevent overwriting the cellTypesByTissueName state with empty
  useEffect(() => {
    if (isLoadingCellTypesByTissueName) return;

    setCellTypesByTissueName(rawCellTypesByTissueName);
  }, [rawCellTypesByTissueName, isLoadingCellTypesByTissueName]);

  const cellTypesByName = useMemo(() => {
    const result: { [name: string]: CellType } = {};
    Object.values(cellTypesByTissueName).forEach((cellTypes) => {
      cellTypes.forEach((cellType) => {
        result[cellType.cellTypeName] = cellType;
      });
    });
    return result;
  }, [cellTypesByTissueName]);

  const mapCellTypesToIDs = useCallback(
    (cellType: string) => cellTypesByName[cellType].id,
    [cellTypesByName]
  );

  const { isLoading: isLoadingFilterDims } = usePrimaryFilterDimensions(); //temp explicit version
  const dispatch = useContext(DispatchContext);
  const [showURLCopyNotification, setShowURLCopyNotification] = useState(0);

  const copyShareUrl = useCallback(() => {
    if (!dispatch) return;

    generateAndCopyShareUrl({
      compare,
      filters: selectedFilters,
      organism: selectedOrganismId,
      genes: selectedGenes,
      cellTypes: filteredCellTypes,
    });

    const filteredCellTypeIDs = filteredCellTypes.map(mapCellTypesToIDs);

    track(EVENTS.WMG_SHARE_CLICKED, {
      organism: selectedOrganismId,
      dataset_filter: selectedFilters.datasets,
      disease_filter: selectedFilters.diseases,
      genes: selectedGenes,
      group_by_option: getCompareOptionNameById(compare),
      self_reported_ethnicity_filter: selectedFilters.ethnicities,
      publication_filter: selectedFilters.publications,
      sex_filter: selectedFilters.sexes,
      tissue_filter: selectedFilters.tissues,
      cell_types_selected: filteredCellTypeIDs,
    });

    setShowURLCopyNotification((prev) => prev + 1);
  }, [
    selectedFilters,
    selectedGenes,
    selectedOrganismId,
    dispatch,
    compare,
    filteredCellTypes,
    mapCellTypesToIDs,
  ]);

  useEffect(() => {
    if (
      isSSR() ||
      isLoadingFilterDims ||
      !dispatch ||
      Object.keys(cellTypesByName).length === 0
    )
      return;
    const { search } = window.location;
    const params = new URLSearchParams(search);
    if (params) {
      // If we later want to display a toast on successful load from url, this function returns true/false
      const loadedState = loadStateFromQueryParams({
        cellTypesByName,
        dispatch,
        params,
        router,
        selectedFilters,
        tissues,
      });

      if (loadedState) {
        const filteredCellTypeIDs =
          loadedState.cellTypes?.map(mapCellTypesToIDs) || [];

        track(EVENTS.WMG_SHARE_LOADED, {
          genes: loadedState.genes,
          organism: loadedState.organism,
          dataset_filter: loadedState.filters.datasets,
          disease_filter: loadedState.filters.diseases,
          group_by_option: getCompareOptionNameById(loadedState.compare),
          self_reported_ethnicity_filter: loadedState.filters.ethnicities,
          publication_filter: loadedState.filters.publications,
          sex_filter: loadedState.filters.sexes,
          cell_types_selected: filteredCellTypeIDs,
        });
      }
    }
  }, [
    isLoadingFilterDims,
    dispatch,
    selectedFilters,
    cellTypesByName,
    compare,
    mapCellTypesToIDs,
    tissues,
    router,
  ]);

  return (
    <>
      {showURLCopyNotification > 0 && (
        <StyledNotificationWrapper>
          <Notification
            key={showURLCopyNotification}
            autoDismiss={5000}
            onClose={noop}
            slideDirection={"left"}
            intent={"info"}
            icon={
              <StyledIcon sdsIcon={"link"} sdsSize={"s"} sdsType={"static"} />
            }
          >
            <StyledNotificationLabel data-testid="share-link-notification">
              Share link copied
            </StyledNotificationLabel>
            <StyledNotificationDetails>
              We regularly expand our single cell data corpus to improve
              results. Downloaded data and figures may differ in the future.
            </StyledNotificationDetails>
          </Notification>
        </StyledNotificationWrapper>
      )}

      <StyledButtonDiv className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
        <StyledLabel>Share</StyledLabel>

        <StyledButtonIcon
          data-testid={"share-button"}
          onClick={copyShareUrl}
          sdsSize="medium"
          sdsType="primary"
          sdsIcon="share"
          disabled={selectedGenes.length === 0}
        />
      </StyledButtonDiv>
    </>
  );
}
