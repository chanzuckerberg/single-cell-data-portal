import { Notification } from "@czi-sds/components";
import { useCallback, useContext, useEffect, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { noop } from "src/common/constants/utils";
import { usePrimaryFilterDimensions } from "src/common/queries/wheresMyGene";
import { isSSR } from "src/common/utils/isSSR";
import { getCompareOptionNameById } from "src/views/WheresMyGene/common/constants";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { StyledButtonIcon } from "../QuickSelect/style";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../SaveExport";
import {
  StyledButtonDiv,
  StyledIcon,
  StyledLabel,
  StyledNotificationDetails,
  StyledNotificationLabel,
  StyledNotificationWrapper,
} from "./style";
import { generateAndCopyShareUrl, loadStateFromQueryParams } from "./utils";

export default function ShareButton(): JSX.Element {
  const state = useContext(StateContext);

  const {
    selectedFilters,
    selectedTissues,
    selectedGenes,
    selectedOrganismId,
    compare,
  } = state;

  const { isLoading: isLoadingFilterDims } = usePrimaryFilterDimensions();
  const dispatch = useContext(DispatchContext);
  const [showURLCopyNotification, setShowURLCopyNotification] = useState(0);

  const copyShareUrl = useCallback(() => {
    if (!dispatch) return;

    generateAndCopyShareUrl({
      compare,
      filters: selectedFilters,
      organism: selectedOrganismId,
      tissues: selectedTissues,
      genes: selectedGenes,
    });

    track(EVENTS.WMG_SHARE_CLICKED, {
      organism: selectedOrganismId,
      dataset_filter: selectedFilters.datasets,
      disease_filter: selectedFilters.diseases,
      genes: selectedGenes,
      group_by_option: getCompareOptionNameById(compare),
      self_reported_ethnicity_filter: selectedFilters.ethnicities,
      publication_filter: selectedFilters.publications,
      sex_filter: selectedFilters.sexes,
      tissues: selectedTissues,
    });

    setShowURLCopyNotification((prev) => prev + 1);
  }, [
    selectedFilters,
    selectedTissues,
    selectedGenes,
    selectedOrganismId,
    dispatch,
    compare,
  ]);

  useEffect(() => {
    if (isSSR() || isLoadingFilterDims || !dispatch) return;
    const { search } = window.location;
    const params = new URLSearchParams(search);
    if (params) {
      // If we later want to display a toast on successful load from url, this function returns true/false
      const loadedState = loadStateFromQueryParams(
        params,
        selectedFilters,
        dispatch
      );

      if (loadedState) {
        track(EVENTS.WMG_SHARE_LOADED, {
          tissues: loadedState.tissues,
          genes: loadedState.genes,
          organism: loadedState.organism,
          dataset_filter: loadedState.filters.datasets,
          disease_filter: loadedState.filters.diseases,
          group_by_option: getCompareOptionNameById(loadedState.compare),
          self_reported_ethnicity_filter: loadedState.filters.ethnicities,
          publication_filter: loadedState.filters.publications,
          sex_filter: loadedState.filters.sexes,
        });
      }
    }
  }, [isLoadingFilterDims, dispatch, selectedFilters]);

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
          disabled={selectedTissues?.length === 0 || selectedGenes.length === 0}
        />
      </StyledButtonDiv>
    </>
  );
}
