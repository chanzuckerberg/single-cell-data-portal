import { Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { useCallback, useContext, useEffect } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { usePrimaryFilterDimensions } from "src/common/queries/wheresMyGene";
import { isSSR } from "src/common/utils/isSSR";
import Toast from "src/views/Collection/components/Toast";
import { getCompareOptionNameById } from "src/views/WheresMyGene/common/constants";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { StyledButtonIcon } from "../QuickSelect/style";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../SaveImage";
import { StyledButtonDiv, StyledLabel } from "./style";
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
  // const [showURLCopyNotification, setShowURLCopyNotification] = useState(0);

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
      tissues: selectedTissues,
      genes: selectedGenes,
    });

    track(EVENTS.WMG_SHARE_CLICKED, {
      dataset_filter: selectedFilters.datasets,
      disease_filter: selectedFilters.diseases,
      genes: selectedGenes,
      group_by_option: getCompareOptionNameById(compare),
      self_reported_ethnicity_filter: selectedFilters.ethnicities,
      sex_filter: selectedFilters.sexes,
      tissues: selectedTissues,
    });

    // setShowURLCopyNotification((prev) => prev + 1);
    Toast.show({
      icon: IconNames.LINK,
      intent: Intent.PRIMARY,
      message: "Share link copied",
      timeout: 1000,
    });
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
          sex_filter: loadedState.filters.sexes,
        });
      }
    }
  }, [isLoadingFilterDims, dispatch, selectedFilters]);

  return (
    <>
      {/* {showURLCopyNotification > 0 && (
        <StyledNotification
          key={showURLCopyNotification}
          autoDismiss={false}
          dismissDirection={"left"}
          intent={"info"}
        >
          Share link copied
        </StyledNotification>
      )} */}

      <StyledButtonDiv className={EXCLUDE_IN_SCREENSHOT_CLASS_NAME}>
        <StyledLabel>Share</StyledLabel>

        <StyledButtonIcon
          data-test-id={"share-button"}
          onClick={copyShareUrl}
          sdsSize="medium"
          sdsType="primary"
          sdsIcon="share"
          disabled={selectedTissues.length === 0 || selectedGenes.length === 0}
        />
      </StyledButtonDiv>
    </>
  );
}
