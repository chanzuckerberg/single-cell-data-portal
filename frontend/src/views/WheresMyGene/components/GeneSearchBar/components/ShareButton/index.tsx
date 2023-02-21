import { Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { Tooltip } from "czifui";
import { useCallback, useContext, useEffect, useState } from "react";
import { usePrimaryFilterDimensions } from "src/common/queries/wheresMyGene";
import { isSSR } from "src/common/utils/isSSR";
import Toast from "src/views/Collection/components/Toast";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { StyledButtonIcon } from "../QuickSelect/style";
import { EXCLUDE_IN_SCREENSHOT_CLASS_NAME } from "../SaveImage";
import { StyledButtonDiv, StyledLabel, StyledNotification } from "./style";
import { generateAndCopyShareUrl, loadStateFromQueryParams } from "./utils";

export default function ShareButton(): JSX.Element {
  const state = useContext(StateContext);
  const { selectedFilters, selectedTissues, selectedGenes } = state;
  const { isLoading: isLoadingFilterDims } = usePrimaryFilterDimensions();
  const dispatch = useContext(DispatchContext);
  // const [showURLCopyNotification, setShowURLCopyNotification] = useState(0);

  const copyShareUrl = useCallback(() => {
    if (!dispatch) return;
    generateAndCopyShareUrl(selectedFilters, selectedTissues, selectedGenes);
    // setShowURLCopyNotification((prev) => prev + 1);
    Toast.show({
      icon: IconNames.LINK,
      intent: Intent.PRIMARY,
      message: "Share link copied",
      timeout: 8000,
    });
  }, [selectedFilters, selectedTissues, selectedGenes, dispatch]);

  useEffect(() => {
    if (isSSR() || isLoadingFilterDims || !dispatch) return;
    const { search } = window.location;
    const params = new URLSearchParams(search);
    if (params) {
      // If we later want to display a toast on successful load from url, this function returns true/false
      loadStateFromQueryParams(params, selectedFilters, dispatch);
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
        />
      </StyledButtonDiv>
    </>
  );
}
