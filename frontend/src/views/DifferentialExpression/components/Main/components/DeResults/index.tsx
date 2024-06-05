import {
  StyledIcon,
  ButtonsWrapper,
  ButtonLabel,
  InstructionsBody,
  InstructionsHeader,
  InstructionsWrapper,
  ResultsHeader,
  ResultsWrapper,
  ResultsHeaderWrapper,
  FlexRow,
} from "./style";

import { StyledSidebarDrawer } from "src/views/WheresMyGeneV2/components/Main/style";
import { DrawerSize } from "@blueprintjs/core";
import SourceData from "./components/SourceData";
import DifferentialExpressionResults from "./components/DifferentialExpressionResults";
import {
  DIFFERENTIAL_EXPRESSION_INSTRUCTIONS_SIDEBAR,
  DIFFERENTIAL_EXPRESSION_RESULTS_DOWNLOAD_BUTTON,
  DIFFERENTIAL_EXPRESSION_SOURCE_DATA_BUTTON,
} from "src/views/DifferentialExpression/common/constants";
import { useConnect } from "./connect";
import { Props } from "./types";

export default function DeResults({ setIsLoading }: Props): JSX.Element {
  const {
    queryGroups,
    isLoading,
    queryGroupsWithNames,
    organismId,
    setSearchQuery,
    setLfcFilter,
    setEffectSizeFilter,
    sortDirection,
    setSortDirection,
    isSourceDatasetSidebarOpen,
    setIsSourceDatasetSidebarOpen,
    sortedAndFilteredResults,
    downloadCSV,
    showEmpty,
    nOverlap,
  } = useConnect({ setIsLoading });

  return (
    <div>
      {!showEmpty ? (
        <ResultsWrapper>
          <ResultsHeaderWrapper>
            <ResultsHeader>Results</ResultsHeader>
            <FlexRow>
              <ButtonsWrapper
                onClick={downloadCSV}
                data-testid={DIFFERENTIAL_EXPRESSION_RESULTS_DOWNLOAD_BUTTON}
                disabled={isLoading}
              >
                <StyledIcon sdsIcon="Download" sdsSize="l" sdsType="static" />
                <ButtonLabel>Download</ButtonLabel>
              </ButtonsWrapper>
              <ButtonsWrapper
                onClick={() =>
                  !isLoading && setIsSourceDatasetSidebarOpen(true)
                }
                data-testid={DIFFERENTIAL_EXPRESSION_SOURCE_DATA_BUTTON}
                disabled={isLoading}
              >
                <StyledIcon sdsIcon="InfoCircle" sdsSize="l" sdsType="static" />
                <ButtonLabel>Source Data</ButtonLabel>
              </ButtonsWrapper>
            </FlexRow>
          </ResultsHeaderWrapper>
          {!!queryGroups &&
            !!queryGroupsWithNames &&
            !!organismId &&
            !isLoading && (
              <DifferentialExpressionResults
                queryGroups={queryGroups}
                queryGroupsWithNames={queryGroupsWithNames}
                organismId={organismId}
                sortedAndFilteredResults={sortedAndFilteredResults}
                nCellsOverlap={nOverlap}
                setSearchQuery={setSearchQuery}
                setLfcFilter={setLfcFilter}
                setEffectSizeFilter={setEffectSizeFilter}
                sortDirection={sortDirection}
                setSortDirection={setSortDirection}
              />
            )}
        </ResultsWrapper>
      ) : (
        <InstructionsWrapper
          data-testid={DIFFERENTIAL_EXPRESSION_INSTRUCTIONS_SIDEBAR}
        >
          <InstructionsHeader>Instructions</InstructionsHeader>
          <InstructionsBody>
            <ol>
              <li>
                Select a cell group of interest within the Cell Group 1 box by
                using the dropdown selectors.
                <br />
                <br />
                To copy the same selection over to Cell Group 2, click the copy
                button to the right of each dropdown in Cell Group 1.
              </li>
              <li>
                Within Cell Group 2, select a group that the cell group of
                interest will be compared to.
              </li>
            </ol>
          </InstructionsBody>
        </InstructionsWrapper>
      )}

      <StyledSidebarDrawer
        position="right"
        isOpen={isSourceDatasetSidebarOpen}
        title="Source Data"
        canEscapeKeyClose={true}
        canOutsideClickClose={true}
        onClose={() => setIsSourceDatasetSidebarOpen(false)}
        size={DrawerSize.SMALL}
      >
        <SourceData />
      </StyledSidebarDrawer>
    </div>
  );
}
