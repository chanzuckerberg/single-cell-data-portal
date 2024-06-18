import {
  Wrapper,
  StepHeader,
  WordPop,
  StepSubHeader,
  CellGroupTitle,
  RunButton,
  RunButtonWrapper,
  FlexRow,
  QuerySelectorWrapper,
  TwoPanelLayout,
  CellCountTitle,
  Spinner,
} from "./style";
import QueryGroupFilters from "./components/Filters";
import Organism from "./components/Organism";
import DeResults from "./components/DeResults";
import Loader from "./components/Loader";
import Method from "./components/Method";
import {
  DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER,
  DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER,
  DIFFERENTIAL_EXPRESSION_FILTERS_LOADING_SPINNER,
  DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT,
  DIFFERENTIAL_EXPRESSION_FIND_GENES_BUTTON,
  DIFFERENTIAL_EXPRESSION_METHOD_INFO_TEXT,
} from "../../common/constants";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { useConnect } from "./connect";

export default function DifferentialExpression(): JSX.Element {
  const {
    isLoading,
    setIsLoading,
    queryGroup1,
    queryGroup2,
    canRunDifferentialExpression,
    handleRunDifferentialExpression,
    nCellsGroup1,
    isLoadingGroup1,
    nCellsGroup2,
    isLoadingGroup2,
  } = useConnect();

  return (
    <TwoPanelLayout>
      <div className="leftPanel">
        <Wrapper>
          {isLoading && <Loader />}

          <QuerySelectorWrapper>
            <StepHeader>
              <WordPop>Differential</WordPop> Expression
            </StepHeader>
            <StepSubHeader
              data-testid={DIFFERENTIAL_EXPRESSION_METHOD_INFO_TEXT}
            >
              Find differentially expressed genes between custom group of cells
              across the CELLxGENE data corpus. For additional help and
              information, read our{" "}
              <a
                href="/docs/01__CellxGene"
                target="_blank"
                rel="noopener noreferrer"
                onClick={() => track(EVENTS.DE_DOCUMENTATION_CLICKED)}
              >
                documentation
              </a>
              .
              <br />
              <br />
              This tool uses Welch&apos;s t-test to identify differentially
              expressed genes between groups of cells in the CELLxGENE Census.
              We recommend using this tool for preliminary investigations and
              following up with more robust methods for formal analysis. Learn
              more about data filtering and normalization{" "}
              <a
                href="/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_3__Gene%20Expression%20Data%20Processing"
                target="_blank"
                rel="noopener noreferrer"
                onClick={() => track(EVENTS.DE_DOCUMENTATION_CLICKED)}
              >
                here
              </a>
              .
            </StepSubHeader>
            <FlexRow>
              <Organism />
              <Method />
            </FlexRow>
            <FlexRow>
              <div data-testid={DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER}>
                <CellGroupTitle>
                  Cell Group 1
                  <CellCountTitle
                    data-testid={DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT}
                  >
                    {isLoadingGroup1 && (
                      <Spinner
                        data-testid={
                          DIFFERENTIAL_EXPRESSION_FILTERS_LOADING_SPINNER
                        }
                      />
                    )}
                    {nCellsGroup1.toLocaleString()} cells
                  </CellCountTitle>
                </CellGroupTitle>
                <QueryGroupFilters
                  key={`query-group-1`}
                  queryGroup={queryGroup1}
                  isQueryGroup1={true}
                />
              </div>

              <div data-testid={DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER}>
                <CellGroupTitle>
                  Cell Group 2
                  <CellCountTitle
                    data-testid={DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT}
                  >
                    {isLoadingGroup2 && (
                      <Spinner
                        data-testid={
                          DIFFERENTIAL_EXPRESSION_FILTERS_LOADING_SPINNER
                        }
                      />
                    )}
                    {nCellsGroup2.toLocaleString()} cells
                  </CellCountTitle>
                </CellGroupTitle>
                <QueryGroupFilters
                  key={`query-group-2`}
                  queryGroup={queryGroup2}
                  isQueryGroup1={false}
                />
              </div>
            </FlexRow>
            <RunButtonWrapper>
              <RunButton
                color="primary"
                size="large"
                variant="contained"
                onClick={handleRunDifferentialExpression}
                disabled={!canRunDifferentialExpression}
                data-testid={DIFFERENTIAL_EXPRESSION_FIND_GENES_BUTTON}
              >
                Find genes
              </RunButton>
            </RunButtonWrapper>
          </QuerySelectorWrapper>
        </Wrapper>
      </div>
      <div className="rightPanel">
        <DeResults setIsLoading={setIsLoading} />
      </div>
    </TwoPanelLayout>
  );
}
