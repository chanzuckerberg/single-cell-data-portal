import { useContext, useEffect, useState } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
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
import { submitQueryGroups } from "src/views/DifferentialExpression/common/store/actions";
import DeResults from "./components/DeResults";
import Loader from "./components/Loader";
import useProcessedQueryGroupFilterDimensions from "./components/common/query_group_filter_dimensions";
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
import { craftPayloadWithQueryGroups } from "./utils";

export default function DifferentialExpression(): JSX.Element {
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [isLoadingGetDeQuery, setIsLoadingGetDeQuery] =
    useState<boolean>(false);

  useEffect(() => {
    setIsLoadingGetDeQuery(isLoadingGetDeQuery);
  }, [isLoadingGetDeQuery]);
  const dispatch = useContext(DispatchContext);
  const { queryGroups } = useContext(StateContext);
  const { queryGroup1, queryGroup2 } = queryGroups;

  // check if any values in queryGroup1 are not empty
  const isQueryGroup1NotEmpty = Object.values(queryGroup1).some(
    (value) => value.length > 0
  );
  const isQueryGroup2NotEmpty = Object.values(queryGroup2).some(
    (value) => value.length > 0
  );
  const canRunDifferentialExpression =
    !isLoading && isQueryGroup1NotEmpty && isQueryGroup2NotEmpty;

  const handleRunDifferentialExpression = () => {
    if (!dispatch) return;
    dispatch(submitQueryGroups());

    track(
      EVENTS.DE_FIND_GENES_CLICKED,
      craftPayloadWithQueryGroups(queryGroups)
    );
  };

  const { n_cells: nCellsGroup1, isLoading: isLoadingGroup1 } =
    useProcessedQueryGroupFilterDimensions(queryGroup1);
  const { n_cells: nCellsGroup2, isLoading: isLoadingGroup2 } =
    useProcessedQueryGroupFilterDimensions(queryGroup2);

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
              information, read our documentation.
              <br />
              <br />
              This tool uses Welch&apos;s t-test to identify differentially
              expressed genes between groups of cell in the CELLxGENE Census.
              While the t-test performs reasonably well on individual datasets{" "}
              <a
                href="https://www.mdpi.com/2073-4425/12/12/1947"
                target="_blank"
                rel="noopener noreferrer"
              >
                [1]
              </a>
              <a
                href="https://www.nature.com/articles/s41467-019-12266-7"
                target="_blank"
                rel="noopener noreferrer"
              >
                [2]
              </a>
              <a
                href="https://www.nature.com/articles/nmeth.4612"
                target="_blank"
                rel="noopener noreferrer"
              >
                [3]
              </a>
              , its performance on concatenated, non-integrated datasets has not
              been extensively evaluated. We recommend using this tool for
              preliminary investigations and following up with a more robust
              method for formal analysis. Learn more about our data filtering
              and normalization{" "}
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
