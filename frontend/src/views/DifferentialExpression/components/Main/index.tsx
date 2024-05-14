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
} from "./style";
import QueryGroupFilters from "./components/Filters";
import Organism from "./components/Organism";
import { submitQueryGroups } from "src/views/DifferentialExpression/common/store/actions";
import DeResults from "./components/DeResults";
import Loader from "./components/Loader";
import useProcessedQueryGroupFilterDimensions from "./components/common/query_group_filter_dimensions";
import Method from "./components/Method";

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
  const isQueryGroup2NotEmpty = Object.values(queryGroup1).some(
    (value) => value.length > 0
  );
  const canRunDifferentialExpression =
    !isLoading && isQueryGroup1NotEmpty && isQueryGroup2NotEmpty;

  const handleRunDifferentialExpression = () => {
    if (!dispatch) return;
    dispatch(submitQueryGroups());
  };

  const { n_cells: nCellsGroup1 } =
    useProcessedQueryGroupFilterDimensions(queryGroup1);
  const { n_cells: nCellsGroup2 } =
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
            <StepSubHeader>
              Find differentially expressed genes between custom group of cells
              across the CELLxGENE data corpus. For additional help and
              information, read our documentation.
              <br />
              <br />
              [Method Information]
            </StepSubHeader>
            <FlexRow>
              <Organism />
              <Method />
            </FlexRow>
            <FlexRow>
              <div>
                <CellGroupTitle>
                  Cell Group 1
                  <CellCountTitle>
                    {" "}
                    {nCellsGroup1.toLocaleString()} cells
                  </CellCountTitle>
                </CellGroupTitle>
                <QueryGroupFilters
                  key={`query-group-1`}
                  queryGroup={queryGroup1}
                  isQueryGroup1={true}
                />
              </div>

              <div>
                <CellGroupTitle>
                  Cell Group 2
                  <CellCountTitle>
                    {" "}
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
