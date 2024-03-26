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
  // CopyButtonWrapper,
  RunButton,
  RunButtonWrapper,
} from "./style";
import QueryGroupFilters from "./components/Filters";
import Organism from "./components/Organism";
import {
  // copyCellGroup1,
  submitQueryGroups,
} from "src/views/DifferentialExpression/common/store/actions";
import DeResults from "./components/DeResults";
import Loader from "./components/Loader";

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
  // const {
  //   queryGroup1: queryGroupWithNames1,
  //   queryGroup2: queryGroupWithNames2,
  // } = queryGroupsWithNames;

  // const handleCopyCellGroup1 = () => {
  //   if (!dispatch) return;
  //   dispatch(copyCellGroup1());
  // };
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

  return (
    <Wrapper>
      <div
        style={{
          width: "100%",
          display: "flex",
          flexDirection: "row",
          columnGap: "120px",
        }}
      >
        {isLoading && <Loader />}
        <div style={{ width: "655px" }}>
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

          <Organism />
          <div style={{ flexDirection: "row", display: "flex" }}>
            <div>
              <CellGroupTitle>Cell Group 1</CellGroupTitle>
              <QueryGroupFilters
                key={`query-group-1`}
                queryGroup={queryGroup1}
                isQueryGroup1={true}
              />
            </div>
            <div>COPY AND INVERT BUTTONS</div>
            <div>
              <CellGroupTitle>
                Cell Group 2
                {/* <CopyButtonWrapper onClick={handleCopyCellGroup1}>
                  Copy Cell Group 1
                </CopyButtonWrapper> */}
              </CellGroupTitle>
              <QueryGroupFilters
                key={`query-group-2`}
                queryGroup={queryGroup2}
                isQueryGroup1={false}
              />
            </div>
          </div>
          <RunButtonWrapper>
            <RunButton
              color="primary"
              size="large"
              variant="contained"
              onClick={handleRunDifferentialExpression}
              disabled={!canRunDifferentialExpression}
            >
              Run differential expression
            </RunButton>
          </RunButtonWrapper>
        </div>
        <DeResults setIsLoading={setIsLoading} />
      </div>
    </Wrapper>
  );
}
