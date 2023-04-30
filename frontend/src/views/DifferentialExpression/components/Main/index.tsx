import { useContext, useState } from "react";
import { DrawerSize } from "@blueprintjs/core";
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
  CopyButtonWrapper,
  RunButton,
  RunButtonWrapper,
  StyledSidebarDrawer,
} from "./style";
import QueryGroupFilters from "./components/Filters";
import Organism from "./components/Organism";
import { copyCellGroup1 } from "src/views/DifferentialExpression/common/store/actions";
import DeResults from "./components/DeResults";

export default function DifferentialExpression(): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const [isSidebarOpen, setIsSidebarOpen] = useState(false);
  const { queryGroups, queryGroupsWithNames } = useContext(StateContext);
  const { queryGroup1, queryGroup2 } = queryGroups;
  const {
    queryGroup1: queryGroupWithNames1,
    queryGroup2: queryGroupWithNames2,
  } = queryGroupsWithNames;

  const handleCopyCellGroup1 = () => {
    if (!dispatch) return;
    dispatch(copyCellGroup1());
  };

  const canRunDifferentialExpression = true;

  const handleRunDifferentialExpression = () => {
    setIsSidebarOpen(true);
  };
  const handleCloseSidebar = () => {
    setIsSidebarOpen(false);
  };

  return (
    <Wrapper>
      <StepHeader>
        <WordPop>Differential</WordPop> Expression
      </StepHeader>
      <StepSubHeader>
        Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod
        tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim
        veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea
        commodo consequat.
      </StepSubHeader>
      <Organism />
      <CellGroupTitle>Cell Group 1</CellGroupTitle>
      <QueryGroupFilters
        key={`query-group-1`}
        queryGroup={queryGroup1}
        queryGroupWithNames={queryGroupWithNames1}
        isQueryGroup1={true}
      />
      <CellGroupTitle>
        Cell Group 2
        <CopyButtonWrapper onClick={handleCopyCellGroup1}>
          Copy Cell Group 1
        </CopyButtonWrapper>
      </CellGroupTitle>
      <QueryGroupFilters
        key={`query-group-2`}
        queryGroup={queryGroup2}
        queryGroupWithNames={queryGroupWithNames2}
        isQueryGroup1={false}
      />
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
      <StyledSidebarDrawer
        position="right"
        isOpen={isSidebarOpen}
        title="Source Data"
        canEscapeKeyClose={true}
        canOutsideClickClose={true}
        onClose={handleCloseSidebar}
        size={DrawerSize.SMALL}
      >
        <DeResults />
      </StyledSidebarDrawer>
    </Wrapper>
  );
}
