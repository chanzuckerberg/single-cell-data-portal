import { useContext, useEffect, useState } from "react";
import TextField from "@mui/material/TextField";
import InputAdornment from "@mui/material/InputAdornment";
import ArrowForwardIosIcon from "@mui/icons-material/ArrowForwardIos";
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
} from "./style";
import QueryGroupFilters from "./components/Filters";
import Organism from "./components/Organism";
import {
  copyCellGroup1,
  selectOrganism,
  setQueryGroupFilters,
  submitQueryGroups,
} from "src/views/DifferentialExpression/common/store/actions";
import DeResults from "./components/DeResults";
import Loader from "./components/Loader";
import { useNaturalLanguageDeQuery } from "src/common/queries/differentialExpression";
import { CircularProgress } from "@mui/material";

export default function DifferentialExpression(): JSX.Element {
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [isLoadingGetDeQuery, setIsLoadingGetDeQuery] =
    useState<boolean>(false);
  const [inputText, setInputText] = useState<string>("");
  const [tempInput, setTempInput] = useState<string>("");

  useEffect(() => {
    setIsLoadingGetDeQuery(isLoadingGetDeQuery);
  }, [isLoadingGetDeQuery]);
  const { organismId } = useContext(StateContext);
  const dispatch = useContext(DispatchContext);
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

  const handleInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setTempInput(event.target.value);
  };

  const handleKeyDown = (event: React.KeyboardEvent) => {
    if (event.key === "Enter") {
      setInputText(tempInput);
    }
  };

  const {
    organism,
    queryCriteria1,
    queryCriteria2,
    queryCriteriaNames1,
    queryCriteriaNames2,
    isLoading: isLoadingGetDeQueryRaw,
  } = useNaturalLanguageDeQuery(inputText);

  useEffect(() => {
    setIsLoadingGetDeQuery(isLoadingGetDeQueryRaw);
  }, [isLoadingGetDeQueryRaw]);

  useEffect(() => {
    if (!dispatch || isLoadingGetDeQuery || organism === "") return;
    if (organismId !== organism) dispatch(selectOrganism(organism));
    dispatch(
      setQueryGroupFilters({
        queryGroup1: queryCriteria1,
        queryGroup2: queryCriteria2,
        queryGroupNames1: queryCriteriaNames1,
        queryGroupNames2: queryCriteriaNames2,
      })
    );
  }, [
    organism,
    organismId,
    queryCriteria1,
    queryCriteria2,
    queryCriteriaNames1,
    queryCriteriaNames2,
    isLoadingGetDeQuery,
    dispatch,
  ]);

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
            Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do
            eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim
            ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut
            aliquip ex ea commodo consequat.
          </StepSubHeader>
          <TextField
            label="Ex: Find me differentially expressed genes between plasma cells in lung and plasma cells in blood"
            variant="outlined"
            fullWidth
            onChange={handleInputChange}
            onKeyDown={handleKeyDown}
            sx={{ height: "60px" }}
            InputProps={{
              endAdornment: (
                <InputAdornment position="end">
                  {isLoadingGetDeQuery ? (
                    <CircularProgress size={20} />
                  ) : (
                    <ArrowForwardIosIcon />
                  )}
                </InputAdornment>
              ),
            }}
          />
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
        </div>
        <DeResults setIsLoading={setIsLoading} />
      </div>
    </Wrapper>
  );
}
