import React, { useContext } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import { StepTitle } from "../../../Steps/style";
import { StepHeader, StepSubHeader, WordPop } from "../StepOne/style";
import {
  AddQueryGroupButton,
  QueryGroupFiltersWrapper,
  ButtonsWrapper,
  BackButton,
  NextButton,
} from "./style";
import { Icon } from "czifui";
import { addQueryGroup } from "src/views/DifferentialExpression/common/store/actions";
import QueryGroupFilters from "./components/Filters";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";

interface Props {
  setStep: (step: number) => void;
}
export default function StepTwo({ setStep }: Props): JSX.Element {
  const state = useContext(StateContext);
  const dispatch = useContext(DispatchContext);

  const { queryGroups, queryGroupsWithNames } = state;

  const handleGoNext = () => {
    setStep(3);
  };
  const handleGoBack = () => {
    setStep(1);
  };
  const handleAddQueryGroup = () => {
    if (!dispatch) return null;
    dispatch(addQueryGroup());
  };
  let canGoNext;
  for (const queryGroup of queryGroups ?? []) {
    canGoNext = false;
    for (const key in queryGroup) {
      if (queryGroup[key as keyof QueryGroup].length > 0) {
        canGoNext = true;
        break;
      }
    }
    if (!canGoNext) break;
  }
  return (
    <div>
      <StepTitle> Step 2 </StepTitle>
      <StepHeader>
        <WordPop>What</WordPop> do you want to find differentially expressed
        genes for?
      </StepHeader>
      <StepSubHeader>
        Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod
        tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim
        veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea
        commodo consequat. Duis aute irure dolor in reprehenderit in voluptate
        velit esse cillum dolore eu fugiat nulla pariatur.
      </StepSubHeader>
      <AddQueryGroupButton
        color="primary"
        size="medium"
        startIcon={<Icon sdsIcon="plus" sdsSize="s" sdsType="button" />}
        onClick={handleAddQueryGroup}
      >
        Query group
      </AddQueryGroupButton>
      <QueryGroupFiltersWrapper>
        {queryGroupsWithNames &&
          queryGroups?.map((queryGroup, index) => {
            return (
              <QueryGroupFilters
                key={`query-group-${index}`}
                queryGroupIndex={index}
                queryGroup={queryGroup}
                queryGroupWithNames={queryGroupsWithNames[index]}
              />
            );
          })}
      </QueryGroupFiltersWrapper>
      <ButtonsWrapper>
        <BackButton
          color="primary"
          size="large"
          variant="contained"
          onClick={handleGoBack}
        >
          Back
        </BackButton>
        <NextButton
          color="primary"
          size="large"
          variant="contained"
          onClick={handleGoNext}
          disabled={!canGoNext}
        >
          Find differentially expressed genes
        </NextButton>
      </ButtonsWrapper>
    </div>
  );
}
