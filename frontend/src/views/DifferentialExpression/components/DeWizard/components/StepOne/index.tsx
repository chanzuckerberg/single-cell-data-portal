import React, { useContext } from "react";
import { StepTitle } from "../../../Steps/style";
import {
  StepHeader,
  StepSubHeader,
  WordPop,
  FiltersWrapper,
  NextButton,
} from "./style";
import Filters from "./components/Filters";
import { StateContext } from "src/views/DifferentialExpression/common/store";

interface Props {
  setStep: (step: number) => void;
}
export default function StepOne({ setStep }: Props): JSX.Element {
  const { selectedFilters } = useContext(StateContext);

  const handleGoNext = () => {
    setStep(2);
  };

  return (
    <div>
      <StepTitle>Step 1</StepTitle>
      <StepHeader>
        <WordPop>Where</WordPop> do you want to find differentially expressed
        genes?
      </StepHeader>
      <StepSubHeader>
        Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod
        tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim
        veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea
        commodo consequat. Duis aute irure dolor in reprehenderit in voluptate
        velit esse cillum dolore eu fugiat nulla pariatur.
        <br />
        <br />
        Filters marked by {"\u002A"} require at least one field.
      </StepSubHeader>
      <FiltersWrapper>
        <Filters />
      </FiltersWrapper>
      <NextButton
        color="primary"
        size="large"
        variant="contained"
        onClick={handleGoNext}
        disabled={!selectedFilters.tissues.length}
      >
        Next: Build query group
      </NextButton>
    </div>
  );
}
