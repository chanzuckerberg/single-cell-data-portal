import React from "react";
import { StepTitle } from "../../../Steps/style";
import {
  StepOneHeader,
  StepOneSubHeader,
  WordPop,
  FiltersWrapper,
  NextButton,
} from "./style";
import Filters from "./components/Filters";

interface Props {
  setStep: (step: number) => void;
}
export default function StepOne({ setStep }: Props): JSX.Element {
  const handleGoNext = () => {
    setStep(2);
  };

  return (
    <div>
      <StepTitle>Step 1</StepTitle>
      <StepOneHeader>
        <WordPop>Where</WordPop> do you want to find differentially expressed
        genes?
      </StepOneHeader>
      <StepOneSubHeader>
        Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod
        tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim
        veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea
        commodo consequat. Duis aute irure dolor in reprehenderit in voluptate
        velit esse cillum dolore eu fugiat nulla pariatur.
      </StepOneSubHeader>
      <FiltersWrapper>
        <Filters />
      </FiltersWrapper>
      <NextButton
        color="primary"
        size="large"
        variant="contained"
        onClick={handleGoNext}
      >
        Next: Build query group
      </NextButton>
    </div>
  );
}
