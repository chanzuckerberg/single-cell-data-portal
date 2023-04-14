import { useContext, useMemo } from "react";
import { StateContext } from "../../common/store";
import {
  Wrapper,
  StepWrapper,
  StepTitle,
  StepText,
  StepOneSelectedFilters,
} from "./style";
import { useAvailableOrganisms } from "src/common/queries/differentialExpression";

interface Props {
  step: number;
}
export default function Steps({ step }: Props): JSX.Element {
  const { organismId, selectedFilterNames } = useContext(StateContext);
  const { data: organisms } = useAvailableOrganisms();
  const organismName = useMemo(() => {
    let result = "";

    if (!organisms || !organismId) return result;

    for (const organism of organisms) {
      if (organism.id === organismId) {
        result = organism.name;
        break;
      }
    }
    return result;
  }, [organisms]);
  const selectedTissues = selectedFilterNames.tissues;
  const selectedTissuesText = selectedTissues.length
    ? `, ${selectedTissues.join(", ")}`
    : "";
  const stepOneSelectedFiltersText = `${organismName}${selectedTissuesText}`;
  return (
    <Wrapper>
      <StepWrapper active={step === 1}>
        <StepTitle>Step 1</StepTitle>
        <StepText>
          Where do you want to find differentially expressed genes?
        </StepText>
        <StepOneSelectedFilters>
          {stepOneSelectedFiltersText}
        </StepOneSelectedFilters>
      </StepWrapper>
      <StepWrapper active={step === 2}>
        <StepTitle>Step 2</StepTitle>
        <StepText>
          What do you want to find differentially expressed genes for?
        </StepText>
      </StepWrapper>
      <StepWrapper active={step === 3}>
        <StepTitle>Step 3</StepTitle>
        <StepText>View differentially expressed genes</StepText>
      </StepWrapper>
    </Wrapper>
  );
}
