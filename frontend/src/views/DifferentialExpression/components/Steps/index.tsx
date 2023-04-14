import { Wrapper, StepWrapper, StepTitle, StepText } from "./style";
interface Props {
  step: number;
}
export default function Steps({ step }: Props): JSX.Element {
  return (
    <Wrapper>
      <StepWrapper active={step === 1}>
        <StepTitle>Step 1</StepTitle>
        <StepText>
          Where do you want to find differentially expressed genes?
        </StepText>
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
