import { Wrapper, StepWrapper, StepTitle, StepText } from "./style";
interface Props {
  step: number;
  setStep: (step: number) => void;
}
export default function Steps({ step, setStep }: Props): JSX.Element {
  return (
    <Wrapper>
      <StepWrapper active={step === 1} onClick={() => setStep(1)}>
        <StepTitle>Step 1</StepTitle>
        <StepText>
          Where do you want to find differentially expressed genes?
        </StepText>
      </StepWrapper>
      <StepWrapper active={step === 2} onClick={() => setStep(2)}>
        <StepTitle>Step 2</StepTitle>
        <StepText>
          What do you want to find differentially expressed genes for?
        </StepText>
      </StepWrapper>
      <StepWrapper active={step === 3} onClick={() => setStep(3)}>
        <StepTitle>Step 3</StepTitle>
        <StepText>View differentially expressed genes</StepText>
      </StepWrapper>
    </Wrapper>
  );
}
