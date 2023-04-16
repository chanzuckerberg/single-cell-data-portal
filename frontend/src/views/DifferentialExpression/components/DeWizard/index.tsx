import { Wrapper } from "./style";
import StepOne from "./components/StepOne";
import StepTwo from "./components/StepTwo";
import StepThree from "./components/StepThree";

interface Props {
  step: number;
  setStep: (step: number) => void;
}

export default function DeWizard({ step, setStep }: Props): JSX.Element {
  const stepComponents = [StepOne, StepTwo, StepThree];

  const renderStep = (Component: React.ElementType) => (
    <Component setStep={setStep} />
  );

  return <Wrapper>{renderStep(stepComponents[step - 1])}</Wrapper>;
}
