import { useState } from "react";
import Steps from "../Steps";
import DeWizard from "../DeWizard";
import { Wrapper } from "./style";

export default function DifferentialExpression(): JSX.Element {
  const [step, setStep] = useState(1);

  return (
    <Wrapper>
      <Steps step={step} />
      <DeWizard step={step} setStep={setStep} />
    </Wrapper>
  );
}
