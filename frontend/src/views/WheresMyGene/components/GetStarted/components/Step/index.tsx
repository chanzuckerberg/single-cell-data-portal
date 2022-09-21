import { ReactNode } from "react";
import { Details, Header, Wrapper } from "./style";

interface Props {
  step: number;
  details: ReactNode;
}

export default function Step({ step, details }: Props): JSX.Element {
  return (
    <Wrapper>
      <Header>STEP {step}</Header>
      <Details>{details}</Details>
    </Wrapper>
  );
}
