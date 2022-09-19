import { ReactNode } from "react";
import {
  Content,
  Details,
  Header,
  Number,
  NumberContent,
  Wrapper,
} from "./style";

interface Props {
  step: number;
  details: ReactNode;
}

export default function Step({ step, header = "", details }: Props): JSX.Element {
  return (
    <Wrapper>
      <Header>STEP {step}</Header>
      <Details>{details}</Details>
    </Wrapper>
  );
}
