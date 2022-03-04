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
  header: string;
  details: ReactNode;
}

export default function Step({ step, header, details }: Props): JSX.Element {
  return (
    <Wrapper>
      <Number>
        <NumberContent>{step}</NumberContent>
      </Number>
      <Content>
        <Header>{header}</Header>
        <Details>{details}</Details>
      </Content>
    </Wrapper>
  );
}
