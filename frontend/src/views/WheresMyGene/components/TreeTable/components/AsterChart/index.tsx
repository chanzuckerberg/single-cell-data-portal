import { Container, FirstPart, SecondPart } from "./style";

interface Props {
  relativeExpression: number;
  proportionalExpression: number;
}

export default function AsterChart({
  relativeExpression,
  proportionalExpression,
}: Props): JSX.Element {
  return (
    <Container>
      <FirstPart
        relativeExpression={relativeExpression}
        proportionalExpression={proportionalExpression}
      >
        <SecondPart
          relativeExpression={relativeExpression}
          proportionalExpression={proportionalExpression}
        />
      </FirstPart>
    </Container>
  );
}
