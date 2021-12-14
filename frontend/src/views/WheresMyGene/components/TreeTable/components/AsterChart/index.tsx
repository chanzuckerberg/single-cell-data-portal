import { Container, FirstPart, SecondPart } from "./style";

interface Props {
  colorValue: number;
  degreeValue: number;
}

export default function AsterChart({
  colorValue,
  degreeValue,
}: Props): JSX.Element {
  return (
    <Container>
      <FirstPart colorValue={colorValue} degreeValue={degreeValue}>
        <SecondPart colorValue={colorValue} degreeValue={degreeValue} />
      </FirstPart>
    </Container>
  );
}
