import { memo } from "react";
import { Container, FirstPart } from "./style";

interface Props {
  colorValue: number;
  degreeValue: number;
}

export default memo(function AsterChart({
  colorValue,
  degreeValue,
}: Props): JSX.Element {
  return (
    <Container>
      <FirstPart colorValue={colorValue} degreeValue={degreeValue} />
    </Container>
  );
});
