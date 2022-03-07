import { Content, Header, LowHigh } from "../../common/style";
import { Dot, Dots, Wrapper } from "./style";

const SIZES = [2, 4, 8, 12, 16];

export default function ExpressedInCells(): JSX.Element {
  return (
    <Wrapper>
      <Header>Expressed in Cells (%)</Header>

      <Content>
        <Dots>
          {SIZES.map((size) => (
            <Dot size={size} key={size} />
          ))}
        </Dots>
        <LowHigh>
          <span>0</span>
          <span>100</span>
        </LowHigh>
      </Content>
    </Wrapper>
  );
}
