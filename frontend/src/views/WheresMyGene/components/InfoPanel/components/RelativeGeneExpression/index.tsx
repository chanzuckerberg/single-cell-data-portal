import { Content, Header, LowHigh } from "../../common/style";
import { Dot, Dots, Wrapper } from "./style";

const COLORS = ["#FFFFC1", "#FDD66D", "#F95929", "#CD0019", "#6C001D"];

export default function RelativeGeneExpression(): JSX.Element {
  return (
    <Wrapper>
      <Header>Relative Gene Expression</Header>
      <Content>
        <Dots>
          {COLORS.map((color) => (
            <Dot color={color} key={color} />
          ))}
        </Dots>
        <LowHigh>
          <span>Low</span>
          <span>High</span>
        </LowHigh>
      </Content>
    </Wrapper>
  );
}
