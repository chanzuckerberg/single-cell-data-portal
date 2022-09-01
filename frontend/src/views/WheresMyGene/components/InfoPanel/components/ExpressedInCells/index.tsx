import { convertPercentageToDiameter } from "../../../HeatMap/utils";
import { Content, Label, LowHigh } from "../../common/style";
import { Dot, Dots, Wrapper } from "./style";

const PERCENTAGES = [0, 0.25, 0.5, 0.75, 1];

export default function ExpressedInCells(): JSX.Element {
  return (
    <Wrapper>
      <Label>Expressed in Cells (%)</Label>

      <Content>
        <Dots>
          {PERCENTAGES.map((percentage) => (
            <Dot
              size={convertPercentageToDiameter(percentage)}
              key={percentage}
            />
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
