import { convertPercentageToDiameter } from "../../../../../WheresMyGeneV2/components/HeatMap/components/Chart/utils";
import { Content, Label, LowHigh } from "../../common/style";
import { Dot, Dots, Wrapper } from "./style";

const PERCENTAGES = [0, 0.25, 0.5, 0.75, 1];

export default function ExpressedInCells(): JSX.Element {
  return (
    <Wrapper id="expressed-in-cells">
      <Label id="expressed-in-cells-label">Expressed in Cells (%)</Label>

      <Content>
        <Dots id="expressed-in-cells-dots">
          {PERCENTAGES.map((percentage) => (
            <Dot
              size={convertPercentageToDiameter(percentage)}
              key={percentage}
              data-testid="expressed-in-cells-dots-size"
            />
          ))}
        </Dots>
        <LowHigh className="low-high">
          <span>0</span>
          <span>100</span>
        </LowHigh>
      </Content>
    </Wrapper>
  );
}
