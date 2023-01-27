import Image from "next/image";
import { Content, Label, LowHigh } from "../../common/style";
import plasmaImage from "./plasma.png";
import { ContentWrapper, Wrapper } from "./style";

const CONTENT_WIDTH_PX = 120;

interface Props {
  isScaled: boolean;
}
export default function RelativeGeneExpression({
  isScaled,
}: Props): JSX.Element {
  return (
    <Wrapper id="relative-gene-expression">
      <Label id="relative-gene-expression-label">Gene Expression</Label>
      <ContentWrapper>
        <Content>
          <Image
            id="visualization-color-scale"
            src={plasmaImage}
            alt="visualization color scale: interpolateMagma(1.0 - meanExpression)"
            width={CONTENT_WIDTH_PX}
          />
          <LowHigh className="low-high">
            <span>0.0</span>
            <span>{isScaled ? "1.0" : "6.0"}</span>
          </LowHigh>
        </Content>
      </ContentWrapper>
    </Wrapper>
  );
}
