import Image from "next/image";
import { Content, Label, LowHigh } from "../../common/style";
import plasmaImage from "./plasma.png";
import { ContentWrapper, Wrapper } from "./style";
import { MAX_EXPRESSION_LABEL_TEST_ID } from "./constants";

const CONTENT_WIDTH_PX = 120;

interface Props {
  isScaled: boolean;
  maxExpression: number;
}
export default function RelativeGeneExpression({
  isScaled,
  maxExpression,
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
          {(isScaled || maxExpression !== -Infinity) && (
            <LowHigh className="low-high">
              <span>0.0</span>
              <span data-testid={MAX_EXPRESSION_LABEL_TEST_ID}>
                {isScaled ? "1.0" : maxExpression.toFixed(2)}
              </span>
            </LowHigh>
          )}
        </Content>
      </ContentWrapper>
    </Wrapper>
  );
}
