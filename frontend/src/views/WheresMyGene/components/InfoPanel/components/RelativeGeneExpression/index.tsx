import Image from "next/image";
import { Content, Label, LowHigh } from "../../common/style";
import plasmaImage from "./plasma.png";
import { ContentWrapper, Wrapper } from "./style";

const CONTENT_WIDTH_PX = 120;

export default function RelativeGeneExpression(): JSX.Element {
  return (
    <Wrapper>
      <Label>Gene Expression</Label>
      <ContentWrapper>
        <Content>
          <Image
            src={plasmaImage}
            alt="visualization color scale: interpolateMagma(1.0 - meanExpression)"
            width={CONTENT_WIDTH_PX}
          />
          <LowHigh>
            <span>0</span>
            <span>1</span>
          </LowHigh>
        </Content>
      </ContentWrapper>
    </Wrapper>
  );
}
