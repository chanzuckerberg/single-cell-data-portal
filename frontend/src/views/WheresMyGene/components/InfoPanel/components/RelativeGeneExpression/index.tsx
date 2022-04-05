import { Checkbox } from "czifui";
import Image from "next/image";
import { Content, Header, LowHigh } from "../../common/style";
import plasmaImage from "./plasma.png";
import { ContentWrapper, StyledFormControlLabel, Wrapper } from "./style";

const CONTENT_WIDTH_PX = 100;

interface Props {
  handleIsScaledChange: () => void;
  isScaled: boolean;
}
export default function RelativeGeneExpression({
  handleIsScaledChange,
  isScaled,
}: Props): JSX.Element {
  return (
    <Wrapper>
      <Header>Gene Expression</Header>
      <ContentWrapper>
        <Content>
          <Image
            src={plasmaImage}
            alt="visualization color scale: interpolatePlasma(1.0 - meanExpression)"
            width={CONTENT_WIDTH_PX}
          />
          <LowHigh>
            <span>0.0</span>
            <span>{isScaled ? "1.0" : "6.0"}</span>
          </LowHigh>
        </Content>
        <StyledFormControlLabel
          control={
            <Checkbox
              inputProps={{
                "aria-label":
                  "is scaled? " + isScaled ? "checked" : "unchecked",
                checked: isScaled,
                hidden: true,
              }}
              onChange={handleIsScaledChange}
              stage={isScaled ? "checked" : "unchecked"}
            />
          }
          label="Scaled"
        />
      </ContentWrapper>
    </Wrapper>
  );
}
