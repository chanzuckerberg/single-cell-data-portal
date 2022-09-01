import { Checkbox, Icon, Tooltip } from "czifui";
import Image from "next/image";
import { Content, Label, LowHigh } from "../../common/style";
import plasmaImage from "./plasma.png";
import {
  ContentWrapper,
  FlexDiv,
  LabelWrapper,
  StyledFormControlLabel,
  Wrapper,
} from "./style";

const CONTENT_WIDTH_PX = 120;

interface Props {
  handleIsScaledChange: () => void;
  isScaled: boolean;
  showScaled?: boolean;
}
export default function RelativeGeneExpression({
  handleIsScaledChange,
  isScaled,
  showScaled,
}: Props): JSX.Element {
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
            <span>0.0</span>
            <span>{isScaled ? "1.0" : "6.0"}</span>
          </LowHigh>
        </Content>
        {showScaled ? (
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
            label={
              <LabelWrapper>
                <span>Scaled</span>
                <Tooltip title="Expression is scaled to the range [0,1]. Scaling is done by assigning the minimum value in the current view to 0 and the max is assigned to 1.">
                  <FlexDiv>
                    <Icon sdsIcon="infoCircle" sdsSize="s" sdsType="static" />
                  </FlexDiv>
                </Tooltip>
              </LabelWrapper>
            }
          />
        ) : (
          ""
        )}
      </ContentWrapper>
    </Wrapper>
  );
}
