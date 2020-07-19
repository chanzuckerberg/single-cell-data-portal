import React, { FC } from "react";
import questionMarkSvg from "src/components/ProjectList/components/Heading/questionMark.svg";
import {
  Info,
  Name,
  QuestionMark,
  QuestionMarkWrapper,
  View,
  Wrapper,
} from "./style";

const TOOLTIP_MESSAGE =
  "A minimally harmonized dataset as well as the originally submitted dataset are available for exploration within cellxgene. The harmonized dataset may have metadata labels that have been changed from the original dataset in order to create consistency across the collection of datasets as a whole.";

const Heading: FC = () => {
  return (
    <Wrapper>
      <Name>Name of dataset</Name>
      <View>
        View dataset in cellxgene
        <QuestionMarkWrapper title={TOOLTIP_MESSAGE}>
          <QuestionMark src={questionMarkSvg} alt="question mark" />
        </QuestionMarkWrapper>
      </View>
      <Info>More information</Info>
    </Wrapper>
  );
};

export default Heading;
