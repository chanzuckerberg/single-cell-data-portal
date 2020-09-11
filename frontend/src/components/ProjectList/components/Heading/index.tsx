import React, { FC } from "react";
import questionMarkSvg from "src/components/ProjectList/components/Heading/questionMark.svg";
import { Download, Info, Name, QuestionMark, View, Wrapper } from "./style";

const TOOLTIP_MESSAGE = `cellxgene augments datasets with a minimal set of metadata fields designed to enable comparisons across datasets. In cases where these columns conflict with author's metadata, the author's columns are prefixed by "original_"`;

const Heading: FC = () => {
  return (
    <Wrapper>
      <Name>Dataset name</Name>
      <View>
        View in cellxgene
        <span title={TOOLTIP_MESSAGE}>
          <QuestionMark src={questionMarkSvg} alt="question mark" />
        </span>
      </View>
      <Download>Download dataset</Download>
      <Info>More information</Info>
    </Wrapper>
  );
};

export default Heading;
