import React, { FC } from "react";
import questionMarkSvg from "src/components/Collections/components/Heading/questionMark.svg";
import { SmallColumn } from "../Dataset/common/style";
import { Name, QuestionMark, Wrapper } from "./style";

const TOOLTIP_MESSAGE = `cellxgene augments datasets with a minimal set of metadata fields designed to enable comparisons across datasets. In cases where these columns conflict with author's metadata, the author's columns are prefixed by "original_"`;

const Heading: FC = () => {
  return (
    <Wrapper>
      <Name>Dataset name</Name>
      <SmallColumn>
        View in cellxgene
        <span title={TOOLTIP_MESSAGE}>
          <QuestionMark src={questionMarkSvg} alt="question mark" />
        </span>
      </SmallColumn>
      <SmallColumn>Download dataset</SmallColumn>
      <SmallColumn>More information</SmallColumn>
    </Wrapper>
  );
};

export default Heading;
