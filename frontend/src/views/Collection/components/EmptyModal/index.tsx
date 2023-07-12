import { H4 } from "@blueprintjs/core";
import { FC, ReactChild } from "react";
import { Border, CenterAlignedDiv } from "./style";

interface Props {
  title: string;
  content?: ReactChild;
  button?: ReactChild;
}

const EmptyModal: FC<Props> = ({ title, content, button }) => {
  return (
    <Border>
      <CenterAlignedDiv>
        <H4>{title}</H4>
        {content}
        {button}
      </CenterAlignedDiv>
    </Border>
  );
};

export default EmptyModal;
