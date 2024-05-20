import React from "react";
import ReactMarkdown from "react-markdown";
import { useConnect } from "./connect";
import { CardWrapper, CardContent, CloseButtonWrapper } from "./style";
import { Props } from "./types";
import { ButtonIcon } from "@czi-sds/components";

const InterpretationCard = ({
  isQueryGroup1,
  differentialExpressionResults,
  setIsLoadingInterpret,
  setIsVisible,
}: Props) => {
  const { analysisMessage, isLoading } = useConnect({
    isQueryGroup1,
    differentialExpressionResults,
    setIsLoadingInterpret,
  });

  return !isLoading ? (
    <CardWrapper>
      <CloseButtonWrapper>
        <ButtonIcon
          sdsType="tertiary"
          sdsIcon="xMark"
          sdsSize="small"
          onClick={() => setIsVisible(false)}
        />
      </CloseButtonWrapper>
      <CardContent>
        <ReactMarkdown>{analysisMessage}</ReactMarkdown>
      </CardContent>
    </CardWrapper>
  ) : null;
};

export default InterpretationCard;
