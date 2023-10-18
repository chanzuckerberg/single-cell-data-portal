import { useRef, useEffect } from "react";
import { StyledText } from "./style";

interface TextProps {
  name: string;
  maxWidth: number;
  isInCorpus: boolean;
  cursor?: string;
}
export default function Text({
  isInCorpus,
  name,
  maxWidth,
  cursor = "pointer",
}: TextProps) {
  const textRef = useRef<SVGTextElement>(null);

  useEffect(() => {
    const textElement = textRef.current;
    let textContent = name;
    if (textElement) {
      while (textElement.getComputedTextLength() > maxWidth) {
        textContent = textContent.slice(0, -1);
        textElement.textContent = textContent + "...";
      }
    }
  }, [name, maxWidth]);

  return (
    <StyledText
      ref={textRef}
      dy=".33em"
      fontSize={9}
      fontWeight={400}
      fontFamily="Arial"
      textAnchor="left"
      isInCorpus={isInCorpus}
      cursor={cursor}
      dx={10}
    >
      {name}
    </StyledText>
  );
}
