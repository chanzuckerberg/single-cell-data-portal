import { useRef, useEffect } from "react";
import { defaultTextColor } from "../../../../common/constants";

interface TextProps {
  name: string;
  maxWidth: number;
  isInCorpus: boolean;
}
export default function Text({ isInCorpus, name, maxWidth }: TextProps) {
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
  }, [name]);

  return (
    <text
      ref={textRef}
      dy=".33em"
      fontSize={9}
      fontFamily="Arial"
      textAnchor="left"
      fontStyle={
        isInCorpus || name.endsWith("cell types") ? "normal" : "italic"
      }
      fontWeight={isInCorpus ? "bold" : "normal"}
      fill={defaultTextColor}
      dx={10}
    >
      {name}
    </text>
  );
}
