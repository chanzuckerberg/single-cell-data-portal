import styled from "@emotion/styled";
import { CSSProperties } from "react";
import { contentWrapper } from "src/components/Layout/style";
import { Container } from "./WheresMyGeneV2/components/HeatMap/style_old";
import { Wrapper } from "./WheresMyGene/style";

interface ViewProps {
  overflow?: CSSProperties["overflow"];
}

export const View = styled.div`
  ${contentWrapper}
  grid-area: content;
  overflow: ${(props: ViewProps) =>
    props.overflow ? props.overflow : undefined};

  &.CLONED {
    overflow: hidden;
    height: fit-content;
    width: fit-content;

    & ${Wrapper}, ${Container} {
      height: fit-content;
      overflow: hidden;
    }
  }
`;
