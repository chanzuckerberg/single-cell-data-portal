import styled from "@emotion/styled";
import { CSSProperties } from "react";
import { layout } from "src/components/common/layout";
import { contentWrapper } from "src/components/Layout/style";
import { Container } from "./WheresMyGene/components/HeatMap/style";
import { Wrapper } from "./WheresMyGene/style";

export const ViewGrid = styled.div`
  ${layout}
  ${contentWrapper}
  display: grid;
  grid-template-columns: repeat(8, 1fr);
  margin-top: 10vh;
`;

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
