import styled from "@emotion/styled";
import { layout } from "src/components/common/layout";
import { contentWrapper } from "src/components/Layout/style";

export const ViewGrid = styled.div`
  ${layout}
  ${contentWrapper}
  display: grid;
  grid-template-columns: repeat(8, 1fr);
  margin-top: 10vh;
`;

interface ViewProps {
  hideOverflow?: boolean;
}

export const View = styled.div`
  ${contentWrapper}
  grid-area: content;
  overflow: ${(props: ViewProps) =>
    props.hideOverflow
      ? "hidden"
      : "auto"}; /* facilitates independent content scrolling for sidebar layout */
`;
