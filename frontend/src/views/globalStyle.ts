import { layout } from "src/components/common/layout";
import { contentWrapper } from "src/components/Layout/style";
import styled from "styled-components";

interface Props {
  isFilterEnabled?: boolean;
}

export const ViewGrid = styled.div<Props>`
  ${layout}
  ${contentWrapper}
  display: grid;
  grid-template-columns: repeat(8, 1fr);
  margin-top: 10vh;
  min-width: ${(props) =>
    props.isFilterEnabled
      ? "unset"
      : undefined}; /* overrides layout min-width specification to facilitate shrink to fit responsiveness */
`;

export const View = styled.div`
  ${contentWrapper}
  grid-area: content;
  overflow: auto; /* facilitates independent content scrolling for sidebar layout */
`;
