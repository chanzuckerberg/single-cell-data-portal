import { PT_TEXT_COLOR } from "src/components/common/theme";
import styled, { css } from "styled-components";

export const DESCRIPTION_LINE_HEIGHT_PX = 18;
export const MAX_LINE_COUNT = 6;

interface Props {
  isEllipsis: boolean;
}

export const CollectionDescription = styled.div`
  grid-area: description;
`;

export const DescriptionText = styled.p<Props>`
  color: ${PT_TEXT_COLOR};
  letter-spacing: -0.1px;
  line-height: ${DESCRIPTION_LINE_HEIGHT_PX}px;
  margin-bottom: 8px;
  ${(props) => {
    return (
      props.isEllipsis &&
      css`
        -webkit-box-orient: vertical;
        display: -webkit-box;
        -webkit-line-clamp: ${MAX_LINE_COUNT};
        overflow: hidden;
      `
    );
  }};
`;
