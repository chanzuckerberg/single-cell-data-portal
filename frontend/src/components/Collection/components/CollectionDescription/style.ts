import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyS } from "@czi-sds/components";
import { spacesXxs } from "src/common/theme";

export const DESCRIPTION_LINE_HEIGHT_PX = 20;
export const MAX_LINE_COUNT = 6;

interface Props extends CommonThemeProps {
  isEllipsis: boolean;
}

export const CollectionDescription = styled.div`
  display: grid;
  gap: ${spacesXxs}px;
  grid-area: description;
  justify-items: flex-start;
`;

export const DescriptionText = styled.p<Props>`
  ${fontBodyS}
  letter-spacing: -0.006em;
  margin-bottom: 0;
  ${(props) => {
    return (
      props.isEllipsis &&
      `
        -webkit-box-orient: vertical;
        display: -webkit-box;
        -webkit-line-clamp: ${MAX_LINE_COUNT};
        overflow: hidden;
      `
    );
  }}
`;
