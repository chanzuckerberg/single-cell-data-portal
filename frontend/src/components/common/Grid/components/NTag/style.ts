import styled from "@emotion/styled";
import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import {
  CommonThemeProps,
  fontBodyXs,
  getColors,
  Tag as SDSTag,
} from "@czi-sds/components";

const gray100 = (props: CommonThemeProps) => getColors(props)?.gray[100];

export const FieldValues = styled.div`
  color: ${GRAY.A};
  white-space: nowrap;
`;

export const ContentWrapper = styled.div`
  display: flex;
  flex-direction: row;
  padding: ${2 * PT_GRID_SIZE_PX}px;
`;

export const ContentColumn = styled.div`
  display: flex;
  flex-direction: column;

  &:not(:last-child) {
    margin-right: ${PT_GRID_SIZE_PX * 3}px;
  }

  min-width: 160px;
`;

export const Tag = styled(SDSTag)`
  &.MuiChip-root {
    background-color: ${gray100};
    margin: 0;
    padding: 2px 8px;

    .MuiChip-label {
      ${fontBodyXs}
      color: #000000;
      letter-spacing: -0.003em;
    }
  }
`;
