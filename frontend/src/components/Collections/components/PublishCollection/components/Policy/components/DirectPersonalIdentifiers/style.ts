import { Icon } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";

export const ContentWrapper = styled.div`
  font-size: 12px;
  width: ${70 * PT_GRID_SIZE_PX}px;
  padding: ${2 * PT_GRID_SIZE_PX}px;
  /* (thuang): Arbitrary number that works for different vh and zoom levels */
  height: 50vh;
  overflow-y: scroll;
`;

export const FirstSentence = styled.div`
  margin-bottom: ${PT_GRID_SIZE_PX}px;
  line-height: ${2 * PT_GRID_SIZE_PX}px;
`;

export const StyledIcon = styled(Icon)`
  margin-right: 8px;
`;

export const BulletWrapper = styled.div`
  display: flex;
  margin-bottom: ${PT_GRID_SIZE_PX}px;
  line-height: 14px;
  align-items: baseline;
`;

export const Text = styled.div`
  color: ${GRAY.A};
`;
