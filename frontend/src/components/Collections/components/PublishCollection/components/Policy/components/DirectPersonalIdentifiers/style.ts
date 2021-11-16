import { Icon } from "@blueprintjs/core";
import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const ContentWrapper = styled.div`
  font-size: 12px;
  width: ${70 * PT_GRID_SIZE_PX}px;
  padding: ${2 * PT_GRID_SIZE_PX}px;
`;

export const FirstSentence = styled.div`
  margin-bottom: ${PT_GRID_SIZE_PX}px;
  line-height: ${2 * PT_GRID_SIZE_PX}px;
`;

export const StyledIcon = styled(Icon)`
  margin: 6px;
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
