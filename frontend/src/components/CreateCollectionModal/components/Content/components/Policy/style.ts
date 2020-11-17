import { Icon } from "@blueprintjs/core";
import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const StyledIcon = styled(Icon)`
  margin: 6px;
`;

export const Wrapper = styled.div`
  margin-top: ${PT_GRID_SIZE_PX}px;
`;

export const BulletWrapper = styled.div`
  display: flex;
  margin-bottom: ${PT_GRID_SIZE_PX}px;
  line-height: ${2 * PT_GRID_SIZE_PX}px;
`;

export const Text = styled.div`
  color: ${GRAY.A};
`;

export const ContentWrapper = styled.div`
  overflow-y: scroll;
  height: ${20 * PT_GRID_SIZE_PX}px;
`;
