import { Icon } from "@blueprintjs/core";
import styled from "@emotion/styled";
import {
  BLUE,
  GRAY,
  LIGHT_GRAY,
  PT_GRID_SIZE_PX,
} from "src/components/common/theme";

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
  padding: ${PT_GRID_SIZE_PX}px;
  border: ${LIGHT_GRAY.C} solid 1px;
  border-radius: 3px;
`;

export const ButtonWrapper = styled.span`
  .bp5-button.bp5-small {
    padding: 0;
  }
`;

export const StyledButton = styled.button`
  color: ${BLUE.B};
  padding: 0;
  border: none;
  background: none;

  &:hover {
    text-decoration: underline;
    cursor: pointer;
  }
`;
