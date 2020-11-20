import { Button, Icon } from "@blueprintjs/core";
import { GRAY, LIGHT_GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
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
  padding: ${PT_GRID_SIZE_PX}px;
  border: ${LIGHT_GRAY.C} solid 1px;
  border-radius: 3px;
`;

export const ButtonWrapper = styled.span`
  .bp3-button.bp3-small {
    padding: 0;
  }
`;

export const StyledButton = styled(Button)`
  && {
    &:hover {
      background: unset !important;
      text-decoration: underline;
    }
  }
`;
