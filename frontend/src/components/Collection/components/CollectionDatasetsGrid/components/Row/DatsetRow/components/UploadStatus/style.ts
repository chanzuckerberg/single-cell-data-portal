import { Intent, Spinner } from "@blueprintjs/core";
import {
  GRAY,
  LIGHT_GRAY,
  PT_GRID_SIZE_PX,
  RED,
} from "src/components/common/theme";
import styled from "styled-components";

const intentColorSwitch = (intent: Intent, border?: boolean) => {
  switch (intent) {
    case Intent.DANGER:
      return RED.C;
    case Intent.NONE:
      if (border) return LIGHT_GRAY.C;
      return GRAY.A;
  }
};

export const StatusContainer = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
`;

interface Props {
  intent?: Intent;
}

export const DatasetStatusTag = styled.div`
  color: ${(props: Props) => intentColorSwitch(props.intent ?? Intent.NONE)};
  border: 1px solid
    ${(props: Props) => intentColorSwitch(props.intent ?? Intent.NONE, true)};
  border-radius: 3px;
  align-self: flex-start;
  width: fit-content;
  padding: ${PT_GRID_SIZE_PX}px;
  display: flex;
  flex-direction: row;
  vertical-align: middle;
  margin-right: ${PT_GRID_SIZE_PX}px;
`;

export const StyledSpinner = styled(Spinner)`
  margin-right: ${PT_GRID_SIZE_PX}px;
`;
