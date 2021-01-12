import { Classes, Intent } from "@blueprintjs/core";
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
  margin-top: ${2 * PT_GRID_SIZE_PX}px;
  display: flex;
  flex-direction: row;
  vertical-align: middle;
  & > .${Classes.SPINNER}, .${Classes.ICON} {
    margin: auto ${PT_GRID_SIZE_PX}px auto 0;
    height: 100%;
  }
`;
