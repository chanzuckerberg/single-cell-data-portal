import { Classes, Colors, Intent } from "@blueprintjs/core";
import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const UploadStatusContainer = styled.div`
  color: ${GRAY.A} !important;
`;

export const DatasetTitleCell = styled.td`
  vertical-align: middle !important;
`;

export const TitleContainer = styled.div`
  display: flex;
  flex-direction: row;
  & > .${Classes.CHECKBOX} {
    margin: 0 ${PT_GRID_SIZE_PX}px 0 0;
  }
`;

const intentColorSwitch = (intent: Intent, border?: boolean) => {
  switch (intent) {
    case Intent.DANGER:
      return Colors.RED3;
    case Intent.NONE:
      if (border) return Colors.LIGHT_GRAY3;
      return Colors.GRAY1;
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
