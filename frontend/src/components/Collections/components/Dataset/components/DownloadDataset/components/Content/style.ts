import { AnchorButton, Button } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { PT_GRID_SIZE_PX } from "src/components/common/theme";

export const Wrapper = styled.div`
  width: 625px;
  height: 300px;
`;

export const DownloadButton = styled(AnchorButton)`
  margin-right: ${PT_GRID_SIZE_PX}px;
`;

export const CancelButton = styled(Button)`
  margin-right: ${PT_GRID_SIZE_PX}px;
`;
