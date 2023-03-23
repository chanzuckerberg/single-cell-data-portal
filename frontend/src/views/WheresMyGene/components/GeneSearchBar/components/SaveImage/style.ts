import { Button } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { PT_GRID_SIZE_PX } from "src/components/common/theme";

export const DownloadButton = styled(Button)`
  margin-right: ${PT_GRID_SIZE_PX}px;
`;
export const StyledDiv = styled.div`
  display: flex;
  flex-direction: column;
  margin: 0;
`;

export const ButtonWrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 0 15px;
`;
