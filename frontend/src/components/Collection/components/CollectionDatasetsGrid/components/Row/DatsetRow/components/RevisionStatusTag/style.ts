import { Tag } from "@blueprintjs/core";
import { PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const StyledTag = styled(Tag)`
  width: max-content;
  min-width: fit-content;
  height: ${3 * PT_GRID_SIZE_PX}px;
  margin-left: ${PT_GRID_SIZE_PX}px;
`;
