import styled from "@emotion/styled";
import { Tag, fontBodyXs } from "@czi-sds/components";
import { chipClasses } from "@mui/material/Chip";

interface StyledTagProps {
  isSingleCellType: boolean;
}
export const StyledTag = styled(Tag)<StyledTagProps>`
  font-weight: 400;
  padding: 2px 8px;
  cursor: ${({ isSingleCellType }) =>
    isSingleCellType ? "pointer" : "default"} !important;
  background-color: rgba(0, 0, 0, 0.05);
  & .${chipClasses.label} {
    color: #000000;
    ${fontBodyXs}
  }
`;
