import styled from "@emotion/styled";
import { Tag, fontBodyXs } from "@czi-sds/components";
import { chipClasses } from "@mui/material/Chip";

export const StyledTag = styled(Tag)`
  font-weight: 400;
  padding: 2px 8px;
  cursor: default !important;
  background-color: rgba(0, 0, 0, 0.05);
  & .${chipClasses.label} {
    color: #000000;
    ${fontBodyXs}
  }
`;
