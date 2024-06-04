import styled from "@emotion/styled";
import {
  ButtonDropdown as SDSButtonDropdown,
  fontBodyS,
} from "@czi-sds/components";
import { grey500, spacesS, spacesXxxs } from "src/common/theme";

export const ButtonDropdown = styled(SDSButtonDropdown)`
  ${fontBodyS}
  color: ${grey500};
  font-weight: 500;
  margin: ${spacesXxxs}px 0 ${spacesS}px;
  padding: 0;
  text-transform: capitalize;
`;
