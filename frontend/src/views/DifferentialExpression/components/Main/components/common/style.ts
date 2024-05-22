import styled from "@emotion/styled";
import { Dropdown, fontBodyS } from "@czi-sds/components";
import { fontWeightSemibold, spacesXxs } from "src/common/theme";

export const Label = styled("div")`
  ${fontBodyS}

  font-weight: ${fontWeightSemibold};
  margin-bottom: ${spacesXxs}px;
`;

export const Wrapper = styled("div")`
  display: flex;
  flex-direction: column;
  max-width: 216px;
`;

export const StyledDropdown = styled(Dropdown)`
  width: 100%;
` as typeof Dropdown;
