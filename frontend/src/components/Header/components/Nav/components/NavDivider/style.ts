import { Divider } from "@mui/material";
import { CommonThemeProps, getColors } from "@czi-sds/components";
import styled from "@emotion/styled";

const gray500 = (props: CommonThemeProps) => getColors(props)?.gray[500];

export const NavDivider = styled(Divider)`
  background-color: ${gray500};
  border: none;
  width: 1px;
`;
