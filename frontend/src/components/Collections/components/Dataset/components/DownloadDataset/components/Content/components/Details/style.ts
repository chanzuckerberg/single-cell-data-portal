import { Icon } from "@czi-sds/components";
import styled from "@emotion/styled";
import { FormControl as MFormControl } from "@mui/material";
import { cornersM, spacesM, spacesS, textPrimary } from "src/common/theme";

export const FormControl = styled(MFormControl)`
  min-height: 240px;
`;

export const SeuratNotice = styled.div`
  align-items: center;
  background-color: #ffefcf;
  border-radius: ${cornersM}px;
  display: flex;
  font-size: 13px;
  color: ${textPrimary};
  gap: ${spacesS}px;
  margin: 0 0 ${spacesM}px;
  padding: ${spacesS}px ${spacesM}px;
`;

export const StyledIcon = styled(Icon)`
  color: #b77800;
`;

export const StyledLink = styled.a`
  color: #b77800;
  text-decoration: underline;
  margin: 0;
`;

export const TextWrapper = styled.span`
  display: inline;
`;
