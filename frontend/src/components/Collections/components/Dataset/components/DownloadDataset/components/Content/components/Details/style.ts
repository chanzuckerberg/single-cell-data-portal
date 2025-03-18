import { fontBodyXxs, fontHeaderS, Icon } from "@czi-sds/components";
import styled from "@emotion/styled";
import { FormControl as MFormControl } from "@mui/material";
import { cornersM, cornersL, spacesM, spacesS, textPrimary, gray100, gray500 } from "src/common/theme";

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

export const NoneSelected = styled.div`
  align-items: center;
  display: flex;
  flex-direction: column;
  gap: ${spacesS}px;
  justify-content: center;
  min-height: 150px;
  border-radius: 8px;
  background-color: ${gray100};
  margin-top: ${spacesS}px;
  h4{
    ${fontHeaderS}  
    font-weight: 600;
    margin-bottom: 0;
  }
  p{
    ${fontBodyXxs}  
    color: ${gray500};
    margin-bottom: 0;
  }
`;