import { fontBodyXs } from "@czi-sds/components";
import styled from "@emotion/styled";
import { Autocomplete, TextField } from "@mui/material";
import { primary400 } from "src/common/theme";
import { FilterOption } from "../../types";

export const Tag = styled.div`
  display: inline-flex;
  align-items: center;
  padding: 4px 8px;
  background-color: ${primary400};
  color: white;
  border-radius: 4px;
  ${fontBodyXs}
  font-weight: 500;
  margin-right: 5px;
  cursor: default;
`;

export const CloseIcon = styled.span`
  display: inline-block;
  margin-left: 8px;
  cursor: pointer;
  &:after {
    content: "×";
  }
`;

export const StyledTextField = styled(TextField)`
  height: 100%;
  width: 300px;
  background-color: white;
`;

export const StyledAutocomplete = styled(
  Autocomplete<FilterOption, true, false, false>
)`
  width: 300px;

  & .MuiAutocomplete-inputRoot[class*="MuiInput-root"] {
    padding: 4px 0px;
  }
  &
    .MuiAutocomplete-inputRoot[class*="MuiOutlinedInput-root"]
    .MuiAutocomplete-input {
    padding: 4px 0px;
  }
  & .MuiInputLabel-root {
    margin-top: -4px;
    z-index: 0;
  }
  & .MuiInputLabel-root.MuiInputLabel-shrink {
    margin-top: 0px;
  }
`;
