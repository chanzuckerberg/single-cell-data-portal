import { fontBodyXs } from "@czi-sds/components";
import styled from "@emotion/styled";
import { Autocomplete, TextField } from "@mui/material";
import { gray400, primary400 } from "src/common/theme";
import { FilterOption } from "../../types";
import { formControlClasses } from "@mui/material/FormControl";
import { inputBaseClasses } from "@mui/material/InputBase";
import { formLabelClasses } from "@mui/material/FormLabel";
import { autocompleteClasses } from "@mui/material/Autocomplete";

console.log();

const Tag = styled.div`
  display: inline-flex;
  align-items: center;
  padding: 4px 8px;
  color: white;
  border-radius: 4px;
  ${fontBodyXs}
  font-weight: 500;
  margin-right: 5px;
  cursor: default;
`;

export const PrimaryTag = styled(Tag)`
  background-color: ${primary400};
`;

export const GrayTag = styled(Tag)`
  background-color: ${gray400};
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

  & .${formControlClasses.root} {
    overflow-x: auto;
    max-width: 300px;
    padding-top: 8px;
    background-color: transparent;
  }
  & .${formLabelClasses.root} {
    margin-top: 8px;
  }
  & .${inputBaseClasses.root} {
    display: flex;
    flex-direction: row;
    align-items: center;
    flex-wrap: nowrap;
    width: fit-content;
    justify-content: space-between;
    min-width: 300px;
    padding: 8px !important;
    background-color: white;
  }
  & .${autocompleteClasses.tag} {
    white-space: nowrap;
    width: fit-content;
  }
  & + .${autocompleteClasses.popper} {
    max-width: 300px;
  }
  & .${autocompleteClasses.endAdornment} {
    position: relative;
    display: flex;
    flex-direction: row;
    right: unset !important;
    order: -1;
  }
  & .${inputBaseClasses.input} {
    width: 45px !important;
    order: -2;
  }
`;
