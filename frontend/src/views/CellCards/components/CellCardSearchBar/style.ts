import styled from "@emotion/styled";
import { Autocomplete } from "@mui/material";
import { fontBodyXs, fontCapsXxxs, getColors } from "czifui";

export const StyledAutocomplete = styled(Autocomplete)`
  height: 32px;

  .MuiInputBase-root {
    padding-right: 8px !important;
  }
  & .MuiAutocomplete-inputRoot[class*="MuiInput-root"] {
    padding: 0px;
  }
  &
    .MuiAutocomplete-inputRoot[class*="MuiOutlinedInput-root"]
    .MuiAutocomplete-input {
    padding: 0px;
  }
  & .MuiInputLabel-root {
    margin-top: -8px;
  }
  & .MuiInputLabel-root.MuiInputLabel-shrink {
    margin-top: 0px;
  }
`;

export const SectionTitle = styled.div`
  ${fontCapsXxxs}
  font-weight: 600;
  margin-bottom: 0px !important;
  cursor: default;
  padding-left: 12px !important;

  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[500]};
      `;
  }}
`;

export const SectionItem = styled.li`
  ${fontBodyXs}
  font-weight: 400;
  cursor: pointer;
  height: 32px;
  padding-left: 20px !important;
  margin-bottom: 0px !important;
`;
