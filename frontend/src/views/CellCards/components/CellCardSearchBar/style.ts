import styled from "@emotion/styled";
import { Autocomplete } from "@mui/material";
import { fontBodyXs, fontCapsXxxs, getColors } from "czifui";

export const StyledAutocomplete = styled(Autocomplete)`
  border-radius: 4px;
  .MuiInputBase-root {
    padding-right: 8px !important;
  }
`;

export const SectionTitle = styled.div`
  ${fontCapsXxxs}
  font-weight: 600;
  margin-bottom: 10px;

  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[500]};
      `;
  }}
`;

export const SectionItem = styled.div`
  ${fontBodyXs}
  font-weight: 400;
  margin-bottom: 12px;
  padding-left: 8px;
  cursor: pointer;
`;
