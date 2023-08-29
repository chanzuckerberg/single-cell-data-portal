import styled from "@emotion/styled";
import { fontBodyXxs } from "@czi-sds/components";
import { gray500 } from "src/common/theme";
import { autocompleteClasses, menuItemClasses, Popper } from "@mui/material";
import { Button, ButtonIcon, MenuItem } from "@czi-sds/components";
import { OFF_WHITE } from "src/common/theme";

export const GENE_SEARCH_BAR_HEIGHT_PX = 32;

export const Container = styled.div`
  height: ${GENE_SEARCH_BAR_HEIGHT_PX}px;
  width: fit-content;
  margin-bottom: 8px;
`;

export const AutocompleteWrapper = styled.div`
  width: 240px;
`;

export const ActionWrapper = styled.div`
  display: flex;
  gap: 16px;
`;

export const Label = styled.label`
  ${fontBodyXxs}

  color: ${gray500};
`;

export const LoadingIndicatorWrapper = styled.div`
  display: flex;
  align-items: center;
`;

export const StyledButtonWrapper = styled.div`
  align-self: center;
`;

export const StyledClearButton = styled(Button)`
  white-space: nowrap;
  font-weight: 500;
`;
