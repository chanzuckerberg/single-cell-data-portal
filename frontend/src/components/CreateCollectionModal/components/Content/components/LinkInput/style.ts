import { Icon } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { GRAY, RED } from "src/components/common/theme";

interface Props {
  warning?: boolean;
}

export const CollectionLink = styled.div`
  column-gap: 16px;
  display: grid;
  grid-template-columns: 160px 1fr 1fr; /* grid columns for collection link selector, name and url */
  position: relative; /* positions close collection link icon */

  /* DOI collection link url */
  [for="DOI"] {
    grid-column: 2 / -1;
  }
`;

export const CloseCollectionLinkIcon = styled(Icon)`
  cursor: pointer;
  padding: 1px;
  position: absolute;
  right: 0;
  top: 0;
`;

export const InputPrefix = styled.span<Props>`
  color: ${(props) => (props.warning ? RED.C : GRAY.A)};
  display: block;
  letter-spacing: -0.1px;
  line-height: 18px;
  padding: 7px 0 7px 8px;
`;

export const HelperText = styled.div<Props>`
  color: ${(props) => (props.warning ? RED.C : GRAY.A)};
  font-size: 12px;
  grid-column: 1 / -1;
  line-height: 15px;
  margin-top: 8px;
`;
