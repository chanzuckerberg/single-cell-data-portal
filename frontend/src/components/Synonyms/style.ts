import styled from "@emotion/styled";
import { fontBodyXs } from "@czi-sds/components";
import { fontWeightRegular, textSecondary } from "src/common/theme";

export const Wrapper = styled.div`
  display: flex;
  padding-bottom: 8px;
  align-items: first baseline;
`;

export const Label = styled.div`
  ${fontBodyXs}

  color: ${textSecondary};
  font-weight: ${fontWeightRegular};
  width: 80px;
`;

export const Synonym = styled.div`
  ${fontBodyXs}

  /**
   * Expands to fill all available horizontal space
   * flex-grow, flex-shrink, and flex-basis
   * Default is: 0 1 auto
   * https://css-tricks.com/snippets/css/a-guide-to-flexbox/#aa-flex
   */
  flex: 1 0 0;
  padding: 4px 8px;
  color: black;
  font-weight: ${fontWeightRegular};
`;
