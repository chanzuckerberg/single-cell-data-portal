import { PT_TEXT_COLOR } from "src/components/common/theme";
import styled from "styled-components";

export const CollectionDetails = styled.div`
  align-items: flex-start;
  display: grid;
  gap: 0 40px;
  grid-template-columns: 1fr 1fr;
  margin: 16px 0 44px;
`;

export const CollectionDescription = styled.div`
  color: ${PT_TEXT_COLOR};
  letter-spacing: -0.1px;
  line-height: 18px;
`;
