import { BLUE, PT_TEXT_COLOR } from "src/components/common/theme";
import styled from "styled-components";

// Collection Title Column
export const CollectionTitleText = styled.a`
  color: ${BLUE.C};
  font-size: 14px;
  font-weight: 500;
  line-height: 18px;
  letter-spacing: -0.1px;
`;

export const DOILink = styled.a`
  color: ${BLUE.C};
  font-style: normal;
  font-weight: normal;
  font-size: 12px;
  line-height: 15px;
`;

export const ContactText = styled.div`
  color: ${PT_TEXT_COLOR};
  font-style: normal;
  font-weight: normal;
  font-size: 12px;
  line-height: 15px;
`;
