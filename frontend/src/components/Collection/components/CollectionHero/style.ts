import { PT_TEXT_COLOR } from "src/components/common/theme";
import styled from "styled-components";

export const CollectionHero = styled.div`
  align-items: flex-start; /* top aligns collection name with action buttons */
  display: flex;
  column-gap: 40px;

  /* collection name */
  h3 {
    color: ${PT_TEXT_COLOR};
    flex: 1;
    letter-spacing: -0.23px;
    margin: 2.5px 0 0; /* facilitates the top alignment but also center alignment of collection name (line height at 25px) with action buttons (height 30px) */
  }
`;
