import styled from "@emotion/styled";
import { PT_TEXT_COLOR } from "src/components/common/theme";

export const ViewCollection = styled.div`
  padding: 40px;
`;

export const CollectionHero = styled.div`
  align-items: flex-start; /* top aligns collection name with action buttons */
  column-gap: 40px;
  display: flex;

  /* collection name */

  h3 {
    color: ${PT_TEXT_COLOR};
    flex: 1;
    letter-spacing: -0.23px;
    margin: 3.5px 0 0; /* facilitates the center alignment of single-line collection name and top alignment of the first line of a multi-line collection name (line height at 25px) with action buttons (height 32px) */
  }
`;

export const CollectionDetail = styled.div`
  align-items: flex-start;
  display: grid;
  grid-template-areas: "description . metadata"; /* grid areas for collection description and metadata with a null cell token (unnamed area) for the allocation of a gutter between the two columns */
  grid-template-columns: 8fr 1fr 8fr; /* grid columns for collection description and metadata (with 1fr allocated to column separation) */
  margin: 16px 0 44px;
`;
