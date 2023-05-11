import styled from "@emotion/styled";
import { fontBodyS, CommonThemeProps, getColors } from "czifui";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { TOP_PADDING_PX } from "../CellCard/style";

export const CellCardsSidebarWrapper = styled.div`
  padding: ${TOP_PADDING_PX}px 0px 32px 40px;
`;

export const StickyWrapper = styled.div`
  display: block;
  position: sticky;
  top: ${HEADER_HEIGHT_PX + TOP_PADDING_PX}px;
`;

export const TableOfContents = styled.div`
  display: flex;
  flex-direction: column;
`;

interface StyledJumpLinkProps extends CommonThemeProps {
  isActive: boolean;
}

export const StyledJumpLink = styled.a<StyledJumpLinkProps>`
  ${fontBodyS}
  ${(props) => {
    const isActive = props.isActive;
    const colors = getColors(props);
    if (isActive) {
      return `
        color: black;
        font-weight: 600;
        border-left: 2px solid #0073ff; // Unsure what color this blue is in SDS
      `;
    } else {
      return `
        color: ${colors?.gray[500]};
        font-weight: 500;
        border-left: 2px solid ${colors?.gray[200]};
      `;
    }
  }}
  padding: 6px 0 6px 16px;

  /* animation */
  transition: border-color 0.1s ease-in-out;
`;
