import styled from "@emotion/styled";

export const CellCardsSidebarWrapper = styled.div`
  padding: 32px 0px 32px 40px;
`;

export const StickyWrapper = styled.div`
  display: block;
  position: sticky;
  top: 80px;
`;

export const TableOfContents = styled.div`
  display: flex;
  flex-direction: column;
`;

export const StyledJumpLink = styled.a`
  ${isActive}

  font-family: Inter;
  letter-spacing: -0.006em;
  line-height: 20px;
  margin: 0;
  padding: 6px 0 6px 16px;
  font-size: 16px;

  /* animation */
  transition: border-color 0.1s ease-in-out;
`;

function isActive({ isActive }: { isActive: boolean }) {
  if (isActive) {
    return `
      color: black;
      font-weight: 600;
      border-left: 2px solid #0073ff;
    `;
  } else {
    return `
      color: #767676;
      font-weight: 500;
      border-left: 2px solid #EAEAEA;
    `;
  }
}
