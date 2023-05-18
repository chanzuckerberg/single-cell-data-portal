import styled from "@emotion/styled";
import { fontHeaderM, fontBodyXxs, fontBodyS, getColors } from "czifui";

export const TableTitle = styled.div`
  ${fontHeaderM}
`;

export const PublicationLinkWrapper = styled.div`
  display: flex;
  column-gap: 2px;
`;

export const TableTitleWrapper = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  margin-top: 48px;
  margin-bottom: 8px;
`;

export const TableUnavailableContainer = styled("div")`
  margin-top: 16px;
  background: #f8f8f8;

  width: 100%;

  height: 120px;

  display: flex;
  flex-direction: column;

  justify-content: center;
  text-align: center;
  border-radius: 8px;
`;

export const TableUnavailableHeader = styled("span")`
  ${fontBodyS}
  color: black;
  font-weight: 500;
`;

export const TableUnavailableDescription = styled("span")`
  ${fontBodyXxs}
  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;
