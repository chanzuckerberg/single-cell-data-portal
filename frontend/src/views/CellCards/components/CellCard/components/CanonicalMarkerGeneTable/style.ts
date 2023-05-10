import styled from "@emotion/styled";
import { fontBodyS, fontHeaderM, getColors } from "czifui";

export const TableTitle = styled.div`
  ${fontHeaderM}
  font-weight: 600;
  margin-bottom: 8px;
`;

export const PublicationLinkWrapper = styled.div`
  display: flex;
  flex-direction: column;
`;

export const TableTitleWrapper = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  margin-top: 48px;
`;

export const WmgLink = styled.a`
  ${fontBodyS}
  font-weight: 500;
  ${(props) => {
    const colors = getColors(props);
    return `color: ${colors?.primary[400]}`;
  }}
`;
