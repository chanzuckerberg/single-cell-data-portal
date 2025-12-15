import styled from "@emotion/styled";

export const Container = styled.div`
  margin: 0;
`;

export const LoadingText = styled.div`
  color: #666;
  font-style: italic;
`;

export const ErrorText = styled.div`
  color: #d32f2f;
  padding: 12px;
  background-color: #ffebee;
  border-radius: 4px;
`;

export const TableContainer = styled.div`
  overflow-x: auto;
  border: 1px solid #e0e0e0;
  border-radius: 4px;
`;

export const StyledTable = styled.table`
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
  margin: 0;

  th,
  td {
    padding: 10px 14px;
    text-align: left;
    border-bottom: 1px solid #e0e0e0;
  }

  th {
    background-color: #f5f5f5;
    font-weight: 600;
  }

  tr:last-child td {
    border-bottom: none;
  }

  tr:hover td {
    background-color: #fafafa;
  }
`;

export const DownloadLink = styled.a`
  color: #1976d2;
  text-decoration: none;
  font-weight: 500;

  &:hover {
    text-decoration: underline;
  }
`;
