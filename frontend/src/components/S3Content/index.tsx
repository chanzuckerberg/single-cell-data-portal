import styled from "@emotion/styled";
import { ReactElement, useEffect, useState } from "react";

const Container = styled.div`
  margin: 16px 0;
`;

const LoadingText = styled.div`
  color: #666;
  font-style: italic;
`;

const ErrorText = styled.div`
  color: #d32f2f;
  padding: 12px;
  background-color: #ffebee;
  border-radius: 4px;
`;

const TableContainer = styled.div`
  overflow-x: auto;
  border: 1px solid #e0e0e0;
  border-radius: 4px;
`;

const StyledTable = styled.table`
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;

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

const DownloadLink = styled.a`
  color: #1976d2;
  text-decoration: none;
  font-weight: 500;

  &:hover {
    text-decoration: underline;
  }
`;

interface S3ContentProps {
  /** Base URL for the CloudFront distribution (e.g., "https://ge-data.cellxgene.cziscience.com") */
  downloadBaseUrl: string;
  /** Name of the JSON manifest file. Defaults to "expression-summary-files.json" */
  manifestFile?: string;
}

async function fetchManifest(
  downloadBaseUrl: string,
  manifestFile: string
): Promise<string[]> {
  const manifestUrl = `${downloadBaseUrl}/${manifestFile}`;
  console.log("S3Content: fetching manifest from", manifestUrl);

  const response = await fetch(manifestUrl);
  console.log(
    "S3Content: response status",
    response.status,
    response.statusText
  );
  console.log("S3Content: response headers", response.headers);

  if (!response.ok) {
    throw new Error(`Failed to fetch manifest: ${response.statusText}`);
  }

  const text = await response.text();
  console.log("S3Content: raw response text", text);

  const data = JSON.parse(text);
  console.log("S3Content: manifest data", data);
  return data.files || [];
}

/**
 * S3Content - Lists files from a JSON manifest with download links
 *
 * Usage in MDX:
 * <S3Content downloadBaseUrl="https://ge-data.cellxgene.cziscience.com" />
 * <S3Content downloadBaseUrl="https://ge-data.cellxgene.cziscience.com" manifestFile="expression-summary-files.json" />
 */
const S3Content = ({
  downloadBaseUrl,
  manifestFile = "expression-summary-files.json",
}: S3ContentProps): ReactElement => {
  const [files, setFiles] = useState<string[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  console.log("S3Content: Component rendering", {
    downloadBaseUrl,
    manifestFile,
  });

  useEffect(() => {
    console.log("S3Content: useEffect triggered", {
      downloadBaseUrl,
      manifestFile,
    });
    setLoading(true);
    setError(null);

    fetchManifest(downloadBaseUrl, manifestFile)
      .then((files) => {
        console.log("S3Content: files fetched", files);
        setFiles(files);
      })
      .catch((err) => {
        console.error("S3Content: fetch error", err);
        setError(
          err instanceof Error ? err.message : "Failed to fetch file manifest"
        );
      })
      .finally(() => setLoading(false));
  }, [downloadBaseUrl, manifestFile]);

  if (loading) {
    return (
      <Container>
        <LoadingText>Loading file listing...</LoadingText>
      </Container>
    );
  }

  if (error) {
    return (
      <Container>
        <ErrorText>Error: {error}</ErrorText>
      </Container>
    );
  }

  if (files.length === 0) {
    return (
      <Container>
        <LoadingText>No files found</LoadingText>
      </Container>
    );
  }

  return (
    <Container>
      <TableContainer>
        <StyledTable>
          <thead>
            <tr>
              <th>File</th>
            </tr>
          </thead>
          <tbody>
            {files.map((filename) => (
              <tr key={filename}>
                <td>
                  <DownloadLink
                    href={`${downloadBaseUrl}/${filename}`}
                    target="_blank"
                    rel="noopener noreferrer"
                    download
                  >
                    {filename}
                  </DownloadLink>
                </td>
              </tr>
            ))}
          </tbody>
        </StyledTable>
      </TableContainer>
    </Container>
  );
};

export default S3Content;
