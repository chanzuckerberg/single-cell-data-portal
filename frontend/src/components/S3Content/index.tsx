import { ReactElement, useEffect, useState } from "react";
import {
  Container,
  DownloadLink,
  ErrorText,
  LoadingText,
  StyledTable,
  TableContainer,
} from "./style";

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
  const response = await fetch(manifestUrl);

  if (!response.ok) {
    throw new Error(`Failed to fetch manifest: ${response.statusText}`);
  }

  const text = await response.text();
  let data;
  try {
    data = JSON.parse(text);
  } catch (err) {
    throw new Error(
      `Manifest file is not valid JSON. Please check the format of "${manifestFile}".` +
        (err instanceof Error ? ` Details: ${err.message}` : "")
    );
  }
  return Array.isArray(data.files) ? data.files : [];
}

/**
 * S3Content - Lists files from a JSON manifest with download links
 *
 * Usage in MDX:
 * <S3Content downloadBaseUrl="https://ge-data.cellxgene.cziscience.com" />
 * <S3Content downloadBaseUrl="https://ge-data.cellxgene.cziscience.com" manifestFile="expression-summary-files.json" />
 */
function S3Content({
  downloadBaseUrl,
  manifestFile = "expression-summary-files.json",
}: S3ContentProps): ReactElement {
  const [files, setFiles] = useState<string[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    setLoading(true);
    setError(null);

    fetchManifest(downloadBaseUrl, manifestFile)
      .then((files) => {
        setFiles(files);
      })
      .catch((err) => {
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
                    download={filename}
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
}

export default S3Content;
