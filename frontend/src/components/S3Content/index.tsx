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

interface S3File {
  key: string;
  size: number;
  lastModified: string;
}

interface S3ContentProps {
  /** The S3 bucket name (e.g., "gene-expression-assets-public-prod") */
  bucket: string;
  /** AWS region (e.g., "us-west-2"). Defaults to "us-west-2" */
  region?: string;
  /** Optional prefix to filter files (e.g., "data/") */
  prefix?: string;
  /** Base URL for download links. If not provided, uses the S3 URL */
  downloadBaseUrl?: string;
}

/**
 * S3Content - Lists files from a public S3 bucket with download links
 *
 * Usage in MDX:
 * <S3Content bucket="gene-expression-assets-public-prod" />
 * <S3Content bucket="gene-expression-assets-public-prod" prefix="data/" downloadBaseUrl="https://ge-data.cellxgene.cziscience.com" />
 */
const S3Content = ({
  bucket,
  region = "us-west-2",
  prefix = "",
  downloadBaseUrl,
}: S3ContentProps): ReactElement => {
  const [files, setFiles] = useState<S3File[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const fetchBucketListing = async () => {
      try {
        setLoading(true);
        setError(null);

        const s3Url = `https://${bucket}.s3.${region}.amazonaws.com/?list-type=2${prefix ? `&prefix=${encodeURIComponent(prefix)}` : ""}`;

        const response = await fetch(s3Url);

        if (!response.ok) {
          throw new Error(
            `Failed to fetch bucket listing: ${response.status} ${response.statusText}`
          );
        }

        const xml = await response.text();
        const parser = new DOMParser();
        const doc = parser.parseFromString(xml, "application/xml");

        // Check for parsing errors
        const parseError = doc.querySelector("parsererror");
        if (parseError) {
          throw new Error("Failed to parse S3 response");
        }

        // Extract file information from XML
        const contents = doc.querySelectorAll("Contents");
        const fileList: S3File[] = [];

        contents.forEach((item) => {
          const key = item.querySelector("Key")?.textContent || "";
          const size = parseInt(
            item.querySelector("Size")?.textContent || "0",
            10
          );
          const lastModified =
            item.querySelector("LastModified")?.textContent || "";

          // Skip "folder" entries (keys ending with /)
          if (key && !key.endsWith("/")) {
            fileList.push({ key, size, lastModified });
          }
        });

        setFiles(fileList);
      } catch (err) {
        setError(
          err instanceof Error ? err.message : "Failed to fetch bucket listing"
        );
      } finally {
        setLoading(false);
      }
    };

    fetchBucketListing();
  }, [bucket, region, prefix]);

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

  const baseUrl =
    downloadBaseUrl || `https://${bucket}.s3.${region}.amazonaws.com`;

  return (
    <Container>
      <TableContainer>
        <StyledTable>
          <thead>
            <tr>
              <th>File</th>
              <th>Size</th>
              <th>Last Modified</th>
            </tr>
          </thead>
          <tbody>
            {files.map((file) => (
              <tr key={file.key}>
                <td>
                  <DownloadLink
                    href={`${baseUrl}/${file.key}`}
                    target="_blank"
                    rel="noopener noreferrer"
                    download
                  >
                    {getFileName(file.key)}
                  </DownloadLink>
                </td>
                <td>{formatFileSize(file.size)}</td>
                <td>{formatDate(file.lastModified)}</td>
              </tr>
            ))}
          </tbody>
        </StyledTable>
      </TableContainer>
    </Container>
  );
};

function getFileName(key: string): string {
  // Extract just the filename from the full key path
  const parts = key.split("/");
  return parts[parts.length - 1];
}

function formatFileSize(bytes: number): string {
  if (bytes === 0) return "0 B";

  const units = ["B", "KB", "MB", "GB", "TB"];
  const k = 1024;
  const i = Math.floor(Math.log(bytes) / Math.log(k));

  return `${parseFloat((bytes / Math.pow(k, i)).toFixed(1))} ${units[i]}`;
}

function formatDate(isoDate: string): string {
  if (!isoDate) return "";

  const date = new Date(isoDate);
  return date.toLocaleDateString("en-US", {
    year: "numeric",
    month: "short",
    day: "numeric",
  });
}

export default S3Content;
