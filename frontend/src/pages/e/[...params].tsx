import { useRouter } from "next/router";
import { useEffect, useState } from "react";
import styled from "@emotion/styled";
import configs from "../../configs/configs";

const Container = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  width: 100%;
  height: 100vh;
  background-color: #f5f5f5;
`;

const Card = styled.div`
  background: white;
  border-radius: 8px;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
  padding: 30px;
  max-width: 600px;
  width: 90%;
  text-align: center;
`;

const Title = styled.h1`
  color: #333;
  margin-top: 0;
  margin-bottom: 20px;
`;

const ErrorMessage = styled.div`
  color: #d32f2f;
  padding: 15px;
  background-color: #ffebee;
  border-left: 4px solid #d32f2f;
  border-radius: 4px;
  margin: 20px 0;
`;

const CellxgeneFrame = styled.iframe`
  width: 100%;
  height: 100vh;
  border: none;
  margin: 0;
  padding: 0;
`;

/**
 * Route handler for explorer URLs: /e/[id].h5ad
 *
 * Launches cellxgene and displays it embedded in an iframe.
 */
export default function ExplorerRoute() {
  const router = useRouter();
  const { params } = router.query;
  const [cellxgeneUrl, setCellxgeneUrl] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    if (!params || params.length === 0) {
      return;
    }

    // Reset state when params change (new dataset selected)
    setCellxgeneUrl(null);
    setError(null);
    setLoading(true);

    const launchCellxgene = async () => {
      try {
        // Extract dataset ID from URL (remove .h5ad extension if present)
        let datasetId = params[0];
        if (datasetId.endsWith(".h5ad")) {
          datasetId = datasetId.replace(".h5ad", "");
        }

        // Call cellxgene-launch endpoint
        const apiUrl = `${configs.API_URL}/dp/v1/datasets/${datasetId}/cellxgene-launch`;

        const response = await fetch(apiUrl);

        if (response.ok) {
          const metadata = await response.json();
          const url = metadata.cellxgene_url;

          // Load cellxgene in iframe after a delay to ensure all resources are ready
          // Cellxgene needs time to fully initialize before serving all static resources
          setTimeout(() => {
            setCellxgeneUrl(url);
            setLoading(false);
          }, 3000);
        } else {
          setError(`Failed to launch cellxgene (${response.status})`);
          setLoading(false);
        }
      } catch (err) {
        const message = err instanceof Error ? err.message : "Unknown error";
        console.error("Error launching cellxgene:", err);
        setError(message);
        setLoading(false);
      }
    };

    launchCellxgene();
  }, [params]);

  if (error) {
    return (
      <Container>
        <Card>
          <Title>Error Loading Cellxgene</Title>
          <ErrorMessage>{error}</ErrorMessage>
        </Card>
      </Container>
    );
  }

  if (loading || !cellxgeneUrl) {
    return (
      <Container>
        <Card>
          <Title>Loading Cellxgene</Title>
          <p style={{ color: "#666" }}>
            Please wait while we prepare your data...
          </p>
        </Card>
      </Container>
    );
  }

  return <CellxgeneFrame src={cellxgeneUrl} title="Cellxgene Visualization" />;
}
