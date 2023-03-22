import React from "react";
import {
  GeneName,
  GeneSummary,
  GeneInfoWrapper,
  GeneSynonyms,
  GeneSynonymsLabel,
  GeneSynonymsWrapper,
  GeneUrl,
  StyledCallout,
} from "./style";
import { useGeneInfo } from "src/common/queries/wheresMyGene";

export interface GeneInfoBarProps {
  geneInfoGene: string;
}

function GeneInfoSideBar({
  geneInfoGene,
}: GeneInfoBarProps): JSX.Element | null {
  const { data, isLoading } = useGeneInfo(geneInfoGene);

  if (isLoading) return <GeneSummary>Loading...</GeneSummary>;

  return (
    <GeneInfoWrapper
      id="gene-info-wrapper"
      data-test-id={`${geneInfoGene}-gene-info`}
    >
      {!data ? (
        <>
          <StyledCallout autoDismiss={false} intent={"error"}>
            Sorry, this gene could not be found on NCBI.
          </StyledCallout>

          <GeneUrl>
            <a
              href={`https://www.google.com/search?q=${geneInfoGene}%20gene`}
              target="_blank"
              rel="noreferrer noopener"
            >
              Search on Google
            </a>
          </GeneUrl>
        </>
      ) : (
        <>
          <GeneName data-test-id="gene-info-header">{data.name}</GeneName>

          {data.show_warning_banner && (
            <>
              <StyledCallout autoDismiss={false} intent={"warning"}>
                NCBI didn't return an exact match for this gene.
              </StyledCallout>
            </>
          )}

          <GeneSummary>{data.summary}</GeneSummary>

          <GeneSynonymsWrapper>
            <GeneSynonymsLabel>Synonyms</GeneSynonymsLabel>
            <GeneSynonyms>{data.synonyms.join(", ")}</GeneSynonyms>
          </GeneSynonymsWrapper>

          <GeneUrl>
            <a href={data.ncbi_url} target="_blank" rel="noreferrer noopener">
              View on NCBI
            </a>
          </GeneUrl>
        </>
      )}
    </GeneInfoWrapper>
  );
}

export default React.memo(GeneInfoSideBar);
