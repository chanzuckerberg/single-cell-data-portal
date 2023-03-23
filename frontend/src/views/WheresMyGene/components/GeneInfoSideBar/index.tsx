import React from "react";
import {
  GeneName,
  GeneSummary,
  GeneSynonyms,
  GeneSynonymsLabel,
  GeneSynonymsWrapper,
  GeneUrl,
  StyledCallout,
} from "./style";
import { useGeneInfo } from "src/common/queries/wheresMyGene";
import { RightSidebarProperties } from "../RightSideBar";

export interface GeneInfoBarProps extends RightSidebarProperties {
  geneInfoGene: string;
}

function GeneInfoSideBar({
  geneInfoGene,
}: GeneInfoBarProps): JSX.Element | null {
  const { data, isLoading } = useGeneInfo(geneInfoGene);

  if (isLoading) return <GeneSummary>Loading...</GeneSummary>;

  return (
    <div id="gene-info-wrapper" data-test-id={`${geneInfoGene}-gene-info`}>
      {!data ? (
        <>
          <StyledCallout autoDismiss={false} intent={"error"}>
            Sorry, this gene could not be found on NCBI.
          </StyledCallout>

          <GeneUrl
            href={`https://www.google.com/search?q=${geneInfoGene}%20gene`}
            target="_blank"
            rel="noreferrer noopener"
          >
            Search on Google
          </GeneUrl>
        </>
      ) : (
        <>
          <GeneName data-test-id="gene-info-header">{data.name}</GeneName>

          {data.show_warning_banner && (
            <>
              <StyledCallout autoDismiss={false} intent={"warning"}>
                NCBI didn&apos;t return an exact match for this gene.
              </StyledCallout>
            </>
          )}

          <GeneSummary>{data.summary}</GeneSummary>

          <GeneSynonymsWrapper>
            <GeneSynonymsLabel>Synonyms</GeneSynonymsLabel>
            <GeneSynonyms>{data.synonyms.join(", ")}</GeneSynonyms>
          </GeneSynonymsWrapper>

          <GeneUrl
            href={data.ncbi_url}
            target="_blank"
            rel="noreferrer noopener"
          >
            View on NCBI
          </GeneUrl>
        </>
      )}
    </div>
  );
}

export default React.memo(GeneInfoSideBar);
