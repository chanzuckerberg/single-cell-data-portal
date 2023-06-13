import React from "react";
import {
  GeneName,
  GeneSummary,
  GeneSynonyms,
  Label,
  GeneSynonymsWrapper,
  GeneUrl,
  OutLinksWrapper,
  StyledCallout,
} from "./style";
import { useGeneInfo } from "src/common/queries/wheresMyGene";
import { RightSidebarProperties } from "../common/RightSideBar";
import { GeneInfo } from "../../views/WheresMyGene/common/types";

const GENE_CARDS_URL = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=";

export interface GeneInfoBarProps extends RightSidebarProperties {
  geneInfoGene: string;
}

function GeneInfoSideBar({
  geneInfoGene,
}: GeneInfoBarProps): JSX.Element | null {
  const { data, isLoading } = useGeneInfo(geneInfoGene);

  if (isLoading) return <GeneSummary>Loading...</GeneSummary>;

  // Component for when gene info cannot be found
  const GeneInfoNoData = (
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
  );

  return (
    <div id="gene-info-wrapper" data-testid={`${geneInfoGene}-gene-info`}>
      {data ? (
        <GeneInfoResult geneInfoGene={geneInfoGene} data={data} />
      ) : (
        GeneInfoNoData
      )}
    </div>
  );
}

function GeneInfoResult({
  geneInfoGene,
  data,
}: {
  geneInfoGene: string;
  data: GeneInfo;
}) {
  const { name, summary, synonyms, show_warning_banner, ncbi_url } = data;

  return (
    <div id="gene-info-wrapper">
      <>
        <GeneName data-testid="gene-info-header">{name}</GeneName>

        {show_warning_banner && (
          <>
            <StyledCallout autoDismiss={false} intent={"warning"}>
              NCBI didn&apos;t return an exact match for this gene.
            </StyledCallout>
          </>
        )}

        <GeneSummary data-testid="gene-info-gene-summary">
          {summary}
        </GeneSummary>

        <GeneSynonymsWrapper data-testid="gene-info-gene-synonyms">
          <Label>Synonyms</Label>
          <GeneSynonyms>{synonyms.join(", ")}</GeneSynonyms>
        </GeneSynonymsWrapper>

        <OutLinksWrapper>
          <GeneUrl
            href={GENE_CARDS_URL + geneInfoGene}
            target="_blank"
            rel="noreferrer noopener"
            data-testid="gene-info-geneCards-link"
          >
            View on GeneCards
          </GeneUrl>

          <span>
            <Label>
              Source:{" "}
              <GeneUrl
                href={ncbi_url}
                target="_blank"
                rel="noreferrer noopener"
                data-testid="gene-info-ncbi-link"
              >
                NCBI
              </GeneUrl>
            </Label>
          </span>
        </OutLinksWrapper>
      </>
    </div>
  );
}

export default React.memo(GeneInfoSideBar);
