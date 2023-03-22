import React from "react";
import {
  GeneName,
  GeneSummary,
  GeneInfoWrapper,
  GeneSynonyms,
  GeneSynonymsLabel,
  GeneSynonymsWrapper,
  GeneUrl,
  WarningBanner,
  InfoButtonWrapper,
} from "./style";
import { useGeneInfo } from "src/common/queries/wheresMyGene";
import { Icon } from "czifui";

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
          <GeneSummary>
            Sorry, this gene could not be found on NCBI.
          </GeneSummary>

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
              <WarningBanner>
                <InfoButtonWrapper>
                  <Icon
                    sdsIcon="exclamationMarkCircle"
                    sdsSize="l"
                    sdsType="static"
                  />
                </InfoButtonWrapper>
                NCBI didn't return an exact match for this gene.
              </WarningBanner>
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
