import React from "react";
import {
  GeneHeader,
  GeneSummary,
  GeneInfoWrapper,
  GeneSynonyms,
  GeneSynonymsLabel,
  GeneSynonymsWrapper,
  GeneUrl,
} from "./style";
import { useGeneInfo } from "src/common/queries/wheresMyGene";

export interface GeneInfoBarProps {
  geneInfoGene: string;
}

function GeneInfoSideBar({
  geneInfoGene,
}: GeneInfoBarProps): JSX.Element | null {
  const { data, isLoading } = useGeneInfo(geneInfoGene);

  if (!data || isLoading) return <>Loading...</>;

  console.log(data);

  return (
    <GeneInfoWrapper
      id="gene-info-wrapper"
      data-test-id={`${geneInfoGene}-gene-info`}
    >
      <GeneHeader data-test-id="gene-info-header">{data.name}</GeneHeader>

      <GeneSummary>{data.summary}</GeneSummary>

      <GeneSynonymsWrapper>
        <GeneSynonymsLabel>Synonyms</GeneSynonymsLabel>
        <GeneSynonyms>{data.synonyms.join(", ")}</GeneSynonyms>
      </GeneSynonymsWrapper>

      <GeneUrl>
        <a href={data.ncbi_url}>View on NCBI</a>
      </GeneUrl>
    </GeneInfoWrapper>
  );
}

export default React.memo(GeneInfoSideBar);
