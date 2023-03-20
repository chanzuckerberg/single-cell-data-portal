import React from "react";
import { GeneHeader, GeneInfoWrapper } from "./style";
import { useGeneInfo } from "src/common/queries/wheresMyGene";

export interface GeneInfoBarProps {
  geneInfoGene: string;
}

function GeneInfoSideBar({
  geneInfoGene,
}: GeneInfoBarProps): JSX.Element | null {
  const { data } = useGeneInfo(geneInfoGene);
  console.log(data);

  return (
    <GeneInfoWrapper
      id="geneinfo_wrapper"
      data-testid={`${geneInfoGene}-gene-info`}
    >
      <GeneHeader
        style={{
          position: "absolute",
        }}
        data-testid="gene-info-header"
      >
        Gene Info
      </GeneHeader>
    </GeneInfoWrapper>
  );
}

export default React.memo(GeneInfoSideBar);
