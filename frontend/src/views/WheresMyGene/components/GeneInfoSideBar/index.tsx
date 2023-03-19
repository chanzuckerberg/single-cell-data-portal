import React from "react";

export interface GeneInfoBarProps {
  geneInfoGene: string;
}

function GeneInfoSideBar({
  geneInfoGene,
}: GeneInfoBarProps): JSX.Element | null {
  return <div>{geneInfoGene}</div>;
}

export default React.memo(GeneInfoSideBar);
