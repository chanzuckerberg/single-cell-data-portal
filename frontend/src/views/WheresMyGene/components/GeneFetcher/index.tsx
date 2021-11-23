import { useEffect } from "react";
import { Gene, RawGeneExpression } from "../../common/types";
import GENE_EXPRESSIONS from "../../mocks/lung_tissue_cond.json";

interface Props {
  name: string;
  onSuccess: (gene: RawGeneExpression) => void;
  onError: (id: string) => void;
  minDelayMS?: number;
  maxDelayMS?: number;
  fetchedGenes: Gene[];
}

const renderedGenes: string[] = [];

export default function GeneFetcher({
  name,
  onSuccess,
  onError,
  minDelayMS = 200,
  maxDelayMS = 500,
}: Props): null {
  useEffect(() => {
    if (renderedGenes.includes(name)) return;

    renderedGenes.push(name);

    const delayMS = Math.floor(minDelayMS + Math.random() * maxDelayMS);

    setTimeout(() => {
      const geneData = (
        GENE_EXPRESSIONS as unknown as RawGeneExpression[]
      ).find((gene) => gene.gene_name === name);

      if (geneData) {
        onSuccess(geneData);
      } else {
        onError(name);
      }
    }, delayMS);
  }, [name, onSuccess, onError, minDelayMS, maxDelayMS]);

  return null;
}
