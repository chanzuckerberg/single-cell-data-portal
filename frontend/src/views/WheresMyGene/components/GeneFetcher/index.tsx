import { useEffect } from "react";
import { API } from "src/common/API";
import { DEFAULT_FETCH_OPTIONS } from "src/common/queries/common";
import { apiTemplateToUrl } from "src/common/utils/apiTemplateToUrl";
import { API_URL } from "src/configs/configs";
import {
  CellTypeGeneExpressionSummaryData,
  Gene,
  GeneExpressionSummary,
  RawCellTypeGeneExpressionSummaryData,
} from "../../common/types";

interface Props {
  name: string;
  onSuccess: (gene: GeneExpressionSummary) => void;
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

    fetchGeneData();

    async function fetchGeneData(): Promise<void> {
      const response = await fetch(
        apiTemplateToUrl(API_URL + API.WMG_GENE, { name }),
        DEFAULT_FETCH_OPTIONS
      );

      // DEBUG
      // DEBUG
      // DEBUG
      // (thuang): Local test data
      // const response = await fetch(
      //   `https://wmg-prototype-data-dev-public.s3.amazonaws.com/lung-tissue-10x-human/genes/${name}.json`
      // );

      const expressions = await response.json();

      setTimeout(() => {
        if (expressions) {
          onSuccess({
            cellTypeGeneExpressionSummaries: expressions.map(
              transformCellTypeGeneExpressionSummaryData
            ),
            name,
          });
        } else {
          onError(name);
        }
      }, delayMS);
    }
  }, [name, onSuccess, onError, minDelayMS, maxDelayMS]);

  return null;
}

function transformCellTypeGeneExpressionSummaryData(
  data: RawCellTypeGeneExpressionSummaryData
): CellTypeGeneExpressionSummaryData {
  const { id, pc, me } = data;

  return {
    id,
    meanExpression: me,
    percentage: pc,
  };
}
