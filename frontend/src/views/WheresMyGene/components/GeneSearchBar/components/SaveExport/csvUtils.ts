import {
  FilterDimensions,
  OntologyTerm,
  COMPARE_OPTION_ID_FOR_AGGREGATED,
} from "src/common/queries/wheresMyGene";
import {
  CompareId,
  getCompareOptionNameById,
} from "src/views/WheresMyGene/common/constants";
import { State } from "src/views/WheresMyGene/common/store";
import { Props } from ".";
import { deserializeCellTypeMetadata } from "../../../HeatMap/utils";
import { generateAndCopyShareUrl } from "../ShareButton/utils";
import { ChartProps } from "src/views/WheresMyGene/common/types";
import { ChartFormat } from "../../../HeatMap/components/Chart/components/Chart/hooks/utils";

const NO_SELECTION_STRING = "No selection";

interface CsvMetadata {
  name: string;
  compareValueName: string;
  viewId: string;
  total_count: number;
}

const UNDERLYING_DATA_CHANGE_MESSAGE =
  "We regularly expand our single cell data corpus to improve results. Downloaded data and figures may differ in the future.";

export function csvHeaders({
  compare,
  availableFilters,
  availableOrganisms,
  selectedFilters,
  selectedGenes,
  selectedOrganismId,
  selectedTissues,
}: {
  compare: CompareId | undefined;
  availableFilters: Partial<FilterDimensions>;
  availableOrganisms: OntologyTerm[] | null | undefined;
  selectedFilters: State["selectedFilters"];
  selectedGenes: Props["selectedGenes"];
  selectedOrganismId: string | null;
  selectedTissues: Props["selectedTissues"];
}) {
  const {
    datasets,
    disease_terms,
    self_reported_ethnicity_terms,
    sex_terms,
    publication_citations,
  } = availableFilters;

  const output: string[][] = [];

  // Metadata as comments
  output.push([`# Created At: ${new Date().toString()}`]);

  // Data change message
  output.push([`# ${UNDERLYING_DATA_CHANGE_MESSAGE}`]);

  // Share URL
  output.push([
    `# Link Generated: ${generateAndCopyShareUrl({
      compare,
      filters: selectedFilters,
      organism: selectedOrganismId,
      tissues: selectedTissues,
      genes: selectedGenes,
      copyToClipboard: false,
    })}`,
  ]);

  // Dataset
  output.push([
    `# Dataset Filter Values: ${
      datasets
        ?.filter((option) => {
          return selectedFilters.datasets.includes(option.id);
        })
        .map((selected) => selected.label)
        .join(", ") || NO_SELECTION_STRING
    }`,
  ]);

  // Disease
  output.push([
    `# Disease Filter Values: ${
      disease_terms
        ?.filter((option) => {
          return selectedFilters.diseases.includes(option.id);
        })
        .map((selected) => selected.name)
        .join(", ") || NO_SELECTION_STRING
    }`,
  ]);

  // Ethnicity
  output.push([
    `# Self-Reported Ethnicity Filter Values: ${
      self_reported_ethnicity_terms
        ?.filter((option) => {
          return selectedFilters.ethnicities.includes(option.id);
        })
        .map((selected) => selected.name)
        .join(", ") || NO_SELECTION_STRING
    }`,
  ]);

  // Publication
  output.push([
    `# Publication Filter Values: ${
      publication_citations
        ?.filter((option) => {
          return selectedFilters.publications.includes(option.id);
        })
        .map((selected) => selected.name)
        .join(", ") || NO_SELECTION_STRING
    }`,
  ]);

  // Sex
  output.push([
    `# Sex Filter Values: ${
      sex_terms
        ?.filter((option) => {
          return selectedFilters.sexes.includes(option.id);
        })
        .map((selected) => selected.name)
        .join(", ") || NO_SELECTION_STRING
    }`,
  ]);

  // Organism
  output.push([
    `# Organism Filter Value: ${availableOrganisms?.find(
      (organism) => organism.id === selectedOrganismId
    )?.name}`,
  ]);

  // Column Names
  if (compare) {
    const compareOptionName = getCompareOptionNameById(compare);
    output.push([
      "Tissue",
      "Cell Type",
      "Cell Count",
      "Tissue Composition",
      compareOptionName,
      "Gene Symbol",
      "Expression",
      "Expression, Scaled",
      "Number of Cells Expressing Genes",
    ]);
  } else {
    output.push([
      "Tissue",
      "Cell Type",
      "Cell Count",
      "Tissue Composition",
      "Gene Symbol",
      "Expression",
      "Expression, Scaled",
      "Number of Cells Expressing Genes",
    ]);
  }

  return output;
}

export function csvGeneExpressionRow({
  metadata,
  tissueName,
  allChartProps,
  geneName,
  compare,
}: {
  metadata: CsvMetadata;
  tissueName: string;
  allChartProps: { [tissue: string]: ChartProps };
  geneName: string;
  compare: CompareId | undefined;
}) {
  const { total_count, name, compareValueName, viewId } = metadata;

  const geneExpression = allChartProps[tissueName].chartData.find(
    (value) => value.id === `${viewId}-${geneName}`
  );

  if (compare) {
    return [
      tissueName,
      name,
      total_count,
      getTissuePercentage(geneExpression),
      compareValueName,
      geneName,
      geneExpression?.meanExpression ?? "",
      geneExpression?.scaledMeanExpression ?? "",
      geneExpression?.expressedCellCount ?? "",
    ];
  } else {
    return [
      tissueName,
      name,
      total_count,
      getTissuePercentage(geneExpression),
      geneName,
      geneExpression?.meanExpression ?? "",
      geneExpression?.scaledMeanExpression ?? "",
      geneExpression?.expressedCellCount ?? "",
    ];
  }
}

function getTissuePercentage(geneExpression: ChartFormat | undefined) {
  if (!geneExpression) return "";

  return Number((geneExpression.tissuePercentage || 0) * 100).toFixed(2) + "%";
}

export function buildCellTypeIdToMetadataMapping(
  tissue: string,
  allChartProps: {
    [tissue: string]: ChartProps;
  }
) {
  const cellTypeIdMapping: {
    [cellTypeId: string]: CsvMetadata[];
  } = {};

  let currentCellTypeName = "";
  for (const cellTypeMetaData of allChartProps[
    tissue
  ].cellTypeMetadata.reverse()) {
    const { id, name, total_count, optionId, viewId } =
      deserializeCellTypeMetadata(cellTypeMetaData);

    if (!cellTypeIdMapping[id]) {
      cellTypeIdMapping[id] = [];
      currentCellTypeName = name;
    }

    cellTypeIdMapping[id].push({
      name: currentCellTypeName,
      compareValueName:
        optionId === COMPARE_OPTION_ID_FOR_AGGREGATED
          ? COMPARE_OPTION_ID_FOR_AGGREGATED
          : name.trim(),
      viewId: viewId,
      total_count: total_count,
    });
  }

  return cellTypeIdMapping;
}

// Gets the date in mmddyy format
export function getCurrentDate() {
  const today = new Date();
  const month = (today.getMonth() + 1).toString().padStart(2, "0");
  const day = today.getDate().toString().padStart(2, "0");
  const year = today.getFullYear().toString().slice(-2);

  return month + day + year;
}
