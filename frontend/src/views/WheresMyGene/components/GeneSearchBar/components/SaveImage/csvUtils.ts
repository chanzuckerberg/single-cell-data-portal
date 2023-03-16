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
import { ChartProps } from "../../../HeatMap/hooks/common/types";
import { deserializeCellTypeMetadata } from "../../../HeatMap/utils";

interface CsvMetadata {
  name: string;
  compareValueName: string;
  viewId: string;
  total_count: number;
}

export function csvHeaders(
  compare: CompareId | undefined,
  availableFilters: Partial<FilterDimensions>,
  availableOrganisms: OntologyTerm[],
  selectedFilters: State["selectedFilters"],
  selectedOrganismId: string | null
) {
  const { datasets, disease_terms, self_reported_ethnicity_terms, sex_terms } =
    availableFilters;

  const output: string[][] = [];

  // Metadata as comments
  output.push([`# ${new Date().toString()}`]);

  // Dataset
  output.push(["# Dataset"]);
  output.push([
    `${
      datasets
        ?.filter((option) => {
          return selectedFilters.datasets.includes(option.id);
        })
        .map((selected) => selected.label)
        .join(", ") || ""
    }`,
  ]);

  // Disease
  output.push(["# Disease"]);
  output.push([
    `${
      disease_terms
        ?.filter((option) => {
          return selectedFilters.diseases.includes(option.id);
        })
        .map((selected) => selected.name)
        .join(", ") || ""
    }`,
  ]);

  // Ethnicity
  output.push(["# Self-Reported Ethnicity"]);
  output.push([
    `${
      self_reported_ethnicity_terms
        ?.filter((option) => {
          return selectedFilters.ethnicities.includes(option.id);
        })
        .map((selected) => selected.name)
        .join(", ") || ""
    }`,
  ]);

  // Sex
  output.push(["# Sex"]);
  output.push([
    `${
      sex_terms
        ?.filter((option) => {
          return selectedFilters.sexes.includes(option.id);
        })
        .map((selected) => selected.name)
        .join(", ") || ""
    }`,
  ]);

  // Organism
  output.push(["# Organism"]);
  output.push([
    `${
      availableOrganisms.find((organism) => organism.id === selectedOrganismId)
        ?.name
    }`,
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

export function csvGeneExpressionRow(
  metadata: CsvMetadata,
  tissue: string,
  allChartProps: { [tissue: string]: ChartProps },
  geneName: string,
  compare: CompareId | undefined
) {
  const { total_count, name, compareValueName, viewId } = metadata;

  const geneExpression = allChartProps[tissue].chartData.find(
    (value) => value.id === `${viewId}-${geneName}`
  );

  if (!compare) {
    return [
      tissue,
      name,
      total_count,
      Number((geneExpression?.tissuePercentage || 0) * 100).toFixed(2) + "%",
      geneName,
      geneExpression?.meanExpression ?? "",
      geneExpression?.scaledMeanExpression ?? "",
      geneExpression?.expressedCellCount ?? "",
    ];
  } else {
    return [
      tissue,
      name,
      total_count,
      Number((geneExpression?.tissuePercentage || 0) * 100).toFixed(2) + "%",
      compareValueName,
      geneName,
      geneExpression?.meanExpression ?? "",
      geneExpression?.scaledMeanExpression ?? "",
      geneExpression?.expressedCellCount ?? "",
    ];
  }
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
    const { id, name, total_count, optionId, viewId, tissue } =
      deserializeCellTypeMetadata(cellTypeMetaData);

    console.log("tissue: " + tissue);
    console.log("id: " + id);
    console.log("name: " + name);
    console.log("total_count: " + total_count);
    console.log("optionId: " + optionId);
    console.log("viewId: " + viewId);
    console.log("---------------------------------------------------");

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

  console.log(cellTypeIdMapping);

  return cellTypeIdMapping;
}
