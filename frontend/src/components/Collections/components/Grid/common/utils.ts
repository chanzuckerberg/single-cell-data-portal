import memoize from "lodash/memoize";
import { Dataset, Ontology } from "src/common/entities";

export const aggregateDatasetsMetadata = memoize(function (
  datasets: Dataset[]
) {
  const datasetsMetadata = datasets.map(extractDatasetMetadata);

  // (thuang): Merge all datasets' metadata
  const result = datasetsMetadata.reduce(
    (prev, metadata) => {
      const { assay, disease, organism, tissue, cell_count } = prev;

      return {
        assay: [...assay, ...metadata.assay],
        cell_count: cell_count + Number(metadata.cell_count),
        disease: [...disease, ...metadata.disease],
        organism: metadata.organism
          ? [...organism, metadata.organism]
          : organism,
        tissue: [...tissue, ...metadata.tissue],
      };
    },
    {
      assay: [] as string[],
      cell_count: 0,
      disease: [] as string[],
      organism: [] as string[],
      tissue: [] as string[],
    }
  );

  return {
    assay: unique(result.assay),
    cell_count: result.cell_count,
    disease: unique(result.disease),
    organism: unique(result.organism),
    tissue: unique(result.tissue),
  };
},
hashFn);

function hashFn(datasets: Dataset[]): string {
  return datasets.reduce((acc, dataset) => {
    return acc + dataset.id;
  }, "");
}

function extractDatasetMetadata(dataset: Dataset) {
  const { tissue, organism, assay, disease, cell_count } = dataset;

  return {
    assay: uniqueOntologies(assay),
    cell_count,
    disease: uniqueOntologies(disease),
    organism: organism?.label,
    tissue: uniqueOntologies(tissue),
  };
}

function uniqueOntologies(ontologies: Ontology[]) {
  if (!ontologies) return [];
  return unique(
    ontologies
      .map((ontology) => {
        return ontology.label;
      })
      .sort()
  );
}

function unique(items: string[]) {
  if (!items) return [];
  return Array.from(new Set(items));
}
