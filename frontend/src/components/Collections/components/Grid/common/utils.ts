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
        cell_count: (cell_count || 0) + Number(metadata.cell_count),
        disease: [...disease, ...metadata.disease],
        organism: [...organism, ...metadata.organism],
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
    assay: caseInsensitiveSort(unique(result.assay)),
    cell_count: result.cell_count,
    disease: caseInsensitiveSort(unique(result.disease)),
    organism: caseInsensitiveSort(unique(result.organism)),
    tissue: caseInsensitiveSort(unique(result.tissue)),
  };
}, hashFn);

function hashFn(datasets: Dataset[]): string {
  return datasets.reduce((acc, { id, cell_count }) => {
    // (thuang): We need `cell_count` as well to aggregate the metadata
    // when they become available
    return acc + id + cell_count;
  }, "");
}

function extractDatasetMetadata(dataset: Dataset) {
  const { tissue, organism, assay, disease, cell_count } = dataset;

  return {
    assay: uniqueOntologies(assay),
    cell_count,
    disease: uniqueOntologies(disease),
    organism: uniqueOntologies(organism),
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

function caseInsensitiveSort(items: string[]): string[] {
  // (thuang): Intl.Collator is used for performance
  // https://stackoverflow.com/a/52369951
  const collator = new Intl.Collator("en", {
    numeric: true,
    sensitivity: "base",
  });

  return items.sort((a, b) => collator.compare(a, b));
}
