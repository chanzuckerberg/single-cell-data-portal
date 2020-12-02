import { Dataset, Ontology } from "src/common/entities";

export function aggregateDatasetsMetadata(datasets: Dataset[]) {
  const datasetsMetadata = datasets.map(extractDatasetMetadata);

  // (thuang): Merge all datasets' metadata
  const result = datasetsMetadata.reduce(
    (prev, metadata) => {
      const { assay, disease, organism, tissue, cell_count } = prev;

      return {
        assay: [...assay, ...metadata.assay],
        cell_count: cell_count + Number(metadata.cell_count),
        disease: [...disease, ...metadata.disease],
        organism: [...organism, metadata.organism],
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
}

function extractDatasetMetadata(dataset: Dataset) {
  const { tissue, organism, assay, disease, cell_count } = dataset;

  return {
    assay: uniqueOntologies(assay),
    cell_count,
    disease: uniqueOntologies(disease),
    organism: organism.label,
    tissue: uniqueOntologies(tissue),
  };
}

function uniqueOntologies(ontologies: Ontology[]) {
  return unique(
    ontologies
      .map((ontology) => {
        return ontology.label;
      })
      .sort()
  );
}

function unique(items: string[]) {
  return Array.from(new Set(items));
}
