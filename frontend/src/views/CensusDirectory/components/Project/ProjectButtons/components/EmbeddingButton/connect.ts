import { useCallback, useEffect, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EmbeddingButtonProps } from "./types";
import { getProjectTier } from "src/views/CensusDirectory/utils";
import { UnionProject } from "../../../types";

// The div contains two lines of the word copy
const NUMBER_OF_EXTRA_LINES = 2;
// There is .5 em of padding on the top and bottom of the snippets
const NUMBER_OF_PADDING_LINES = 1;
// Total amount of padding around the highlighted line
const LINE_HIGHLIGHT_BACKGROUND_PADDING = 8;

const MAINTAINED_PYTHON_NOTEBOOK_LINK =
  "https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_access_maintained_embeddings.html";
const MAINTAINED_R_NOTEBOOK_LINK =
  "https://chanzuckerberg.github.io/cellxgene-census/r/articles/census_access_maintained_embeddings.html";
const HOSTED_PYTHON_NOTEBOOK_LINK =
  "https://chanzuckerberg.github.io/cellxgene-census/notebooks/api_demo/census_embedding.html";

function pythonCodeSnippet(project: UnionProject, uri: string): string {
  const censusVersion = project.census_version;
  const organism = project.experiment_name;
  const measurement = project.measurement_name;

  return project.tier === "maintained"
    ? `import cellxgene_census

census = cellxgene_census.open_soma(census_version="${project.census_version}")
adata = cellxgene_census.get_anndata(
    census,
    organism = "${organism}",
    measurement_name = "${measurement}",
    obs_value_filter = "tissue_general == 'central nervous system'",
    obsm_layers = ["${project.obsm_layer}"]
)`
    : `import cellxgene_census
from cellxgene_census.experimental import get_embedding

embedding_uri = \\
    "${uri}"
census = cellxgene_census.open_soma(census_version="${censusVersion}")

adata = cellxgene_census.get_anndata(
    census,
    organism = "${organism}",
    measurement_name = "${measurement}",
    obs_value_filter = "tissue_general == 'central nervous system'",
)
embeddings = get_embedding("${censusVersion}", embedding_uri, adata.obs["soma_joinid"].to_numpy())
adata.obsm["emb"] = embeddings`;
}

function rCodeSnippet(project: UnionProject): string {
  const censusVersion = project.census_version;
  const organism = project.experiment_name;

  return project.tier === "maintained"
    ? `library("cellxgene.census")
library("Seurat")

census <- open_soma(census_version = "${censusVersion}")
seurat_obj <- get_seurat(
  census,
  organism = "${organism}",
  obs_value_filter = "tissue_general == 'central nervous system'",
  obs_column_names = c("cell_type"),
  obsm_layers = c("${project.obsm_layer}")
)`
    : "";
}

export const useConnect = ({ project }: EmbeddingButtonProps) => {
  const [isOpen, setIsOpen] = useState(false);
  const [isCopied, setIsCopied] = useState(false);
  const [language, setLanguage] = useState<string>("python");
  const [uriTopPosition, setURITopPosition] = useState<number>(-1);
  const [lineHeight, setLineHeight] = useState<number>(-1);

  const projectTier = getProjectTier(project);

  const handleButtonClick = useCallback(() => {
    if (!isOpen)
      track(EVENTS.CENSUS_EMBEDDING_CLICKED, {
        project: project.title,
        category: projectTier,
      });
    setIsOpen(!isOpen);
  }, [isOpen, projectTier, project.title]);

  const uri = `s3://cellxgene-contrib-public/contrib/cell-census/soma/${project.census_version}/${project.id}`;

  const codeSnippet =
    language === "python"
      ? pythonCodeSnippet(project, uri)
      : rCodeSnippet(project);

  const codeSnippetRef = useCallback(
    (node: HTMLDivElement) => {
      if (node !== null) {
        const lines = node.innerText.split("\n");

        const newLineHeight =
          node.clientHeight /
          (lines.length - NUMBER_OF_EXTRA_LINES + NUMBER_OF_PADDING_LINES);

        // index of the URI line
        const lineIndex = lines.findIndex((line: string) => line.includes(uri));

        setURITopPosition(
          newLineHeight * (lineIndex + 1 + NUMBER_OF_PADDING_LINES) +
            LINE_HIGHLIGHT_BACKGROUND_PADDING / 2
        );
        setLineHeight(newLineHeight + LINE_HIGHLIGHT_BACKGROUND_PADDING);
      }
    },
    [uri]
  );

  const [notebookLink, setNotebookLink] = useState("");
  useEffect(() => {
    if (projectTier === "maintained") {
      if (language === "python") {
        setNotebookLink(MAINTAINED_PYTHON_NOTEBOOK_LINK);
      } else {
        setNotebookLink(MAINTAINED_R_NOTEBOOK_LINK);
      }
    } else {
      if (language === "python") {
        setNotebookLink(HOSTED_PYTHON_NOTEBOOK_LINK);
      }
    }
  }, [language, projectTier]);

  const handleCopyMouseEnter = () => setIsCopied(false);

  return {
    isOpen,
    isCopied,
    language,
    codeSnippet,
    projectTier,
    uri,
    uriTopPosition,
    lineHeight,
    notebookLink,
    codeSnippetRef,
    setLanguage,
    handleButtonClick,
    handleCopyMouseEnter,
  };
};
