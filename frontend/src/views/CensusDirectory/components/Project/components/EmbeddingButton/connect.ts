import { useCallback, useState } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EmbeddingButtonProps } from "./types";
import { getProjectTier } from "../../../../utils";
import { Project } from "src/common/queries/censusDirectory";
import { StaticProject } from "census-projects.json";

// The div contains two lines of the word copy
const NUMBER_OF_EXTRA_LINES = 2;
// There is .5 em of padding on the top and bottom of the snippets
const NUMBER_OF_PADDING_LINES = 1;
// Total amount of padding around the highlighted line
const LINE_HIGHLIGHT_BACKGROUND_PADDING = 8;

function pythonCodeSnippet(project: StaticProject | Project): string {
  const censusVersion = project.census_version;
  const organism = project.experiment_name;
  const measurement = project.measurement_name;

  const uri = `"s3://cellxgene-contrib-archive/contrib/cell-census/${project.id}"`;

  return project.tier === "maintained"
    ? `
    import cellxgene_census

    # Census version should be gotten from the metadata
    census = cellxgene_census.open_soma(census_version="${project.census_version}")
    adata = cellxgene_census.get_anndata(
        census,
        organism = "${organism}",
        measurement_name = "${measurement}",
        obs_value_filter = "tissue == 'tongue'",
        obsm_layers = "${project.obs_matrix}"
    )`
    : `
  import cellxgene_census
  from cellxgene_census.experimental import get_embedding

  embedding_uri = ${uri}

  # Census version should be gotten from the metadata
  census = cellxgene_census.open_soma(census_version="${censusVersion}")

  adata = cellxgene_census.get_anndata(
      census,
      organism = "${organism}",
      measurement_name = "${measurement}",
      obs_value_filter = "tissue == 'tongue'",
  )

  embeddings = get_embedding("${censusVersion}", embedding_uri, adata.obs["soma_joinid"])

  adata.obsm["emb"] = embeddings`;
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

  const codeSnippet = language === "python" ? pythonCodeSnippet(project) : "";

  // These can be derived from the static S3 namespace + the accessor_id or will be a static url provided in json blob
  const uri = `s3://cellxgene-contrib-archive/contrib/cell-census/${project.id}`;

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
          newLineHeight * lineIndex +
            NUMBER_OF_PADDING_LINES +
            LINE_HIGHLIGHT_BACKGROUND_PADDING / 2
        );
        setLineHeight(newLineHeight + LINE_HIGHLIGHT_BACKGROUND_PADDING);
      }
    },
    [uri]
  );

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
    codeSnippetRef,
    setLanguage,
    handleButtonClick,
    handleCopyMouseEnter,
  };
};
