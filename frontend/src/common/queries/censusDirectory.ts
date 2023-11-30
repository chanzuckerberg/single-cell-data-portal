import { useQuery } from "react-query";
import { ENTITIES } from "./entities";
import { apiTemplateToUrl } from "../utils/apiTemplateToUrl";

export interface Project {
  census_version: string;
  experiment_name: string;
  measurement_name: string;
  submission_date: string;
  n_features: number;
  data_type: string;
  revised_by?: string;
  model_link?: string;
  last_updated?: string;
  title: string;
  description: string;
  contact_name: string;
  contact_affiliation: string;
  contact_email: string;
  DOI?: string;
  publication_info?: string;
  publication_link?: string;
  n_cells?: number;
  n_genes?: number;
  n_columns?: number;
}

export interface ProjectResponse {
  [accessionID: string]: Project;
}

async function fetchProjects(): Promise<ProjectResponse> {
  // const url = API_URL + API.CENSUS_DIRECTORY;

  // const response = await fetch(url, {
  //   ...DEFAULT_FETCH_OPTIONS,
  //   ...JSON_BODY_FETCH_OPTIONS,
  // });

  // const result = await response.json();
  // if (!response.ok) throw result;

  // return result.data;
  const data: ProjectResponse = {
    accession_ID_1: {
      census_version: "2099-08-11",
      experiment_name: "homo_sapiens",
      measurement_name: "RNA",
      submission_date: "2023-01-11",
      n_features: 2,
      data_type: "obs_embedding",
      revised_by: "accession_ID_2",
      DOI: "10.3352/jeehp.2013.10.3",
      title: "brief title",
      description: "...",
      contact_name: "...",
      contact_email: "...",
      contact_affiliation: "...",
    },
    accession_ID_2: {
      census_version: "2099-10-11",
      submission_date: "2023-01-11",
      last_updated: "2023-11-02",
      n_features: 2,
      data_type: "obs_embedding",
      measurement_name: "RNA",
      experiment_name: "homo_sapiens",
      model_link: "http://example.com",
      contact_affiliation: "Turing Institute for Biomedical Machine Learning",
      contact_name: "Haotian Cu",
      contact_email: "mailto:name@example.com",
      DOI: "10.3352/jeehp.2013.10.3",
      title: "BioAI",
      description:
        "Ligula imperdiet eget et enim id morbi. Pretium diam risus placerat felis vulputate adipiscing sed integer. Mauris commodo risus scelerisque tempus mi venenatis egestas. Sed at scelerisque vulputate egestas vulputate condimentum libero tempus convallis. Nulla id eget fringilla ultrices pellentesque nunc faucibus condimentum. Ornare porta eget porttitor cum arcu id ultricies id. Massa interdum orci risus arcu mattis massa. Amet metus nibh enim nam pellentesque sagittis diam id quam.",
    },
  };

  Object.entries(data).forEach(
    async ([id, project]: [
      keyof ProjectResponse,
      ProjectResponse[keyof ProjectResponse],
    ]) => {
      if (!project.DOI) return;

      // include a mailto: query param to insure reliable service
      const url = apiTemplateToUrl(
        "https://api.crossref.org/works/{DOI}?mailto=cellxgene@cziscience.com",
        {
          DOI: encodeURIComponent(project.DOI),
        }
      );
      try {
        const response = await fetch(url);

        const result = await response.json();
        if (!response.ok) throw result;

        const publication_info = parseCrossRefResponse(result);

        data[id].publication_info = publication_info;
        data[id].publication_link = result.URL;
      } catch (error) {
        console.log(error);
      }
    }
  );

  return data;
}

export const USE_PROJECTS = {
  entities: [ENTITIES.CENSUS_DIRECTORY_PROJECTS],
  id: "census-directory-projects",
};

export function useProjects() {
  return useQuery<ProjectResponse>([USE_PROJECTS], fetchProjects);
}

function parseCrossRefResponse({ message }: any) {
  const author = message.author[0].family;
  const journal = message["short-container-title"];
  const year = message.issued["date-parts"][0][0];

  return `${author} et al. (${year}) ${journal}`;
}
