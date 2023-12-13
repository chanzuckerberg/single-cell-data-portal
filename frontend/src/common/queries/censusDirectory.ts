import { useQuery } from "react-query";
import { ENTITIES } from "./entities";
import { apiTemplateToUrl } from "../utils/apiTemplateToUrl";
import { CENSUS_SPOTLIGHT_DATA_URL } from "src/configs/configs";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { API } from "../API";

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
  primary_contact: Contact;
  additional_contacts: Contact[];
  DOI?: string;
  publication_info?: string;
  publication_link?: string;
  n_cells?: number;
  n_genes?: number;
  n_columns?: number;
}

export interface Contact {
  name: string;
  email: string;
  affiliation: string;
}

export interface ProjectResponse {
  [accessionID: string]: Project;
}

async function fetchProjects(): Promise<ProjectResponse | undefined> {
  const response = await fetch(
    CENSUS_SPOTLIGHT_DATA_URL + API.CENSUS_SPOTLIGHT_MANIFEST,
    {
      ...DEFAULT_FETCH_OPTIONS,
      ...JSON_BODY_FETCH_OPTIONS,
    }
  );

  try {
    const result = await response.json();
    if (!response.ok) throw result;

    const data = result as ProjectResponse;
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

        const response = await fetch(url);

        const result = await response.json();
        if (!response.ok) throw result;

        const publication_info = parseCrossRefResponse(result);

        data[id].publication_info = publication_info;
        data[id].publication_link = result.URL;
      }
    );
    return data;
  } catch (error) {
    console.log(error);
  }
}

export const USE_PROJECTS = {
  entities: [ENTITIES.CENSUS_DIRECTORY_PROJECTS],
  id: "census-directory-projects",
};

export function useProjects() {
  return useQuery<ProjectResponse | undefined>([USE_PROJECTS], fetchProjects);
}

function parseCrossRefResponse({ message }: any) {
  const author = message.author[0].family;
  const journal = message["short-container-title"];
  const year = message.issued["date-parts"][0][0];

  return `${author} et al. (${year}) ${journal}`;
}
