import { useQuery } from "react-query";
import { DEFAULT_FETCH_OPTIONS, JSON_BODY_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";

async function fetchProjects() {
  // TMP LOCALDEV ONLY BE SERVER
  // TMP LOCALDEV ONLY BE SERVER
  // TMP LOCALDEV ONLY BE SERVER
  // TMP LOCALDEV ONLY BE SERVER
  console.log("fetching projects");

  const response = await fetch("http://localhost:5005/census-directory", {
    ...DEFAULT_FETCH_OPTIONS,
    ...JSON_BODY_FETCH_OPTIONS,
  });

  const result = await response.json();
  if (!response.ok) throw result;

  return result.data;
}

export const USE_PROJECTS = {
  entities: [ENTITIES.CENSUS_DIRECTORY_PROJECTS],
  id: "census-directory-projects",
};

// maybe we should have a general census directory project type with specific typing for each tier?
export interface CensusDirectoryProject {
  title: string;
  description: string;
  submitter: string;
  contact_name: string;
  contact_email: string;
  doi: string;
  publication_info: string;
  publication_link: string;
  model_link: string;
  embeddings_link: string;
  tier: number;
  census_release: string;
  experiment_name: string;
  measurement_name: string;
  contribution_type: string;
  notebook_links: string[];
}

export function useProjects() {
  return useQuery<CensusDirectoryProject[]>([USE_PROJECTS], fetchProjects);
}
