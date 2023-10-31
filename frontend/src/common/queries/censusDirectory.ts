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

// assumed tier 2 for now
export interface CensusDirectoryProject {
  id: string;
  obsoleted_by: string;
  title: string;
  description: string;
  submitter: string;
  contact_name: string;
  contact_email: string;
  doi: string;
  model_link: string;
  census_release: string;
  organism_name: string;
  measurement_name: string;
  embeddings_link: string;
  contribution_type: string;
  notebook_links: string[]; // committed directly
}

export function useProjects() {
  return useQuery<CensusDirectoryProject[]>([USE_PROJECTS], fetchProjects);
}
