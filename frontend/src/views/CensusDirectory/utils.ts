import { ProjectProps } from "./components/Project/types";

export const getProjectTier = (project: ProjectProps["project"]) => {
  let projectTier = "hosted";
  if ("tier" in project) projectTier = project.tier;

  return projectTier;
};
