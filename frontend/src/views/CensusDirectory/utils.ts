import { ProjectProps } from "./components/Project/types";

export const getProjectTier = (project: ProjectProps["project"]) => {
  let projectTier = "hosted";
  if ("tier" in project)
    switch (project.tier ?? 2) {
      case 1:
        projectTier = "maintained";
        break;
      case 2:
        projectTier = "hosted";
        break;
      case 3:
        projectTier = "community";
        break;
    }

  return projectTier;
};
