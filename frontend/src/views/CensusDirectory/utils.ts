import { UnionProject } from "./components/Project/types";

export const getProjectTier = (project: UnionProject) => {
  let projectTier = "hosted";
  if ("tier" in project) projectTier = project.tier;

  return projectTier;
};

export type ClobberedProjects = [Partial<UnionProject>, UnionProject[]][];

function isEqual(obj1: any, obj2: any): boolean {
  Object.keys(obj1).forEach((key) => {
    if (typeof obj1[key] === "object" && typeof obj2[key] === "object") {
      if (!isEqual(obj1[key], obj2[key])) {
        return false;
      }
    } else if (obj1[key] !== obj2[key]) {
      return false;
    }
  });
  return true;
}

// Iterates over all projects and returns a list of projects that have the same title
// and the overlapping data itself
export const clobberAndDifferentiateProjectMetadata = (
  projects: UnionProject[]
): ClobberedProjects => {
  const clobberedProjects: ClobberedProjects = [];
  // get unique project titles and the projects that have that title
  const titlesToProjects = projects.reduce(
    (acc, project) => {
      if (!project.title) return acc;
      if (!acc[project.title]) acc[project.title] = [];
      acc[project.title].push(project);
      return acc;
    },
    {} as Record<string, UnionProject[]>
  );

  Object.values(titlesToProjects).forEach((projects) => {
    if (projects.length > 1) {
      const clobberedProject = projects.reduce(
        (acc, project) => {
          (
            Object.keys(project) as unknown as Array<keyof UnionProject>
          ).forEach((key) => {
            if (
              acc[key] &&
              typeof acc[key] === "object" &&
              typeof project[key] === "object"
            ) {
              if (!isEqual(acc[key], project[key])) {
                delete acc[key];
              }
            } else if (acc[key] && acc[key] !== project[key]) {
              delete acc[key];
            }
          });
          return acc;
        },
        { ...projects[0] }
      );

      clobberedProjects.push([clobberedProject, projects]);
    } else {
      clobberedProjects.push([projects[0], projects]);
    }
  });

  return clobberedProjects;
};
