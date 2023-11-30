import { getProjectTier } from "../../utils";
import { ProjectProps } from "./types";
import notebookLinks from "census-notebook-links.json";

export const useConnect = ({ project, id }: ProjectProps) => {
  const date = new Date(
    project.last_updated || project.submission_date || ""
  ).toLocaleDateString("en-US", {
    dateStyle: "long",
  });

  const projectNotebookLinks: [string, string][] | undefined =
    "notebook_links" in project
      ? project.notebook_links
      : notebookLinks[id ?? ""];

  const projectTier = getProjectTier(project);

  return { date, projectNotebookLinks, projectTier };
};
