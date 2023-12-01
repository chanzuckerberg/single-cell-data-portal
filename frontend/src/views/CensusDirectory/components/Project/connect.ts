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

  let authorsString = "";
  const primary_affiliation = project.primary_contact?.affiliation || "";
  const affiliations = new Map<string, string[]>();
  affiliations.set(primary_affiliation, [project.primary_contact?.name || ""]);
  project.additional_contacts?.forEach((contact) => {
    const affiliationNames = affiliations.get(contact.affiliation) || [];
    affiliationNames.push(contact.name);
    affiliations.set(contact.affiliation, affiliationNames);
  });

  let index = 0;
  affiliations.forEach((names, affiliation) => {
    if (index > 0) authorsString += " · ";
    index++;

    if (names.length > 1) {
      const last = names.pop();
      authorsString += `${names.join(", ")} & ${last} at ${affiliation}`;
    } else {
      authorsString += `${names[0]} at ${affiliation}`;
    }
  });

  return { date, projectNotebookLinks, projectTier, authorsString };
};
