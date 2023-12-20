/* eslint-disable sonarjs/cognitive-complexity */
import { ProjectProps } from "./types";
import notebookLinks from "census-notebook-links.json";

const DEFAULT_EPOCH_TIME = 0;

export const useConnect = ({ clobberedProjects }: ProjectProps) => {
  const sharedProject = clobberedProjects[0];

  let date = new Date(
    sharedProject.last_updated ||
      sharedProject.submission_date ||
      DEFAULT_EPOCH_TIME
  );

  if (date.getDate() === DEFAULT_EPOCH_TIME) {
    clobberedProjects[1].forEach((project) => {
      const projectDate = new Date(
        project.last_updated || project.submission_date || DEFAULT_EPOCH_TIME
      );
      // check if the project date is more recent than the current date
      if (projectDate > date) {
        date = projectDate;
      }
    });
  }

  const formattedDate = date.toLocaleDateString("en-US", {
    dateStyle: "long",
  });

  const projectNotebookLinks: [string, string][] =
    ("notebook_links" in sharedProject
      ? sharedProject.notebook_links
      : notebookLinks[sharedProject.id ?? ""]) ?? [];

  if (projectNotebookLinks.length === 0) {
    clobberedProjects[1].forEach((project) => {
      projectNotebookLinks?.push(
        ...(("notebook_links" in project
          ? project.notebook_links
          : notebookLinks[project.id ?? ""]) ?? [])
      );
    });
  }

  const projectTier = sharedProject.tier;

  let authorsString = "";
  const primary_affiliation = sharedProject.primary_contact?.affiliation || "";
  const affiliations = new Map<string, string[]>();
  affiliations.set(primary_affiliation, [
    sharedProject.primary_contact?.name || "",
  ]);
  sharedProject.additional_contacts?.forEach((contact) => {
    const affiliationNames = affiliations.get(contact.affiliation) || [];
    affiliationNames.push(contact.name);
    affiliations.set(contact.affiliation, affiliationNames);
  });

  affiliations.forEach((names, affiliation) => {
    authorsString.length > 0 && (authorsString += " · ");

    if (names.length > 1) {
      const last = names.pop();
      authorsString += `${names.join(", ")} & ${last} at ${affiliation}`;
    } else {
      authorsString += `${names[0]} at ${affiliation}`;
    }
  });

  return {
    date: formattedDate,
    projectNotebookLinks,
    projectTier,
    authorsString,
    sharedProject,
  };
};
