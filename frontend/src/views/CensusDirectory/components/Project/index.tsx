import { type Project } from "src/common/queries/censusDirectory";
import {
  ProjectContainer,
  ProjectDetails,
  ProjectTitle,
  ProjectSubmitter,
  ProjectDescription,
  DetailsContainer,
} from "../../style";
import DetailItem from "../DetailItem";

import { ProjectProps } from "./types";
import { useConnect } from "./connect";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import ProjectButtons from "./ProjectButtons";

const Project = ({ clobberedProjects }: ProjectProps) => {
  const {
    date,
    projectNotebookLinks,
    projectTier,
    authorsString,
    sharedProject,
  } = useConnect({
    clobberedProjects,
  });

  return (
    <ProjectContainer key={sharedProject.title}>
      <ProjectDetails>
        <ProjectTitle>{sharedProject.title}</ProjectTitle>
        <ProjectSubmitter>{authorsString}</ProjectSubmitter>
        <ProjectDescription>{sharedProject.description}</ProjectDescription>
        <DetailsContainer>
          <DetailItem
            label="contact"
            link={sharedProject.primary_contact?.email}
            onClick={() => {
              track(EVENTS.CENSUS_CONTACT_CLICKED, {
                project: sharedProject.title,
                contact: sharedProject.primary_contact?.name,
              });
            }}
          >
            {sharedProject.primary_contact?.name}
          </DetailItem>
          <DetailItem
            label="publication"
            link={sharedProject.publication_link}
            onClick={() => {
              track(EVENTS.CENSUS_PUBLICATION_CLICKED, {
                publication: sharedProject.publication_info,
                project: sharedProject.title,
              });
            }}
          >
            {sharedProject.publication_info}
          </DetailItem>
          <DetailItem label="Last Updated">{date}</DetailItem>
        </DetailsContainer>
        <DetailsContainer>
          <DetailItem label="Census Version">
            {sharedProject.census_version}
          </DetailItem>
          <DetailItem label="organism">
            {sharedProject.experiment_name}
          </DetailItem>
          <DetailItem label="measurement">
            {sharedProject.measurement_name}
          </DetailItem>
          <DetailItem label="embedding">{sharedProject.data_type}</DetailItem>
          {projectNotebookLinks?.map((link) => (
            <DetailItem
              label="notebook"
              link={link[1]}
              key={link[1]}
              onClick={() => {
                track(EVENTS.CENSUS_NOTEBOOK_CLICKED, {
                  project: sharedProject.title,
                  category: projectTier,
                  notebook: link[0],
                });
              }}
            >
              {link[0]}
            </DetailItem>
          ))}
          <DetailItem label="cells">{sharedProject.n_cells}</DetailItem>
          <DetailItem label="genes">{sharedProject.n_genes}</DetailItem>
          <DetailItem label="columns">{sharedProject.n_columns}</DetailItem>
        </DetailsContainer>
      </ProjectDetails>
      <ProjectButtons clobberedProjects={clobberedProjects} />
    </ProjectContainer>
  );
};
export default Project;
