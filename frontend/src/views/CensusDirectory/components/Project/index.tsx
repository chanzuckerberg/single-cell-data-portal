import Link from "next/link";
import { type Project } from "src/common/queries/censusDirectory";
import {
  ProjectContainer,
  ProjectDetails,
  ProjectTitle,
  ProjectSubmitter,
  ProjectDescription,
  DetailsContainer,
  ProjectButtons,
  StyledButton,
} from "../../styles";
import DetailItem from "../DetailItem";

import EmbeddingButton from "../EmbeddingButton";

import { ProjectProps } from "./types";
import { useConnect } from "./connect";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

const Project = ({ project, id }: ProjectProps) => {
  const { date, projectNotebookLinks, projectTier } = useConnect({
    project,
    id,
  });

  return (
    <ProjectContainer key={project.title}>
      <ProjectDetails>
        <ProjectTitle>{project.title}</ProjectTitle>
        <ProjectSubmitter>{project.contact_affiliation}</ProjectSubmitter>
        <ProjectDescription>{project.description}</ProjectDescription>
        <DetailsContainer>
          <DetailItem
            label="contact"
            link={project.contact_email}
            onClick={() => {
              track(EVENTS.CENSUS_CONTACT_CLICKED, {
                project: project.title,
                contact: project.contact_name,
              });
            }}
          >
            {project.contact_name}
          </DetailItem>
          <DetailItem
            label="publication"
            link={project.publication_link}
            onClick={() => {
              track(EVENTS.CENSUS_PUBLICATION_CLICKED, {
                publication: project.publication_info,
                project: project.title,
              });
            }}
          >
            {project.publication_info}
          </DetailItem>
          <DetailItem label="Last Updated">{date}</DetailItem>
        </DetailsContainer>
        <DetailsContainer>
          <DetailItem label="Census Version">
            {project.census_version}
          </DetailItem>
          <DetailItem label="experiment">{project.experiment_name}</DetailItem>
          <DetailItem label="measurement">
            {project.measurement_name}
          </DetailItem>
          <DetailItem label="embedding">{project.data_type}</DetailItem>
          {projectNotebookLinks?.map((link) => (
            <DetailItem
              label="notebook"
              link={link[1]}
              key={link[1]}
              onClick={() => {
                track(EVENTS.CENSUS_NOTEBOOK_CLICKED, {
                  project: project.title,
                  category: projectTier,
                  notebook: link[0],
                });
              }}
            >
              {link[0]}
            </DetailItem>
          ))}
          <DetailItem label="cells">{project.n_cells}</DetailItem>
          <DetailItem label="genes">{project.n_genes}</DetailItem>
          <DetailItem label="columns">{project.n_columns}</DetailItem>
        </DetailsContainer>
      </ProjectDetails>
      <ProjectButtons>
        {"project_page" in project && !!project.project_page && (
          <Link href={project.project_page}>
            <StyledButton
              sdsType="secondary"
              sdsStyle="square"
              onClick={() => {
                track(EVENTS.CENSUS_PROJECT_LINK_CLICKED, {
                  project: project.title,
                  category: projectTier,
                });
              }}
            >
              Project Page
            </StyledButton>
          </Link>
        )}
        {projectTier === "hosted" && <EmbeddingButton project={project} />}
        {!!project.model_link && (
          <Link href={project.model_link}>
            <StyledButton
              sdsType="primary"
              sdsStyle="square"
              onClick={() => {
                track(EVENTS.CENSUS_MODEL_CLICKED, {
                  project: project.title,
                  category: projectTier,
                });
              }}
            >
              Model
            </StyledButton>
          </Link>
        )}
      </ProjectButtons>
    </ProjectContainer>
  );
};
export default Project;
