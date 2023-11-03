import { StaticProject } from "census-projects.json";
import Link from "next/link";
import { Project, ProjectResponse } from "src/common/queries/censusDirectory";
import {
  ProjectContainer,
  ProjectDetails,
  ProjectTitle,
  ProjectSubmitter,
  ProjectDesctiption,
  DetailsContainer,
  ProjectButtons,
  StyledButton,
} from "../styles";
import DetailItem from "./DetailItem";
interface Props {
  project: StaticProject | Project;
  id?: keyof ProjectResponse;
}

import EmbeddingButton from "./EmbeddingButton";

import notebookLinks from "census-notebook-links.json";

const Project = ({ project, id }: Props) => {
  const date = new Date(
    project.last_updated || project.submission_date || ""
  ).toLocaleDateString("en-US", {
    dateStyle: "long",
  });

  const projectNotebookLinks: [string, string][] | undefined =
    "notebook_links" in project
      ? project.notebook_links
      : notebookLinks[id ?? ""];

  return (
    <ProjectContainer key={project.title}>
      <ProjectDetails>
        <ProjectTitle>{project.title}</ProjectTitle>
        <ProjectSubmitter>{project.contact_affiliation}</ProjectSubmitter>
        <ProjectDesctiption>{project.description}</ProjectDesctiption>
        <DetailsContainer>
          <DetailItem label="contact" link={project.contact_email}>
            {project.contact_name}
          </DetailItem>
          <DetailItem label="publication" link={project.publication_link}>
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
            <DetailItem label="notebook" link={link[1]} key={link[1]}>
              {link[0]}
            </DetailItem>
          ))}
        </DetailsContainer>
      </ProjectDetails>
      <ProjectButtons>
        <EmbeddingButton />
        {project.model_link && (
          <Link href={project.model_link}>
            <StyledButton sdsType="primary" sdsStyle="square">
              Model
            </StyledButton>
          </Link>
        )}
      </ProjectButtons>
    </ProjectContainer>
  );
};
export default Project;
