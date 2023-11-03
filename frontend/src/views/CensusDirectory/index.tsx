import { fontCapsXxxs } from "@czi-sds/components";
import styled from "@emotion/styled";
import Link from "next/link";
import React from "react";
import {
  Project,
  ProjectResponse,
  useProjects,
} from "src/common/queries/censusDirectory";
import { fontWeightSemibold, gray400, spacesXxs } from "src/common/theme";

import staticProjects, { type StaticProject } from "census-projects.json";
import notebookLinks from "census-notebook-links.json";
import {
  ProjectContainer,
  ProjectDetails,
  ProjectTitle,
  ProjectSubmitter,
  ProjectDesctiption,
  DetailsContainer,
  ProjectButtons,
  StyledButton,
  Content,
  Header,
  DirectoryDescription,
  TierContainer,
  TierTitle,
  TierDescription,
} from "./styles";
import EmbeddingsButton from "./components/EmbeddingButton";

function DetailItem(props: { label: string; children: string; link?: string }) {
  const ItemContainer = styled.div`
    display: flex;
    flex-direction: column;
    gap: ${spacesXxs}px;
  `;

  const ItemLabel = styled.div`
    ${fontCapsXxxs}
    font-weight: ${fontWeightSemibold};
    font-feature-settings:
      "clig" off,
      "liga" off;
    color: ${gray400};
  `;

  return (
    <ItemContainer>
      <ItemLabel>{props.label}</ItemLabel>
      {props.link ? (
        <Link href={props.link} passHref>
          {props.children}
        </Link>
      ) : (
        props.children
      )}
    </ItemContainer>
  );
}

function CensusDirectory() {
  const { data: projects } = useProjects();

  const hostedProjects = Object.entries(projects || {}).filter(
    ([_, project]) => !project.revised_by
  );

  const communityProjects = Object.values(staticProjects).filter(
    (project) => project.tier === 1
  );
  const maintainedProjects = Object.values(staticProjects).filter(
    (project) => project.tier === 3
  );

  const renderDetailItem = (label: string, value?: string, link?: string) => {
    return value ? (
      <DetailItem label={label} link={link}>
        {value}
      </DetailItem>
    ) : null;
  };

  const renderProject = (
    project: StaticProject | Project,
    id?: keyof ProjectResponse
  ) => (
    <ProjectContainer key={project.title}>
      <ProjectDetails>
        <ProjectTitle>{project.title}</ProjectTitle>
        <ProjectSubmitter>{project.contact_affiliation}</ProjectSubmitter>
        <ProjectDesctiption>{project.description}</ProjectDesctiption>
        <DetailsContainer>
          {renderDetailItem(
            "contact",
            project.contact_name,
            project.contact_email
          )}
          {renderDetailItem(
            "publication",
            project.publication_info,
            project.publication_link
          )}
          {renderDetailItem(
            "Last Updated",
            //convert date to month, day, year
            new Date(
              project.last_updated || project.submission_date || ""
            ).toLocaleDateString("en-US", {
              dateStyle: "long",
            })
          )}
        </DetailsContainer>
        <DetailsContainer>
          {renderDetailItem("Census Version", project.census_version)}
          {renderDetailItem("experiment", project.experiment_name)}
          {renderDetailItem("measurement", project.measurement_name)}
          {renderDetailItem("embedding", project.data_type)}
          {"notebook_links" in project
            ? project.notebook_links?.map((link) =>
                renderDetailItem("notebook", link[0], link[1])
              )
            : id &&
              notebookLinks[id]?.map((link) =>
                renderDetailItem("notebook", link[0], link[1])
              )}
        </DetailsContainer>
      </ProjectDetails>
      <ProjectButtons>
        <EmbeddingsButton />
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
  return (
    <Content>
      <Header>Census Directory</Header>
      <DirectoryDescription>
        This page features models and integrated embeddings of the Census data
        corpus, organized by CELL×GENE’s level of involvement with their
        maintenance and availability. If you’d like to have your project
        featured here, please
        <Link href="mailto:cellxgene@chanzuckerberg.com">get in touch</Link>!
      </DirectoryDescription>
      {maintainedProjects.length > 0 && (
        <TierContainer>
          <TierTitle>CELL×GENE Maintained Projects</TierTitle>
          <TierDescription>
            These projects are actively maintained and regularly updated by
            CELLxGENE in close collaboration with their creators. Embeddings are
            accessible via the Census API; models are available via
            CELLxGENE-maintained links.
          </TierDescription>
          {maintainedProjects.map((project) => renderProject(project))}
        </TierContainer>
      )}
      {hostedProjects.length > 0 && (
        <TierContainer>
          <TierTitle>CELLxGENE Hosted Projects</TierTitle>
          <TierDescription>
            CELLxGENE makes these projects available, but does not actively
            maintain or update them. Embeddings are accessible via the Census
            API; models are available via CELLxGENE-maintained links.
          </TierDescription>
          {hostedProjects.map(([id, project]) => renderProject(project, id))}
        </TierContainer>
      )}
      {communityProjects.length > 0 && (
        <TierContainer>
          <TierTitle>Community Projects</TierTitle>
          <TierDescription>
            The community has also developed many wonderful projects using
            Census data. While CELLxGENE does not directly host or maintain
            these projects, we’re excited to showcase them here.
          </TierDescription>
          {communityProjects.map((project) => renderProject(project))}
        </TierContainer>
      )}
    </Content>
  );
}

export default CensusDirectory;
