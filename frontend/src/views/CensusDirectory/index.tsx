import Link from "next/link";
import React from "react";
import { useProjects } from "src/common/queries/censusDirectory";

import staticProjects from "census-projects.json";
import {
  Content,
  Header,
  DirectoryDescription,
  TierContainer,
  TierTitle,
  TierDescription,
} from "./styles";
import Project from "./components/Project";

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

  return (
    <Content>
      <Header>Models & Embeddings Using Census</Header>
      <DirectoryDescription>
        This page features models and integrated embeddings of the Census data
        corpus, organized by CELL×GENE’s level of involvement with their
        maintenance and availability. <br />
        <br /> If you’d like to have your project featured here, please{" "}
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
          {maintainedProjects.map((project) => (
            <Project key={project.title} project={project} />
          ))}
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
          {hostedProjects.map(([id, project]) => (
            <Project key={id} id={id} project={project} />
          ))}
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
          {communityProjects.map((project) => (
            <Project key={project.title} project={project} />
          ))}
        </TierContainer>
      )}
    </Content>
  );
}

export default CensusDirectory;
