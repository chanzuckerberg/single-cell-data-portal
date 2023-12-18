import Link from "next/link";
import React from "react";
import {
  type Project as ProjectType,
  useProjects,
} from "src/common/queries/censusDirectory";

import staticProjects from "census-projects.json";
import {
  Content,
  Header,
  DirectoryDescription,
  TierContainer,
  TierTitle,
  TierDescription,
} from "./style";
import Project from "./components/Project";

function CensusDirectory() {
  const { data: projects } = useProjects();

  const hostedProjects = Object.entries(
    projects ?? ({} as ProjectType[])
  ).filter(([_, project]) => !project.revised_by);

  const communityProjects = Object.values(staticProjects).filter(
    (project) => project.tier === "community"
  );
  const maintainedProjects = Object.values(staticProjects).filter(
    (project) => project.tier === "maintained"
  );

  return (
    <Content>
      <Header>Census Models</Header>
      <DirectoryDescription>
        <p>
          This page features models and integrated embeddings of the Census data
          corpus, organized by CELL×GENE’s level of involvement with their
          maintenance and availability. These models are breaking new ground and
          will continue to improve. We encourage you to try out these models and
          provide feedback!
        </p>
        <p>
          {/* TODO: add link to notebooks once available */}
          Please <Link href="">see these tutorials</Link> for usage details.
        </p>
        <p>
          If you’d like to have your project featured here, please{" "}
          <Link href="mailto:cellxgene@chanzuckerberg.com">get in touch</Link>.
        </p>
      </DirectoryDescription>
      {maintainedProjects.length > 0 && (
        <TierContainer>
          <TierTitle>CELL×GENE Maintained Projects</TierTitle>
          <TierDescription>
            These models and their output embeddings are maintained and
            regularly re-trained by CELL×GENE in close collaboration with their
            creators. Embeddings are accessible via the Census API;
            corresponding models are available via CELL×GENE-maintained links.
            <br />
            Please{" "}
            <Link href="mailto:cellxgene@chanzuckerberg.com">
              contact the CELL×GENE team with feedback
            </Link>
            .
          </TierDescription>
          {maintainedProjects.map((project) => (
            <Project key={project.title} project={project} />
          ))}
        </TierContainer>
      )}
      {hostedProjects.length > 0 && (
        <TierContainer>
          <TierTitle>CELL×GENE Hosted Projects</TierTitle>
          <TierDescription>
            CELL×GENE makes these embeddings directly available through the
            Census API, but does not actively maintain or update them.
            Corresponding models are accessible via external links (when
            available).
            <br />
            For issues accessing these embeddings, please{" "}
            <Link href="mailto:cellxgene@chanzuckerberg.com">
              contact the CELL×GENE team
            </Link>
            . For feedback on the embeddings themselves, please contact the
            creators.
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
            Census data. While CELL×GENE does not directly host or maintain
            these projects, we’re excited to showcase them here.
            <br />
            Please contact their creators with questions or feedback.
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
