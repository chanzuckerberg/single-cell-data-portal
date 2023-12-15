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
import { clobberAndDifferentiateProjectMetadata } from "./utils";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";

function CensusDirectory() {
  const { data: projects } = useProjects();

  const hostedProjects = clobberAndDifferentiateProjectMetadata(
    Object.values(projects ?? ({} as ProjectType[])).filter(
      (project) => !project.revised_by
    )
  );

  const communityProjects = clobberAndDifferentiateProjectMetadata(
    Object.values(staticProjects).filter(
      (project) => project.tier === "community"
    )
  );
  const maintainedProjects = clobberAndDifferentiateProjectMetadata(
    Object.values(staticProjects).filter(
      (project) => project.tier === "maintained"
    )
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
          Please{" "}
          <a
            href="https://chanzuckerberg.github.io/cellxgene-census/examples.html"
            target="__blank"
            rel="noopener noreferrer"
            onClick={() => {
              track(EVENTS.CENSUS_MODELS_TUTORIALS_CLICKED);
            }}
          >
            see these tutorials
          </a>{" "}
          for usage details.
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
          {maintainedProjects.map((clobberedProjects) => (
            <Project
              key={clobberedProjects[0].id}
              clobberedProjects={clobberedProjects}
            />
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
          {hostedProjects.map((clobberedProjects) => (
            <Project
              key={clobberedProjects[0].id}
              clobberedProjects={clobberedProjects}
            />
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
          {communityProjects.map((clobberedProjects) => (
            <Project
              key={clobberedProjects[0].id}
              clobberedProjects={clobberedProjects}
            />
          ))}
        </TierContainer>
      )}
    </Content>
  );
}

export default CensusDirectory;
