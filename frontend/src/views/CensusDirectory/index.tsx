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
          <Link href="mailto:soma@chanzuckerberg.com">get in touch</Link>.
        </p>
      </DirectoryDescription>
      {maintainedProjects.length > 0 && (
        <TierContainer>
          <TierTitle>CELL×GENE Collaboration Projects</TierTitle>
          <TierDescription>
            These models and their output embeddings are ongoing collaborations.
            CZI and the partner labs are improving the models as the Census
            resource grows. Embeddings are accessible via the Census API;
            corresponding models are available for download.
            <br />
            Please{" "}
            <Link href="mailto:soma@chanzuckerberg.com">
              contact the CELL×GENE team with feedback
            </Link>
            .
          </TierDescription>
          {maintainedProjects.map((clobberedProjects, key) => {
            return (
              <Project
                key={clobberedProjects[0].id || key}
                clobberedProjects={clobberedProjects}
              />
            );
          })}
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
            <Link href="mailto:soma@chanzuckerberg.com">
              contact the CELL×GENE team
            </Link>
            . For feedback on the embeddings themselves, please contact the
            creators.
          </TierDescription>
          {hostedProjects.map((clobberedProjects, key) => (
            <Project
              key={clobberedProjects[0].id || key}
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
          {communityProjects.map((clobberedProjects, key) => (
            <Project
              key={clobberedProjects[0].id || key}
              clobberedProjects={clobberedProjects}
            />
          ))}
        </TierContainer>
      )}
    </Content>
  );
}

export default CensusDirectory;
