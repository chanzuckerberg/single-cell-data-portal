import {
  Button,
  fontBodyS,
  fontHeaderL,
  fontHeaderM,
  fontHeaderXl,
  fontHeaderXxl,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import React from "react";

function CensusDirectory() {
  const Content = styled.div`
    display: flex;
    flex-direction: column;
    margin: 56px 160px;
  `;

  const Header = styled.h1`
    margin-bottom: 8px;
    font-weight: 700;
    ${fontHeaderXxl}
  `;

  const Paragraph = styled.p`
    font-weight: 400;
    ${fontBodyS}
  `;

  const TierTitle = styled.h3`
    margin-bottom: 8px;
    font-weight: 600;
    ${fontHeaderXl}
  `;

  const ProjectTitle = styled.h4`
    font-weight: 600;
    margin-bottom: 8px;
    ${fontHeaderL}
  `;

  const ProjectSubmitter = styled.h4`
    font-weight: 600;
    margin-bottom: 8px;
    ${fontHeaderM}
  `;

  const ProjectContainer = styled.div`
    display: flex;
    flex-direction: row;
    justify-content: space-between;
  `;
  const ProjectButtons = styled.div`
    display: flex;
    flex-direction: row;
    gap: 8px;
  `;
  const ProjectDetails = styled.div``;
  const SubmissionDetails = styled.div``;
  const DataDetails = styled.div``;

  return (
    <Content>
      <Header>Census Directory</Header>
      <Paragraph>
        Habitant sit tristique pharetra at. Quis ipsum morbi pharetra venenatis
        amet purus aliquam nunc. Mi feugiat elementum nec sagittis enim. Turpis
        purus mi tellus leo in vestibulum enim varius. Ut dictum lobortis in
        non. Sed rhoncus enim pharetra pulvinar semper faucibus ut at sapien.
        Parturient pharetra amet sit facilisis sagittis quis. Dignissim
        fermentum consectetur fames vulputate semper neque est non pharetra.
        Amet et elementum neque turpis hac bibendum ac id ipsum.
      </Paragraph>
      <TierTitle>Census Partners</TierTitle>
      <ProjectContainer>
        <ProjectDetails>
          <ProjectTitle> BioAI</ProjectTitle>
          <ProjectSubmitter>
            Turing Institute for Biomedical Machine Learning
          </ProjectSubmitter>
          <Paragraph>
            Ligula imperdiet eget et enim id morbi. Pretium diam risus placerat
            felis vulputate adipiscing sed integer. Mauris commodo risus
            scelerisque tempus mi venenatis egestas. Sed at scelerisque
            vulputate egestas vulputate condimentum libero tempus convallis.
            Nulla id eget fringilla ultrices pellentesque nunc faucibus
            condimentum. Ornare porta eget porttitor cum arcu id ultricies id.
            Massa interdum orci risus arcu mattis massa. Amet metus nibh enim
            nam pellentesque sagittis diam id quam.
          </Paragraph>
          <SubmissionDetails></SubmissionDetails>
          <DataDetails></DataDetails>
        </ProjectDetails>
        <ProjectButtons>
          <Button sdsType="secondary" sdsStyle="square">
            Embedding
          </Button>

          <Button sdsType="primary" sdsStyle="square">
            Model
          </Button>
        </ProjectButtons>
      </ProjectContainer>
      <TierTitle>Tier 2</TierTitle>
      <TierTitle>Tier 3</TierTitle>
    </Content>
  );
}

export default CensusDirectory;
