import {
  Button,
  fontBodyS,
  fontCapsXxxs,
  fontHeaderL,
  fontHeaderXl,
  fontHeaderXxl,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import React from "react";

function DetailItem(props: { label: string; children: string }) {
  const ItemContainer = styled.div`
    display: flex;
    flex-direction: column;
  `;
  const ItemLabel = styled.div`
    ${fontCapsXxxs}
    font-weight: 600;
    font-feature-settings:
      "clig" off,
      "liga" off;
    color: #959595;
  `;

  return (
    <ItemContainer>
      <ItemLabel>{props.label}</ItemLabel>
      {props.children}
    </ItemContainer>
  );
}

function CensusDirectory() {
  const Content = styled.div`
    display: flex;
    flex-direction: column;
    margin: 56px 160px;
  `;

  const Header = styled.h1`
    ${fontHeaderXxl}
    margin-bottom: 8px;
    font-weight: 700;
  `;

  const Paragraph = styled.p`
    ${fontBodyS}
    font-weight: 400;
    margin-bottom: 0;
  `;

  const DirectoryDescription = styled(Paragraph)`
    margin-bottom: 80px;
  `;

  const TierContainer = styled.div`
    margin-bottom: 120px;
  `;

  const TierTitle = styled.h3`
    ${fontHeaderXl}
    margin-bottom: 8px;
    font-weight: 600;
  `;

  const TierDescription = styled.p`
    ${fontBodyS}
    color: #767676;
    font-weight: 400;
    margin-bottom: 0;
  `;

  const ProjectTitle = styled.h4`
    ${fontHeaderL}
    font-weight: 600;
    margin-bottom: 8px;
  `;

  const ProjectSubmitter = styled.h4`
    ${fontBodyS}
    font-weight: 600;
    margin-bottom: 8px;
  `;

  const ProjectContainer = styled.div`
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    margin-top: 20px;
  `;
  const ProjectButtons = styled.div`
    display: flex;
    flex-direction: row;
    gap: 8px;
  `;
  const ProjectDetails = styled.div`
    display: flex;
    flex-direction: column;
  `;
  const DetailsContainer = styled.div`
    display: flex;
    flex-direction: row;
    gap: 24px;
    margin-top: 16px;
  `;

  return (
    <Content>
      <Header>Census Directory</Header>
      <DirectoryDescription>
        Habitant sit tristique pharetra at. Quis ipsum morbi pharetra venenatis
        amet purus aliquam nunc. Mi feugiat elementum nec sagittis enim. Turpis
        purus mi tellus leo in vestibulum enim varius. Ut dictum lobortis in
        non. Sed rhoncus enim pharetra pulvinar semper faucibus ut at sapien.
        Parturient pharetra amet sit facilisis sagittis quis. Dignissim
        fermentum consectetur fames vulputate semper neque est non pharetra.
        Amet et elementum neque turpis hac bibendum ac id ipsum.
      </DirectoryDescription>
      <TierContainer>
        <TierTitle>Census Partners</TierTitle>
        <TierDescription>
          Ut nisi non lorem adipiscing. Orci tellus quisque quam ac purus vitae.
          Aliquet quis egestas viverra nulla quis lectus adipiscing.
        </TierDescription>
        <ProjectContainer>
          <ProjectDetails>
            <ProjectTitle>BioAI</ProjectTitle>
            <ProjectSubmitter>
              Turing Institute for Biomedical Machine Learning
            </ProjectSubmitter>
            <Paragraph>
              Ligula imperdiet eget et enim id morbi. Pretium diam risus
              placerat felis vulputate adipiscing sed integer. Mauris commodo
              risus scelerisque tempus mi venenatis egestas. Sed at scelerisque
              vulputate egestas vulputate condimentum libero tempus convallis.
              Nulla id eget fringilla ultrices pellentesque nunc faucibus
              condimentum. Ornare porta eget porttitor cum arcu id ultricies id.
              Massa interdum orci risus arcu mattis massa. Amet metus nibh enim
              nam pellentesque sagittis diam id quam.
            </Paragraph>
            <DetailsContainer>
              <DetailItem label="contact">Haotian Cui</DetailItem>
              <DetailItem label="publication">
                Cul et al. (2023) bioRxiv
              </DetailItem>
            </DetailsContainer>
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
      </TierContainer>
      <TierContainer>
        <TierTitle>Tier 2</TierTitle>
        <TierDescription>
          Ut nisi non lorem adipiscing. Orci tellus quisque quam ac purus vitae.
          Aliquet quis egestas viverra nulla quis lectus adipiscing.
        </TierDescription>
        <ProjectContainer>
          <ProjectDetails>
            <ProjectTitle>scGPT</ProjectTitle>
            <ProjectSubmitter>
              WangLab at University of Toronto
            </ProjectSubmitter>
            <Paragraph>
              Ligula imperdiet eget et enim id morbi. Pretium diam risus
              placerat felis vulputate adipiscing sed integer. Mauris commodo
              risus scelerisque tempus mi venenatis egestas. Sed at scelerisque
              vulputate egestas vulputate condimentum libero tempus convallis.
              Nulla id eget fringilla ultrices pellentesque nunc faucibus
              condimentum. Ornare porta eget porttitor cum arcu id ultricies id.
              Massa interdum orci risus arcu mattis massa. Amet metus nibh enim
              nam pellentesque sagittis diam id quam.
            </Paragraph>
            <DetailsContainer>
              <DetailItem label="contact">Haotian Cui</DetailItem>
              <DetailItem label="publication">
                Cul et al. (2023) bioRxiv
              </DetailItem>
            </DetailsContainer>
            <DetailsContainer>
              <DetailItem label="Census Version">LTS 1.0</DetailItem>
              <DetailItem label="experiment">homo_sapiens</DetailItem>
              <DetailItem label="measurement">RNA</DetailItem>
              <DetailItem label="embedding">obs</DetailItem>
            </DetailsContainer>
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
      </TierContainer>
      <TierContainer>
        <TierTitle>Tier 3</TierTitle>
        <TierDescription>
          Ut nisi non lorem adipiscing. Orci tellus quisque quam ac purus vitae.
          Aliquet quis egestas viverra nulla quis lectus adipiscing.
        </TierDescription>
        <ProjectContainer>
          <ProjectDetails>
            <ProjectTitle>OpenCell ML</ProjectTitle>
            <ProjectSubmitter>CZ BioHub </ProjectSubmitter>
            <Paragraph>
              Ligula imperdiet eget et enim id morbi. Pretium diam risus
              placerat felis vulputate adipiscing sed integer. Mauris commodo
              risus scelerisque tempus mi venenatis egestas. Sed at scelerisque
              vulputate egestas vulputate condimentum libero tempus convallis.
              Nulla id eget fringilla ultrices pellentesque nunc faucibus
              condimentum. Ornare porta eget porttitor cum arcu id ultricies id.
              Massa interdum orci risus arcu mattis massa. Amet metus nibh enim
              nam pellentesque sagittis diam id quam.
            </Paragraph>
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
      </TierContainer>
    </Content>
  );
}

export default CensusDirectory;
