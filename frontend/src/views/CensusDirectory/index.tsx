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
import {
  fontWeightBold,
  fontWeightMedium,
  fontWeightRegular,
  fontWeightSemibold,
  gray400,
  spacesDefault,
  spacesL,
  spacesXl,
  spacesXxs,
  textSecondary,
} from "src/common/theme";

function DetailItem(props: { label: string; children: string }) {
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
      {props.children}
    </ItemContainer>
  );
}

function CensusDirectory() {
  const Content = styled.div`
    box-sizing: content-box;
    padding: 0px 40px;
    display: flex;
    flex-direction: column;
    margin: 80px auto;
    max-width: 1200px;
  `;

  const Header = styled.h1`
    ${fontHeaderXxl}
    margin-bottom: ${spacesDefault}px;
    font-weight: ${fontWeightBold};
  `;

  const Paragraph = styled.p`
    ${fontBodyS}
    font-weight: ${fontWeightRegular};
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
    margin-bottom: ${spacesDefault}px;
    font-weight: ${fontWeightSemibold};
  `;

  const TierDescription = styled.p`
    ${fontBodyS}
    color: ${textSecondary};
    font-weight: ${fontWeightRegular};
    margin-bottom: 0;
  `;

  const ProjectTitle = styled.h4`
    ${fontHeaderL}
    font-weight: ${fontWeightSemibold};
    margin-bottom: ${spacesDefault}px;
  `;

  const ProjectSubmitter = styled.h4`
    ${fontBodyS}
    font-weight: ${fontWeightSemibold};
    margin-bottom: ${spacesDefault}px;
  `;

  const ProjectDesctiption = styled(Paragraph)`
    max-width: 85ch;
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
    gap: ${spacesDefault}px;
  `;
  const ProjectDetails = styled.div`
    display: flex;
    flex-direction: column;
  `;
  const DetailsContainer = styled.div`
    display: flex;
    flex-direction: row;
    gap: ${spacesXl}px;
    margin-top: ${spacesL}px;
  `;

  const StyledButton = styled(Button)`
    font-weight: ${fontWeightMedium};
    min-width: 80px;
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
            <ProjectDesctiption>
              Ligula imperdiet eget et enim id morbi. Pretium diam risus
              placerat felis vulputate adipiscing sed integer. Mauris commodo
              risus scelerisque tempus mi venenatis egestas. Sed at scelerisque
              vulputate egestas vulputate condimentum libero tempus convallis.
              Nulla id eget fringilla ultrices pellentesque nunc faucibus
              condimentum. Ornare porta eget porttitor cum arcu id ultricies id.
              Massa interdum orci risus arcu mattis massa. Amet metus nibh enim
              nam pellentesque sagittis diam id quam.
            </ProjectDesctiption>
            <DetailsContainer>
              <DetailItem label="contact">Haotian Cui</DetailItem>
              <DetailItem label="publication">
                Cul et al. (2023) bioRxiv
              </DetailItem>
            </DetailsContainer>
          </ProjectDetails>
          <ProjectButtons>
            <StyledButton sdsType="secondary" sdsStyle="square">
              Embedding
            </StyledButton>
            <StyledButton sdsType="primary" sdsStyle="square">
              Model
            </StyledButton>
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
            <ProjectDesctiption>
              Ligula imperdiet eget et enim id morbi. Pretium diam risus
              placerat felis vulputate adipiscing sed integer. Mauris commodo
              risus scelerisque tempus mi venenatis egestas. Sed at scelerisque
              vulputate egestas vulputate condimentum libero tempus convallis.
              Nulla id eget fringilla ultrices pellentesque nunc faucibus
              condimentum. Ornare porta eget porttitor cum arcu id ultricies id.
              Massa interdum orci risus arcu mattis massa. Amet metus nibh enim
              nam pellentesque sagittis diam id quam.
            </ProjectDesctiption>
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
            <StyledButton sdsType="secondary" sdsStyle="square">
              Embedding
            </StyledButton>
            <StyledButton sdsType="primary" sdsStyle="square">
              Model
            </StyledButton>
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
            <ProjectDesctiption>
              Ligula imperdiet eget et enim id morbi. Pretium diam risus
              placerat felis vulputate adipiscing sed integer. Mauris commodo
              risus scelerisque tempus mi venenatis egestas. Sed at scelerisque
              vulputate egestas vulputate condimentum libero tempus convallis.
              Nulla id eget fringilla ultrices pellentesque nunc faucibus
              condimentum. Ornare porta eget porttitor cum arcu id ultricies id.
              Massa interdum orci risus arcu mattis massa. Amet metus nibh enim
              nam pellentesque sagittis diam id quam.
            </ProjectDesctiption>
          </ProjectDetails>
          <ProjectButtons>
            <StyledButton sdsType="secondary" sdsStyle="square">
              Embedding
            </StyledButton>
            <StyledButton sdsType="primary" sdsStyle="square">
              Model
            </StyledButton>
          </ProjectButtons>
        </ProjectContainer>
      </TierContainer>
    </Content>
  );
}

export default CensusDirectory;
