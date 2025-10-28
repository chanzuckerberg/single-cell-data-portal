import React from "react";
import { Link } from "@czi-sds/components";
import styled from "@emotion/styled";
import { primary400, primary500 } from "src/common/theme";

const StyledLink = styled(Link)`
  color: ${primary400};
  text-decoration: underline dotted;
  text-underline-offset: 2px;

  &:hover {
    color: ${primary500};
    text-decoration: underline solid;
  }
`;

interface FeatureCard {
  icon: "Upload" | "Rocket" | "SearchLinesHorizontal3";
  title: string;
  description: string | React.ReactNode;
  buttonText: string;
  href: string;
}

export const FEATURE_CARDS: FeatureCard[] = [
  {
    icon: "Upload",
    title: "Upload your data",
    description: (
      <span>
        Upload your own single cell data to the Platform to run with CZI-hosted
        AI models. Your uploaded data remains{" "}
        <StyledLink
          target="_blank"
          rel="noopener noreferrer"
          href="https://virtualcellmodels.cziscience.com/docs/ai-workspace/data-policy"
        >
          private and secure
        </StyledLink>
        .
      </span>
    ),
    buttonText: "Learn about Upload",
    href: "https://virtualcellmodels.cziscience.com/docs/ai-workspace/user-guide/data-requirements",
  },
  {
    icon: "Rocket",
    title: "Run Workflows",
    description: `Choose from our pre-built workflow pipelines to generate
          analyses with your own data on CZI-hosted compute. Customize
          steps for your chosen workflows, save your preferences, and
          re-run workflows later with new data.`,
    buttonText: "Learn about Workflows",
    href: "https://virtualcellmodels.cziscience.com/docs/ai-workspace/user-guide/workflows",
  },
  {
    icon: "SearchLinesHorizontal3",
    title: "Explore and Annotate Results",
    description: `View workflow results with the embedding viewer, annotate your findings, or download your results to continue analysis in your own system later.`,
    buttonText: "Learn About Analysis",
    href: "https://virtualcellmodels.cziscience.com/docs/ai-workspace/user-guide/results",
  },
];
