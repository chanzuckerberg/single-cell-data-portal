interface FeatureCard {
  icon: "Upload" | "Rocket" | "SearchLinesHorizontal3";
  title: string;
  description: string;
  buttonText: string;
  href: string;
}

export const FEATURE_CARDS: FeatureCard[] = [
  {
    icon: "Upload",
    title: "Upload your own data",
    description:
      "Upload your own single cell data to enable analysis with the Platform's workflows. Data is private to you.",
    buttonText: "Learn about Upload",
    href: "https://virtualcellmodels.cziscience.com/docs/ai-workspace/user-guide/data-requirements",
  },
  {
    icon: "Rocket",
    title: "Run Workflows",
    description:
      "Choose from our pre-built workflow pipelines to generate no-code. Customize specific steps to ensure your results are relevant to your research. Save workflows to re-run later with new data.",
    buttonText: "Learn about Workflows",
    href: "https://virtualcellmodels.cziscience.com/docs/ai-workspace/user-guide/workflows",
  },
  {
    icon: "SearchLinesHorizontal3",
    title: "Explore and Annotate Results",
    description:
      "Explore your workflow results in Explorer embed viewer. Annotate your findings. Download to continue your analysis in your own system.",
    buttonText: "Learn About Analysis",
    href: "https://virtualcellmodels.cziscience.com/docs/ai-workspace/user-guide/results",
  },
];
