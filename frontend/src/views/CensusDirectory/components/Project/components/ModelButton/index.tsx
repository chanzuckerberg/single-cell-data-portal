import { StaticProject } from "census-projects.json";
import Link from "next/link";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { Project } from "src/common/queries/censusDirectory";
import { StyledButton } from "src/views/CensusDirectory/style";
import Toast from "src/views/Collection/components/Toast";

const ModelButton = ({ project }: { project: StaticProject | Project }) => {
  if (!project.model_link) return null;
  return project.model_link.startsWith("s3") ? (
    <StyledButton
      sdsType="secondary"
      sdsStyle="square"
      onClick={() => {
        track(EVENTS.CENSUS_MODEL_CLICKED, {
          project: project.title,
          category: project.tier,
        });
        // copy URI to clipboard
        navigator.clipboard.writeText(project.model_link || "");
        // show toast to verify copy
        Toast.show({
          message: "Model URI copied to clipboard",
          intent: "success",
        });
      }}
    >
      Copy Model URI
    </StyledButton>
  ) : (
    <Link href={project.model_link}>
      <StyledButton
        sdsType="secondary"
        sdsStyle="square"
        onClick={() => {
          track(EVENTS.CENSUS_MODEL_CLICKED, {
            project: project.title,
            category: project.tier,
          });
        }}
      >
        Model Page
      </StyledButton>
    </Link>
  );
};

export default ModelButton;
