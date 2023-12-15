import Link from "next/link";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import Toast from "src/views/Collection/components/Toast";
import { UnionProject } from "../../../types";
import { StyledButton } from "../../style";
import { ClobberedProjects } from "src/views/CensusDirectory/utils";

const ModelButton = ({
  project,
  uniqueMetadata,
}: {
  project: UnionProject;
  uniqueMetadata?: ClobberedProjects[number][0];
}) => {
  if (!project.model_link) return null;
  return project.model_link.startsWith("s3") ? (
    <StyledButton
      sdsType="secondary"
      sdsStyle="square"
      onClick={() => {
        track(EVENTS.CENSUS_MODEL_COPIED, {
          project: project.title,
          category: project.tier,
          ...uniqueMetadata,
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
