import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ClobberedProjects } from "src/views/CensusDirectory/utils";
import ModelButton from "./components/ModelButton";
import EmbeddingButton from "./components/EmbeddingButton";
import { ButtonsColumn, ButtonsRow, StyledButton } from "./style";
import DetailItem from "../../DetailItem";

const IGNORE_DIFFERENT_METADATA_KEYS = ["model_link", "id"];
const ATTRIBUTE_TO_LABEL: Record<string, string> = {
  experiment_name: "Organism",
  n_cells: "Cells",
};

const ProjectButtons = ({
  clobberedProjects,
}: {
  clobberedProjects: ClobberedProjects[number];
}) => {
  const sharedProject = clobberedProjects[0];
  if (clobberedProjects[1].length === 1) {
    return (
      <ButtonsRow>
        {"project_page" in sharedProject && !!sharedProject.project_page && (
          <a
            href={sharedProject.model_link}
            target="_blank"
            rel="noopener noreferrer"
          >
            <StyledButton
              sdsType="secondary"
              sdsStyle="square"
              onClick={() => {
                track(EVENTS.CENSUS_PROJECT_LINK_CLICKED, {
                  project: sharedProject.title,
                  category: sharedProject.tier,
                });
              }}
            >
              Project Page
            </StyledButton>
          </a>
        )}
        <ModelButton project={clobberedProjects[1][0]} />
        <EmbeddingButton project={clobberedProjects[1][0]} />
      </ButtonsRow>
    );
  }

  return (
    <ButtonsColumn>
      {clobberedProjects[1].map((project) => {
        const uniqueMetadata = Object.fromEntries(
          Object.entries(project).filter(([key, value]) => {
            return !(key in sharedProject) && value;
          })
        );

        return (
          <ButtonsRow key={project.id ?? project.experiment_name}>
            {Object.entries(uniqueMetadata)
              .filter(([key]) => !IGNORE_DIFFERENT_METADATA_KEYS.includes(key))
              .map(([key, value]) => {
                return (
                  <DetailItem key={key} label={ATTRIBUTE_TO_LABEL[key]}>
                    {value}
                  </DetailItem>
                );
              })}
            <ModelButton project={project} uniqueMetadata={uniqueMetadata} />
            <EmbeddingButton
              project={project}
              uniqueMetadata={uniqueMetadata}
            />
          </ButtonsRow>
        );
      })}
    </ButtonsColumn>
  );
};

export default ProjectButtons;
