import React, { useEffect, useState } from "react";
import {
  CellGuideCardDescription,
  ChatGptTooltipSubtext,
  ChatGptTooltipText,
  DescriptionHeader,
  Source,
  SourceLink,
  StyledTooltip,
  Wrapper,
} from "./style";
import { useDescription, useClDescription } from "src/common/queries/cellGuide";
import Link from "../common/Link";
import { StyledLink } from "../common/Link/style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { CELL_GUIDE_CORRECTION_SURVEY_LINK } from "src/common/constants/airtableLinks";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { StyledIconImage } from "../common/HelpTooltip/style";
import {
  CELL_GUIDE_CARD_CL_DESCRIPTION,
  CELL_GUIDE_CARD_GPT_DESCRIPTION,
  CELL_GUIDE_CARD_GPT_TOOLTIP_LINK,
} from "src/views/CellGuide/components/CellGuideCard/components/Description/constants";

interface DescriptionProps {
  cellTypeName: string;
  cellTypeId: string;
}
export default function Description({
  cellTypeId,
  cellTypeName,
}: DescriptionProps): JSX.Element {
  const [descriptionGpt, setDescriptionGpt] = useState<string>("");
  const [descriptionCl, setDescriptionCl] = useState<string>("");

  const [timerId, setTimerId] = useState<NodeJS.Timer | null>(null); // For chatgpt hover event

  const { data: rawDescriptionGpt } = useDescription(cellTypeId);
  const { data: rawDescriptionCl } = useClDescription(cellTypeId);

  useEffect(() => {
    if (rawDescriptionGpt) setDescriptionGpt(rawDescriptionGpt);
    else setDescriptionGpt("");
    if (rawDescriptionCl) setDescriptionCl(rawDescriptionCl);
    else setDescriptionCl("");
  }, [rawDescriptionGpt, rawDescriptionCl]);

  const copyHandler = () => {
    if (!window) return;

    const selectedText = window.getSelection()?.toString().trim();

    if (selectedText !== "") {
      track(EVENTS.CG_COPY_CELL_TYPE_DESCRIPTION);
    }
  };

  return (
    <Wrapper>
      {descriptionCl && (
        <CellGuideCardDescription
          data-testid={CELL_GUIDE_CARD_CL_DESCRIPTION}
          onCopy={copyHandler}
        >
          {descriptionCl}
          <Source>
            <SourceLink>
              {"Source: "}
              <Link
                url={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeId.replace(
                  ":",
                  "_"
                )}`}
                label={"Cell Ontology"}
              />
            </SourceLink>
          </Source>
        </CellGuideCardDescription>
      )}
      <br />
      <CellGuideCardDescription
        data-testid={CELL_GUIDE_CARD_GPT_DESCRIPTION}
        onCopy={copyHandler}
      >
        <DescriptionHeader>Experimental Description</DescriptionHeader>
        {descriptionGpt}
        <Source>
          <div>
            <em>
              We&apos;re still validating ChatGPT descriptions with our
              Biocurator team. Once a description is validated, we&apos;ll add
              references and a validation icon. If you believe a description is
              inaccurate, please{" "}
              <a
                href={CELL_GUIDE_CORRECTION_SURVEY_LINK}
                target="_blank"
                rel="noreferrer noopener"
                onClick={() => {
                  track(EVENTS.SUBMIT_CORRECTION_CLICKED, {
                    cell_type_name: cellTypeName,
                  });
                }}
              >
                submit a correction
              </a>
              .
            </em>
          </div>

          <SourceLink>
            {"Source: ChatGPT "}
            <StyledTooltip
              sdsStyle="dark"
              leaveDelay={0}
              placement="left"
              width="wide"
              arrow
              slotProps={{
                tooltip: {
                  style: {
                    maxWidth: 550, // This is needed because SDS bug where width prop doesn't affect dark sdsStyle
                    textAlign: "start", // Also dark sdsStyle aligns content to center by default
                  },
                },
              }}
              title={
                <div>
                  <ChatGptTooltipText>
                    {`This summary on "${cellTypeName}" was generated with ChatGPT, powered by the GPT4 model. Keep in mind that ChatGPT may occasionally present information that is not entirely accurate. For transparency, the prompts used to generate this summary are shared below. CZI is currently offering this as a pilot feature and we may update or change this feature at our discretion.`}
                  </ChatGptTooltipText>
                  <br />
                  <ChatGptTooltipSubtext>
                    System role: You are a knowledgeable cell biologist that has
                    professional experience writing and curating accurate and
                    informative descriptions of cell types.
                    <br />
                    <br />
                    {`User role: I am making a knowledge-base about cell types. Each cell type is a term from the Cell Ontology and will have its own page with a detailed description of that cell type and its function. Please write me a description for "${cellTypeName}". Please return only the description and no other dialogue. The description should include information about the cell type's function. The description should be at least three paragraphs long.`}
                  </ChatGptTooltipSubtext>
                </div>
              }
            >
              <StyledLink
                href={"https://platform.openai.com/docs/models/gpt-4"}
                target="_blank"
                data-testid={CELL_GUIDE_CARD_GPT_TOOLTIP_LINK}
                onMouseOver={() => {
                  const id = setTimeout(() => {
                    track(EVENTS.CG_CHAT_GPT_HOVER);
                  }, 2000);
                  setTimerId(id);
                }}
                onMouseOut={() => {
                  if (timerId) {
                    clearTimeout(timerId);
                    setTimerId(null);
                  }
                }}
              >
                <StyledIconImage src={questionMarkIcon} />
              </StyledLink>
            </StyledTooltip>
          </SourceLink>
        </Source>
      </CellGuideCardDescription>
    </Wrapper>
  );
}
