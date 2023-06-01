import React, { useEffect, useState } from "react";
import { CellCardDescription, Source, Wrapper } from "./style";
import { useDescription, useClDescription } from "src/common/queries/cellCards";
import { Tooltip } from "@czi-sds/components";
import Link from "../common/Link";
import { StyledLink } from "../common/Link/style";

export const CELL_CARD_CL_DESCRIPTION = "cell-card-cl-description";
export const CELL_CARD_GPT_DESCRIPTION = "cell-card-gpt-description";
export const CELL_CARD_GPT_TOOLTIP_LINK = "cell-card-gpt-tooltip-link";

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

  const { data: rawDescriptionGpt } = useDescription(cellTypeId);
  const { data: rawDescriptionCl } = useClDescription(cellTypeId);

  useEffect(() => {
    if (rawDescriptionGpt) setDescriptionGpt(rawDescriptionGpt);
    else setDescriptionGpt("");
    if (rawDescriptionCl) setDescriptionCl(rawDescriptionCl);
    else setDescriptionCl("");
  }, [rawDescriptionGpt, rawDescriptionCl]);

  return (
    <Wrapper>
      {descriptionCl && (
        <CellCardDescription data-testid={CELL_CARD_CL_DESCRIPTION}>
          {descriptionCl}
          <Source>
            {"Source: "}
            <Link
              url={"http://obofoundry.org/ontology/cl.html"}
              label={"Cell Ontology"}
            />
          </Source>
        </CellCardDescription>
      )}
      <br />
      <CellCardDescription data-testid={CELL_CARD_GPT_DESCRIPTION}>
        {descriptionGpt}
        <Source>
          {"Source: "}

          <Tooltip
            leaveDelay={0}
            placement="left"
            width="wide"
            arrow
            title={
              <div>
                {`This summary on \"${cellTypeName}\" was generated with ChatGPT, powered by the GPT3.5 Turbo model. Keep in mind that ChatGPT may occasionally present information that is not entirely accurate. For transparency, the prompts used to generate this summary are shared below.`}
                <br />
                <br />
                <b>System role</b>
                <blockquote>
                  <i>
                    You are a knowledgeable cell biologist that has professional
                    experience writing and curating accurate and informative
                    descriptions of cell types.
                  </i>
                </blockquote>
                <b>User role</b>
                <blockquote>
                  <i>{`I am making a knowledge-base about cell types. Each cell type is a term from the Cell Ontology and will have its own page with a detailed description of that cell type and its function. Please write me a description for "${cellTypeName}". Please return only the description and no other dialogue. The description should include information about the cell type's function. The description should be at least three paragraphs long.`}</i>
                </blockquote>
              </div>
            }
          >
            <StyledLink
              href={"https://platform.openai.com/docs/models/gpt-3-5"}
              target="_blank"
              data-testid={CELL_CARD_GPT_TOOLTIP_LINK}
            >
              ChatGPT
            </StyledLink>
          </Tooltip>
        </Source>
      </CellCardDescription>
    </Wrapper>
  );
}
