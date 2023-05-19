import React, { useEffect, useState } from "react";
import { useRouter } from "next/router";
import { CellCardDescription, Source, Wrapper } from "./style";
import { useDescription, useClDescription } from "src/common/queries/cellCards";
import { Tooltip } from "czifui";
import Link from "../common/Link";
import { StyledLink } from "../common/Link/style";

interface DescriptionProps {
  cellTypeName: string;
}
export default function Description({
  cellTypeName,
}: DescriptionProps): JSX.Element {
  const [descriptionGpt, setDescriptionGpt] = useState<string>("");
  const [descriptionCl, setDescriptionCl] = useState<string>("");

  const router = useRouter();
  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";

  const { data: rawDescriptionGpt } = useDescription(cellTypeId);
  const { data: rawDescriptionCl } = useClDescription(cellTypeId);

  useEffect(() => {
    if (rawDescriptionGpt) setDescriptionGpt(rawDescriptionGpt);
    else setDescriptionGpt("");
    if (rawDescriptionCl) setDescriptionCl(rawDescriptionCl);
    else setDescriptionCl("");
  }, [cellTypeIdRaw, rawDescriptionGpt, rawDescriptionCl]);

  return (
    <Wrapper>
      {descriptionCl && (
        <CellCardDescription>
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
      <CellCardDescription>
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
            <StyledLink>ChatGPT</StyledLink>
          </Tooltip>
        </Source>
      </CellCardDescription>
    </Wrapper>
  );
}
