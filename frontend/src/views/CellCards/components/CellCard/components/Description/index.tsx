import React, { useEffect, useState } from "react";
import { useRouter } from "next/router";
import { CellCardDescription, Wrapper } from "./style";
import { useDescription, useClDescription } from "src/common/queries/cellCards";
import { Tooltip } from "czifui";

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
          <b>
            <i>
              The below description is sourced from the{" "}
              <a
                href={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}`}
                target="_blank"
              >
                Cell Ontology
              </a>
              {".\n\n"}
            </i>
          </b>
          {descriptionCl}
        </CellCardDescription>
      )}
      <br />
      <CellCardDescription>
        <b>
          <i>
            {
              "The below summary is produced by ChatGPT. ChatGPT may produce inaccurate information about people, places, or facts. "
            }
          </i>
        </b>
        {"\n"}
        <Tooltip
          placement="right"
          width="wide"
          arrow
          title={
            <div>
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
          <a>View prompt</a>
        </Tooltip>
        {"\n\n"}
        {descriptionGpt}
      </CellCardDescription>
    </Wrapper>
  );
}
