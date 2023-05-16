import React, { useEffect, useState } from "react";
import { useRouter } from "next/router";
import { CellCardDescription } from "./style";
import { useDescription } from "src/common/queries/cellCards";

export default function Description(): JSX.Element {
  const [descriptionGpt, setDescriptionGpt] = useState<string>("");
  const [descriptionOls, setDescriptionOls] = useState<string>("");
  const [descriptionOlsReference, setDescriptionOlsReference] =
    useState<string>("");

  const router = useRouter();
  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";

  const { data: rawDescriptionGpt } = useDescription(cellTypeId);

  useEffect(() => {
    if (rawDescriptionGpt) setDescriptionGpt(rawDescriptionGpt);
    else setDescriptionGpt("");

    // hardcoding this for now
    const olsUrl = `
      https://www.ebi.ac.uk/ols4/api/ontologies/hcao/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}
    `;
    if (cellTypeIdRaw) {
      fetch(olsUrl).then(async (res) => {
        const data = await res.json();
        if (data) {
          if (data?.obo_definition_citation) {
            const descriptionData = data.obo_definition_citation[0];
            const description = descriptionData.definition;
            const reference = descriptionData.oboXrefs[0].url;
            if (reference) setDescriptionOlsReference(reference);
            else setDescriptionOlsReference("");
            setDescriptionOls(description);
          } else {
            setDescriptionOls("");
            setDescriptionOlsReference("");
          }
        }
      });
    }
  }, [cellTypeIdRaw, rawDescriptionGpt]);

  return (
    <>
      {descriptionOls && (
        <CellCardDescription>
          <b>
            <i>
              The below summary is sourced from the{" "}
              <a href="https://www.ebi.ac.uk/ols/ontologies/cl" target="_blank">
                Ontology Lookup Service
              </a>{" "}
              (OLS v4).{"\n\n"}
            </i>
          </b>
          {descriptionOls}
          {" ("}
          {descriptionOlsReference && (
            <a href={descriptionOlsReference} target="_blank">
              {"Reference"}
            </a>
          )}
          {")"}
        </CellCardDescription>
      )}
      <br />
      <CellCardDescription>
        <b>
          <i>
            {
              "The below summary is produced by ChatGPT. ChatGPT may produce inaccurate information about people, places, or facts. Alternate sources can be selected above.\n\n"
            }
          </i>
        </b>
        {descriptionGpt}
      </CellCardDescription>
    </>
  );
}
