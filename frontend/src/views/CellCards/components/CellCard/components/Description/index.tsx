import React, { useEffect, useState } from "react";
import { useRouter } from "next/router";
import { CellCardDescription } from "./style";
import { allCellTypeDescriptions } from "src/views/CellCards/common/fixtures";

// enum of available descriptions
type DescriptionOptions = "GPT3.5" | "Wikipedia" | "OLS v4";

interface Props {
  selectedDescription: DescriptionOptions;
  descriptions: {
    descriptionGpt: string;
    descriptionWiki: string;
    descriptionOls: string;
    descriptionOlsReference: string;
  };
  setDescriptions: {
    setDescriptionGpt: (descriptionGpt: string) => void;
    setDescriptionWiki: (descriptionWiki: string) => void;
    setDescriptionOls: (descriptionOls: string) => void;
    setDescriptionOlsReference: (descriptionOlsReference: string) => void;
  };
}

export default function Description({
  selectedDescription,
  descriptions,
  setDescriptions,
}: Props): JSX.Element {
  const router = useRouter();
  const {
    descriptionGpt,
    descriptionWiki,
    descriptionOls,
    descriptionOlsReference,
  } = descriptions;
  const {
    setDescriptionGpt,
    setDescriptionWiki,
    setDescriptionOls,
    setDescriptionOlsReference,
  } = setDescriptions;

  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";

  useEffect(() => {
    // hardcoding this for now.
    const olsUrl = `
      https://www.ebi.ac.uk/ols4/api/ontologies/hcao/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}
    `;
    if (cellTypeIdRaw) {
      // set all descriptions
      setDescriptionGpt(
        allCellTypeDescriptions[
          cellTypeId as keyof typeof allCellTypeDescriptions
        ]
      );
      fetch(`/api/scrape?cellTypeId=${cellTypeIdRaw}`).then(async (res) => {
        const data = await res.json();
        if (data.content !== "Try OLS") {
          setDescriptionWiki(data.content);
        } else {
          setDescriptionWiki("");
        }
      });
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
  }, [cellTypeIdRaw]);

  return (
    <>
      {selectedDescription === "GPT3.5" && (
        <CellCardDescription>
          <b>
            <i>
              {
                "These summaries are produced by ChatGPT. ChatGPT may produce inaccurate information about people, places, or facts. Alternate sources can be selected above.\n\n"
              }
            </i>
          </b>
          {descriptionGpt}
        </CellCardDescription>
      )}
      {descriptionOls && selectedDescription === "OLS v4" && (
        <CellCardDescription>
          {descriptionOls}{" "}
          {descriptionOlsReference && (
            <a href={descriptionOlsReference} target="_blank">
              {" (Reference)"}{" "}
            </a>
          )}
        </CellCardDescription>
      )}
      {descriptionWiki && selectedDescription === "Wikipedia" && (
        <CellCardDescription>{descriptionWiki}</CellCardDescription>
      )}
    </>
  );
}
