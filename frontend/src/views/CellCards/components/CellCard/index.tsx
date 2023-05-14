import React, { useState, useEffect } from "react";
import { useRouter } from "next/router";
import {
  Wrapper,
  CellCardName,
  CellCardHeader,
  StyledTag,
  Divider,
  CellCardsView,
} from "./style";
import { useCellTypesById } from "src/common/queries/cellCards";
import { SelectChangeEvent } from "@mui/material/Select";
import Description from "./components/Description";
import CanonicalMarkerGeneTable from "./components/CanonicalMarkerGeneTable";
import EnrichedGenesTable from "./components/EnrichedGenesTable";
import SourceDataTable from "./components/SourceDataTable";
import CellCardSidebar, {
  INTRO_SECTION_ID,
} from "./components/CellCardSidebar";
import DropdownSelect from "./components/common/DropdownSelect";
import OntologyDagView from "./components/OntologyDagView";

// enum of available descriptions
type DescriptionOptions = "GPT3.5" | "Wikipedia" | "OLS v4";
const availableDescriptions: DescriptionOptions[] = [
  "GPT3.5",
  "Wikipedia",
  "OLS v4",
];

export default function CellCard(): JSX.Element {
  const router = useRouter();

  // descriptions
  const [selectedDescription, setSelectedDescription] =
    useState<DescriptionOptions>("GPT3.5");
  const [descriptionGpt, setDescriptionGpt] = useState<string>("");
  const [descriptionWiki, setDescriptionWiki] = useState<string>("");
  const [descriptionOls, setDescriptionOls] = useState<string>("");
  const [descriptionOlsReference, setDescriptionOlsReference] =
    useState<string>("");

  // cell type id
  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";
  const cellTypesById = useCellTypesById() ?? {};
  const cellTypeName = cellTypesById[cellTypeId] ?? "";

  const available = availableDescriptions.filter((description) => {
    if (description === "GPT3.5") {
      return descriptionGpt !== "";
    } else if (description === "Wikipedia") {
      return descriptionWiki !== "";
    } else if (description === "OLS v4") {
      return descriptionOls !== "";
    }
  });

  const handleChange = (event: SelectChangeEvent) => {
    setSelectedDescription(event.target.value as DescriptionOptions);
  };

  useEffect(() => {
    if (
      (selectedDescription === "OLS v4" && descriptionOls === "") ||
      (selectedDescription === "Wikipedia" && descriptionWiki === "")
    )
      setSelectedDescription("GPT3.5");
  }, [cellTypeId, selectedDescription, descriptionWiki, descriptionOls]);

  // For testing purposes in the prototype, we will have multiple kinds of descriptions
  // This codepath will be simplified once we settle on a particular description.
  const descriptions = {
    descriptionGpt,
    descriptionWiki,
    descriptionOls,
    descriptionOlsReference,
  };
  const setDescriptions = {
    setDescriptionGpt,
    setDescriptionWiki,
    setDescriptionOls,
    setDescriptionOlsReference,
  };

  return (
    <>
      <style>
        {`
          /* Hack because main has a global overflow CSS prop which interferes with sticky sidebar and scroll listener */
          main {
            overflow: unset !important;
          }
          html {
            height: unset !important;
          }
        `}
      </style>
      <CellCardsView>
        {/* Flex item left */}
        <Wrapper>
          <CellCardHeader id={INTRO_SECTION_ID}>
            <CellCardName>
              {cellTypeName.charAt(0).toUpperCase() + cellTypeName.slice(1)}
            </CellCardName>
            <StyledTag
              label={cellTypeId}
              sdsType="secondary"
              sdsStyle="square"
              color="gray"
              hover={false}
            />

            <DropdownSelect
              handleChange={handleChange}
              options={available}
              selectedOption={selectedDescription}
            />
          </CellCardHeader>
          <Divider />
          <Description
            selectedDescription={selectedDescription}
            descriptions={descriptions}
            setDescriptions={setDescriptions}
          />
          <OntologyDagView cellTypeId={cellTypeId} width={1000} height={500} />
          <CanonicalMarkerGeneTable cellTypeId={cellTypeId} />
          <EnrichedGenesTable cellTypeId={cellTypeId} />
          <SourceDataTable cellTypeId={cellTypeId} />
        </Wrapper>

        {/* Flex item right */}
        <CellCardSidebar />
      </CellCardsView>
    </>
  );
}
