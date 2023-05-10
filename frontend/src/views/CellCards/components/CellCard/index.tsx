import React, { useState } from "react";
import { useRouter } from "next/router";
import CellCardSearchBar from "../CellCardSearchBar";
import {
  Wrapper,
  SearchBarWrapper,
  CellCardName,
  CellCardHeader,
  StyledTag,
  Divider,
} from "./style";
import { useCellTypesById } from "src/common/queries/cellCards";
import Select, { SelectChangeEvent } from "@mui/material/Select";
import MenuItem from "@mui/material/MenuItem";
import Description from "./components/Description";
import CanonicalMarkerGeneTable from "./components/CanonicalMarkerGeneTable";

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
  const cellTypesById = useCellTypesById();
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
    <Wrapper>
      <SearchBarWrapper>
        <CellCardSearchBar />
      </SearchBarWrapper>
      <CellCardHeader>
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
        <Select
          labelId="dropdown-label"
          id="dropdown"
          value={selectedDescription}
          onChange={handleChange}
          label="Select description"
        >
          {available.map((description) => (
            <MenuItem value={description}>{description}</MenuItem>
          ))}
        </Select>
      </CellCardHeader>
      <Divider />
      <Description
        selectedDescription={selectedDescription}
        descriptions={descriptions}
        setDescriptions={setDescriptions}
      />
      <CanonicalMarkerGeneTable cellTypeId={cellTypeId} />
    </Wrapper>
  );
}
