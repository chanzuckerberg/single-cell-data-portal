import React, { ReactElement, useEffect, useState } from "react";
import { useRouter } from "next/router";
import CellCardSearchBar from "../CellCardSearchBar";
import {
  Wrapper,
  SearchBarWrapper,
  CellCardName,
  CellCardHeader,
  StyledTag,
  Divider,
  CellCardDescription,
  TableTitle,
  TableTitleWrapper,
  WmgLink,
  PublicationLinkWrapper,
} from "./style";
import { useCellTypesById } from "src/common/queries/cellCards";
import {
  allCellTypeDescriptions,
  allCellTypeMarkerGenes,
} from "src/views/CellCards/common/fixtures";
import Table from "./components/Table";
import Select, { SelectChangeEvent } from "@mui/material/Select";
import MenuItem from "@mui/material/MenuItem";

const Link = ({ title, url }: { title: string; url: string }) => {
  return (
    <a href={url} target="_blank">
      {title}
    </a>
  );
};

// enum of available descriptions
type DescriptionOptions = "GPT3.5" | "Wikipedia" | "OLS v4";
const availableDescriptions: DescriptionOptions[] = [
  "GPT3.5",
  "Wikipedia",
  "OLS v4",
];

interface TableRow {
  symbol: ReactElement;
  name: string;
  publications: ReactElement | string;
  organ: string;
  tissue: string;
}
const tableColumns: Array<keyof TableRow> = [
  "symbol",
  "name",
  "organ",
  "tissue",
  "publications",
];

export default function CellCard(): JSX.Element {
  const router = useRouter();
  const [selectedDescription, setSelectedDescription] =
    useState<DescriptionOptions>("GPT3.5");

  const [descriptionGpt, setDescriptionGpt] = useState<string>("");
  const [descriptionWiki, setDescriptionWiki] = useState<string>("");
  const [descriptionOls, setDescriptionOls] = useState<string>("");
  const [descriptionOlsReference, setDescriptionOlsReference] =
    useState<string>("");
  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";
  const cellTypesById = useCellTypesById();
  const cellTypeName = cellTypesById[cellTypeId] ?? "";

  const tableRows: TableRow[] = [];
  if (cellTypeId in allCellTypeMarkerGenes) {
    const genes =
      allCellTypeMarkerGenes[cellTypeId as keyof typeof allCellTypeMarkerGenes];
    for (const markerGene of genes) {
      const publications = markerGene.publication.split(";;");
      const publicationTitles = markerGene.publication_titles.split(";;");
      const publicationLinks = (
        <PublicationLinkWrapper>
          {publications.map((publication, index) => {
            if (publication && publicationTitles[index]) {
              return (
                <Link
                  title={publicationTitles[index]}
                  url={`https://doi.org/${publication}`}
                />
              );
            }
          })}
        </PublicationLinkWrapper>
      );

      tableRows.push({
        symbol: (
          <Link
            title={`${markerGene.symbol}`}
            url={`https://www.genecards.org/cgi-bin/carddisp.pl?gene=${markerGene.symbol}`}
          />
        ),
        name: markerGene.name,
        publications: publicationLinks,
        organ: markerGene.tissue_general,
        tissue: markerGene.tissue_specific,
      });
    }
  }
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

  const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");
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
      {selectedDescription === "GPT3.5" && (
        <CellCardDescription>
          <i>
            {
              "ChatGPT may produce inaccurate information about people, places, or facts.\n\n"
            }
          </i>
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
      <TableTitleWrapper>
        <TableTitle>Marker Genes</TableTitle>
        {tableRows.length > 0 && (
          <WmgLink
            href={`https://cellxgene.cziscience.com/gene-expression?tissues=lung&genes=${genesForShareUrl}&ver=2`}
            target="_blank"
          >
            Open in Gene Expression
          </WmgLink>
        )}
      </TableTitleWrapper>
      {tableRows.length ? (
        <Table<TableRow> columns={tableColumns} rows={tableRows} />
      ) : (
        <div>Canonical marker genes are not available yet.</div>
      )}
    </Wrapper>
  );
}
