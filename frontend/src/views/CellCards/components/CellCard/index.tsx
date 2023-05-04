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
} from "./style";
import { useCellTypesById } from "src/common/queries/cellCards";
import { allCellTypeDescriptions } from "../CellCardSearchBar/fixture";
import Table from "./components/Table";

const Link = ({ title, url }: { title: string; url: string }) => {
  return (
    <a href={url} target="_blank">
      {title}
    </a>
  );
};

interface TableRow {
  symbol: string;
  name: string;
  publication: ReactElement;
}
const tableColumns: Array<keyof TableRow> = ["symbol", "name", "publication"];
const tableRows: TableRow[] = [
  {
    symbol: "JCHAIN",
    name: "Joining chain Of multimeric IgA and IgM",
    publication: (
      <Link
        title="Gene Card"
        url="https://www.genecards.org/cgi-bin/carddisp.pl?gene=JCHAIN"
      />
    ),
  },
  {
    symbol: "MZB1",
    name: "Marginal zone B and B1 cell specific protein",
    publication: (
      <Link
        title="Gene Card"
        url="https://www.genecards.org/cgi-bin/carddisp.pl?gene=MZB1"
      />
    ),
  },
  {
    symbol: "IGKC",
    name: "Immunoglobulin kappa constant",
    publication: (
      <Link
        title="Gene Card"
        url="https://www.genecards.org/cgi-bin/carddisp.pl?gene=IGKC"
      />
    ),
  },
  {
    symbol: "IGHG1",
    name: "Immunoglobulin heavy constant gamma 1 (G1m marker)",
    publication: (
      <Link
        title="Gene Card"
        url="https://www.genecards.org/cgi-bin/carddisp.pl?gene=IGHG1"
      />
    ),
  },
  {
    symbol: "IGHA1",
    name: "Immunoglobulin heavy constant alpha 1",
    publication: (
      <Link
        title="Gene Card"
        url="https://www.genecards.org/cgi-bin/carddisp.pl?gene=IGHA1"
      />
    ),
  },
  {
    symbol: "DERL3",
    name: "Derlin 3",
    publication: (
      <Link
        title="Gene Card"
        url="https://www.genecards.org/cgi-bin/carddisp.pl?gene=DERL3"
      />
    ),
  },
];

export default function CellCard(): JSX.Element {
  const router = useRouter();
  const [descriptionGpt, setDescriptionGpt] = useState<string>("");
  const [descriptionWiki, setDescriptionWiki] = useState<string>("");
  const [descriptionOls, setDescriptionOls] = useState<string>("");
  const { cellTypeId: cellTypeIdRaw } = router.query;
  const cellTypeId = (cellTypeIdRaw as string)?.replace("_", ":") ?? "";
  const cellTypesById = useCellTypesById();
  const cellTypeName = cellTypesById[cellTypeId] ?? "";

  useEffect(() => {
    // const olsUrl = `
    //   https://www.ebi.ac.uk/ols4/api/v2/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}?lang=en
    // `;
    const olsUrl = `
     https://www.ebi.ac.uk/ols/api/ontologies/cl/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}
    `;
    if (cellTypeIdRaw) {
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
          if (data?.annotation?.definition) {
            setDescriptionOls(data.annotation.definition);
          } else {
            setDescriptionOls("");
          }
        }
      });
    }
  }, [cellTypeIdRaw]);

  // if (wordsToLink && description) {
  //   const regex = new RegExp(Object.keys(wordsToLink).join('|'), 'gi');

  //   // Custom split function that includes the matched phrases in the resulting array
  //   const splitWithMatches = (inputText: string, inputRegex: RegExp) => {
  //     const result = [];
  //     let match;
  //     let lastIndex = 0;

  //     while ((match = inputRegex.exec(inputText)) !== null) {
  //       result.push(inputText.slice(lastIndex, match.index));
  //       result.push(match[0]);
  //       lastIndex = match.index + match[0].length;
  //     }

  //     result.push(inputText.slice(lastIndex));
  //     return result;
  //   };

  //   const parts = splitWithMatches(description, regex);

  //   linkedText = parts.map((part, index) => {
  //     const url = (wordsToLink[part as keyof typeof wordsToLink] ?? "").replace(":", "_");

  //     return url ? (
  //       <Link key={index} href={`${ROUTES.CELL_CARDS}/${url}`}>
  //         {part}
  //       </Link>
  //     ) : (
  //       <React.Fragment key={index}>{part}</React.Fragment>
  //     );
  //   });
  // } else {
  //   linkedText = "This Cell Card is not available yet. Please check back later."
  // }

  const genesForShareUrl = tableRows.map((row) => row.symbol).join("%2C");

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
      </CellCardHeader>
      <Divider />
      <CellCardDescription>
        <i>{"Generated by GPT:\n"}</i>
        {descriptionGpt}
      </CellCardDescription>
      <br />
      {descriptionOls && (
        <CellCardDescription>
          <i>{"From OLS:\n"}</i>
          {descriptionOls}
        </CellCardDescription>
      )}
      <br />
      {descriptionWiki && (
        <CellCardDescription>
          <i>{"From Wikipedia:\n"}</i>
          {descriptionWiki}
        </CellCardDescription>
      )}
      <TableTitleWrapper>
        <TableTitle>Marker Genes</TableTitle>
        <WmgLink
          href={`https://cellxgene.cziscience.com/gene-expression?tissues=lung&genes=${genesForShareUrl}&ver=2`}
          target="_blank"
        >
          Open in Gene Expression
        </WmgLink>
      </TableTitleWrapper>
      <Table<TableRow> columns={tableColumns} rows={tableRows} />
    </Wrapper>
  );
}
