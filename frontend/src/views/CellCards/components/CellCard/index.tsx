import React, { useEffect, useState } from "react";
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
} from "./style";
import { useCellTypesById } from "src/common/queries/cellCards";

export default function CellCard(): JSX.Element {
  const router = useRouter();
  const [description, setDescription] = useState<string>("");
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
      fetch(`/api/scrape?cellTypeId=${cellTypeIdRaw}`).then(async (res) => {
        const data = await res.json();
        if (data.content !== "Try OLS") {
          setDescription(data.content);
        } else {
          const response = await fetch(olsUrl);
          const dataOls = await response.json();
          if (dataOls) {
            console.log(dataOls);
            if (dataOls?.annotation?.definition)
              setDescription(dataOls.annotation.definition);
            else setDescription("No description available for this cell type.");
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
      <CellCardDescription>{description}</CellCardDescription>
    </Wrapper>
  );
}
