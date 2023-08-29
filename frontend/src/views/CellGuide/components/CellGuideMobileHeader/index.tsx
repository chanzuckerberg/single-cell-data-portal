import { useState } from "react";
import CellGuideCardSearchBar from "../CellGuideCardSearchBar";
import { ButtonIcon } from "@czi-sds/components";
import {
  MobileHeader,
  MobileHeaderWrapper,
  MobilePageNavWrapper,
  MobileSearchBarWrapper,
  SearchBarWrapper,
  StyledCancelButton,
  StyledTitle,
} from "./style";
import { MobileSearchTint } from "../CellGuideCard/style";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

interface Props {
  title: string;
  pageNav: JSX.Element | null;
  openSearch?: boolean; // Used on landing page
  top?: number; // Used on landing page
}

const CellGuideMobileHeader = ({
  title = "",
  pageNav,
  openSearch = false,
  top = HEADER_HEIGHT_PX,
}: Props) => {
  const [searchIsOpen, setSearchIsOpen] = useState(openSearch);
  const [pageNavIsOpen, setPageNavIsOpen] = useState(false);

  const search = (
    <MobileSearchBarWrapper
      onBlur={() => {
        // Hide the input box on blur
        setSearchIsOpen(false);
      }}
    >
      <SearchBarWrapper>
        <CellGuideCardSearchBar autoFocus />
      </SearchBarWrapper>

      <StyledCancelButton
        onClick={() => {
          setSearchIsOpen(false);
        }}
      >
        Cancel
      </StyledCancelButton>
    </MobileSearchBarWrapper>
  );

  const navigation = (
    <>
      {/* Flex Item Left */}
      <div id="cellguide-search-icon">
        <ButtonIcon
          sdsIcon="search"
          id="cellguide-search-icon"
          onClick={() => {
            setSearchIsOpen(true);
          }}
        />
      </div>

      {/* Flex Item Middle */}
      <StyledTitle id="cellguide-topic">{title}</StyledTitle>

      {/* Flex Item Right */}
      <div id="cellguide-nav-dropdown">
        <ButtonIcon
          sdsIcon={pageNavIsOpen ? "chevronUp" : "chevronDown"}
          onClick={() => {
            setPageNavIsOpen(!pageNavIsOpen);
          }}
        />
      </div>
    </>
  );

  return (
    <>
      {/* Screen tint when mobile search is open */}
      {searchIsOpen && <MobileSearchTint />}

      <MobileHeaderWrapper top={top}>
        {/* CellGuide Header - First Row */}
        <MobileHeader>{searchIsOpen ? search : navigation}</MobileHeader>

        {/* CellGuide Page Nav - Second Row */}
        <MobilePageNavWrapper
          onClick={() => {
            // Close nav when clicking on a section
            setPageNavIsOpen(false);
          }}
        >
          {pageNavIsOpen && pageNav}
        </MobilePageNavWrapper>
      </MobileHeaderWrapper>
    </>
  );
};

export default CellGuideMobileHeader;
