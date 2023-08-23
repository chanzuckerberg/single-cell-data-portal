import { useState } from "react";
import CellGuideCardSearchBar from "../CellGuideCardSearchBar";
import { ButtonIcon } from "@czi-sds/components";
import {
  MobileHeader,
  MobileHeaderWrapper,
  MobilePageNavWrapper,
  MobileSearchBarWrapper,
  StyledTitle,
} from "./style";
import { MobileSearchTint } from "../CellGuideCard/style";

interface Props {
  title: string;
  pageNav: JSX.Element | null;
}

const CellGuideMobileHeader = ({ title = "", pageNav }: Props) => {
  const [searchIsOpen, setSearchIsOpen] = useState(false);
  const [pageNavIsOpen, setPageNavIsOpen] = useState(false);

  const search = (
    <MobileSearchBarWrapper
      onBlur={() => {
        // Hide the input box on blur
        setSearchIsOpen(false);
      }}
    >
      <CellGuideCardSearchBar autoFocus />
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
          sdsIcon="chevronDown"
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

      <MobileHeaderWrapper>
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
