import { useContext, useState } from "react";
import CellGuideCardSearchBar from "../CellGuideCardSearchBar";
import { ButtonIcon } from "@czi-sds/components";
import {
  MobileHeader,
  MobileHeaderWrapper,
  MobilePageNavWrapper,
  MobileSearchBarWrapper,
  StyledTitle,
} from "./style";
import { StateContext } from "../../common/store";

const CellGuideMobileHeader = () => {
  const { cellGuideTitle, cellGuideNav, skinnyMode } = useContext(StateContext);

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
      <StyledTitle id="cellguide-topic">{cellGuideTitle}</StyledTitle>

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
      {skinnyMode && (
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
            {pageNavIsOpen && cellGuideNav}
          </MobilePageNavWrapper>
        </MobileHeaderWrapper>
      )}
    </>
  );
};

export default CellGuideMobileHeader;
