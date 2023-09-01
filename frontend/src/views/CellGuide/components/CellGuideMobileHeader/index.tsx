import { useState, useRef, useEffect } from "react";
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
  pageNavIsOpen?: boolean;
  setPageNavIsOpen?: (isOpen: boolean) => void;
  openSearch?: boolean; // Used on landing page
  top?: number; // Used on landing page
}

const CellGuideMobileHeader = ({
  title = "",
  pageNav,
  openSearch = false,
  top = HEADER_HEIGHT_PX,
  pageNavIsOpen,
  setPageNavIsOpen,
}: Props) => {
  const [searchIsOpen, setSearchIsOpen] = useState(openSearch);

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
          onClick={() => setPageNavIsOpen && setPageNavIsOpen(!pageNavIsOpen)}
        />
      </div>
    </>
  );

  const navRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      let targetElement = event.target as Element;
      // this is a hack to prevent the page nav from closing when clicking on child dropdowns
      do {
        if (
          targetElement &&
          targetElement.classList?.contains("MuiPopper-root")
        ) {
          return;
        }
        // Go up the DOM
        targetElement = targetElement.parentNode as Element;
      } while (targetElement);

      if (navRef.current && !navRef.current.contains(event.target as Node)) {
        setPageNavIsOpen && setPageNavIsOpen(false);
      }
    }

    // Bind the event listener
    document.addEventListener("mousedown", handleClickOutside);
    return () => {
      // Unbind the event listener on clean up
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [navRef, setPageNavIsOpen]);

  return (
    <>
      {/* Screen tint when mobile search is open */}
      {searchIsOpen && <MobileSearchTint />}

      <MobileHeaderWrapper top={top}>
        {/* CellGuide Header - First Row */}
        <MobileHeader>{searchIsOpen ? search : navigation}</MobileHeader>

        {/* CellGuide Page Nav - Second Row */}
        <MobilePageNavWrapper ref={navRef}>
          {pageNavIsOpen && pageNav}
        </MobilePageNavWrapper>
      </MobileHeaderWrapper>
    </>
  );
};

export default CellGuideMobileHeader;
