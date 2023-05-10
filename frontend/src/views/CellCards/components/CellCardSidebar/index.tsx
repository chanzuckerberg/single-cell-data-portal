import { useEffect, useState } from "react";
import { noop } from "src/common/constants/utils";
import { isSSR } from "src/common/utils/isSSR";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { SearchBarWrapper, TOP_PADDING_PX } from "../CellCard/style";
import CellCardSearchBar from "../CellCardSearchBar";
import {
  CellCardsSidebarWrapper,
  StickyWrapper,
  StyledJumpLink,
  TableOfContents,
} from "./style";

export default function CellCardSidebar({
  headings,
}: {
  headings: { id: string; title: string }[];
}): JSX.Element {
  const [activeSection, setActiveSection] = useState("intro");
  const [sectionOffsets, setSectionOffsets] = useState<{
    [id: string]: number;
  }>({});

  const globalLayoutWrapper = !isSSR()
    ? document.getElementById("global-layout-wrapper")
    : null;

  useEffect(() => {
    // map header id => section's scroll offset on page
    const sectionOffsets = headings.reduce<{
      [id: string]: number;
    }>((agg, heading) => {
      const target = document.getElementById(heading.id);
      return {
        ...agg,
        [heading.id]: target?.offsetTop || 0,
      };
    }, {});

    setSectionOffsets(sectionOffsets);

    const handleScroll = () => {
      console.log("triggered");
      let currentScrollPosition = globalLayoutWrapper?.scrollTop || 0;
      currentScrollPosition += HEADER_HEIGHT_PX + TOP_PADDING_PX;

      const sectionOffsetsArray = Object.keys(sectionOffsets);

      sectionOffsetsArray.forEach((id, i) => {
        const sectionOffset = sectionOffsets[id];
        const nextSectionOffset = sectionOffsets[sectionOffsetsArray[i + 1]];

        // Check if current scroll position is within a section but not before the next section
        if (
          currentScrollPosition >= sectionOffset &&
          nextSectionOffset &&
          currentScrollPosition < nextSectionOffset
        ) {
          setActiveSection(id);
        }
      });
    };

    window.addEventListener("scroll", handleScroll, true);

    // Unmount cleanup
    return () => {
      window.removeEventListener("scroll", handleScroll, true);
    };
  }, [globalLayoutWrapper?.scrollTop, headings]);

  return (
    <>
      <CellCardsSidebarWrapper>
        <StickyWrapper>
          <SearchBarWrapper>
            <CellCardSearchBar />
          </SearchBarWrapper>

          <nav>
            <TableOfContents>
              {headings.map((heading) => (
                <StyledJumpLink
                  key={heading.id}
                  onClick={() => {
                    // Jump to section with keeping in mind header and padding widths
                    // Jumping like this will also trigger the scroll listener
                    if (globalLayoutWrapper) {
                      if (heading.id === "intro") {
                        // If scrolling to intro just scroll to top
                        globalLayoutWrapper.scrollTop = 0;
                      } else {
                        globalLayoutWrapper.scrollTop =
                          sectionOffsets[heading.id] - TOP_PADDING_PX;
                      }
                    }

                    // Prevents race condition from scroll listener being triggered on section jump
                    setTimeout(() => {
                      setActiveSection(heading.id);
                    }, 50);
                  }}
                  style={{
                    scrollMarginTop: TOP_PADDING_PX + HEADER_HEIGHT_PX,
                    fontFamily: "Inter",
                    letterSpacing: "-0.006em",
                    lineHeight: "20px",
                    color: `${
                      heading.id === activeSection ? "black" : "#767676"
                    }`,
                    fontWeight: `${
                      heading.id === activeSection ? "600" : "500"
                    }`,
                    borderLeft: `2px solid ${
                      heading.id === activeSection ? "#0073ff" : "#EAEAEA"
                    }`,
                  }}
                >
                  {heading.title}
                </StyledJumpLink>
              ))}
            </TableOfContents>
          </nav>
        </StickyWrapper>
      </CellCardsSidebarWrapper>
    </>
  );
}
