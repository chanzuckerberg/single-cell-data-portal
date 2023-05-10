import { useEffect, useState } from "react";
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
    const handleScroll = () => {
      setSectionOffsets(
        headings.reduce<{
          [id: string]: number;
        }>((agg, heading) => {
          const target = document.getElementById(heading.id);
          return {
            ...agg,
            [heading.id]: target?.offsetTop || 0,
          };
        }, {})
      );

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
  }, [globalLayoutWrapper, headings, sectionOffsets]);

  return (
    <>
      <CellCardsSidebarWrapper>
        <StickyWrapper>
          <SearchBarWrapper>
            <CellCardSearchBar />
          </SearchBarWrapper>

          <nav>
            {/* <h2>Table of Contents</h2> */}
            <TableOfContents>
              {headings.map((heading) => (
                <StyledJumpLink
                  key={heading.id}
                  // href={`#${heading.id}`}
                  onClick={() => {
                    setActiveSection(heading.id);
                    // Jump to section with keeping in mind header and padding widths
                    if (globalLayoutWrapper) {
                      if (heading.id === "intro") {
                        // If scrolling to intro just scroll to top
                        globalLayoutWrapper.scrollTop = 0;
                      } else {
                        globalLayoutWrapper.scrollTop =
                          sectionOffsets[heading.id] - TOP_PADDING_PX;
                      }
                    }
                  }}
                  style={{
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
