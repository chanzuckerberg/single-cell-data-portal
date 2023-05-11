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

interface SectionOffset {
  [id: string]: number;
}

interface Section {
  id: string;
  title: string;
}

export const INTRO_SECTION_ID = "intro-section";
export const MARKER_GENES_SECTION_ID = "marker-genes-section";
export const HIGHLY_EXPRESSED_GENES_SECTION_ID =
  "highly-expressed-genes-section";
export const SOURCE_DATA_SECTION_ID = "source-data-section";

const SECTIONS: Section[] = [
  { id: INTRO_SECTION_ID, title: "Intro" },
  { id: MARKER_GENES_SECTION_ID, title: "Marker Genes" },
  { id: HIGHLY_EXPRESSED_GENES_SECTION_ID, title: "Highly Expressed Genes" },
  { id: SOURCE_DATA_SECTION_ID, title: "Source Data" },
];

let scrollTimeout: number | null = null;

export default function CellCardSidebar(): JSX.Element {
  const [activeSectionId, setActiveSection] = useState(INTRO_SECTION_ID);
  const [targetActiveSectionId, setTargetActiveSection] =
    useState(INTRO_SECTION_ID);

  const globalLayoutWrapper = !isSSR()
    ? document.getElementById("global-layout-wrapper")
    : null;

  // map header id => section's scroll offset on page
  function calculateSectionOffsets(sections: Section[]) {
    return sections.reduce<SectionOffset>((agg, section) => {
      const target = document.getElementById(section.id);
      return {
        ...agg,
        [section.id]: target?.offsetTop || 0,
      };
    }, {});
  }

  useEffect(() => {
    const handleScroll = () => {
      if (globalLayoutWrapper) {
        const isScrollAtBottom =
          (globalLayoutWrapper.scrollTop + globalLayoutWrapper.clientHeight ??
            0) >= globalLayoutWrapper.scrollHeight;
        if (scrollTimeout !== null) {
          window.clearTimeout(scrollTimeout);
        }

        // Set a timeout to run after scrolling ends
        scrollTimeout = window.setTimeout(() => {
          if (isScrollAtBottom) {
            setActiveSection(SECTIONS[SECTIONS.length - 1].id);
            setTargetActiveSection(SECTIONS[SECTIONS.length - 1].id);
          } else {
            setActiveSection(targetActiveSectionId);
            setTargetActiveSection(targetActiveSectionId);
          }
        }, 100);

        if (activeSectionId === targetActiveSectionId) {
          const sectionOffsets = calculateSectionOffsets(SECTIONS);

          let currentScrollPosition = globalLayoutWrapper.scrollTop || 0;

          // Computes the offset taking into account the header and padding height
          currentScrollPosition += HEADER_HEIGHT_PX + TOP_PADDING_PX;

          // Calculate what section we're currently on
          const sectionOffsetsArray = Object.keys(sectionOffsets);
          const newActiveSectionId =
            sectionOffsetsArray.find((id, i) => {
              const sectionOffset = sectionOffsets[id];
              const nextSectionOffset =
                sectionOffsets[sectionOffsetsArray[i + 1]]; // Get the offset of the next section

              return (
                (currentScrollPosition >= sectionOffset &&
                  nextSectionOffset &&
                  currentScrollPosition < nextSectionOffset) ||
                !nextSectionOffset
              );
            }) || "";
          if (isScrollAtBottom) {
            setActiveSection(
              sectionOffsetsArray[sectionOffsetsArray.length - 1]
            );
            setTargetActiveSection(
              sectionOffsetsArray[sectionOffsetsArray.length - 1]
            );
          } else {
            setActiveSection(newActiveSectionId);
            setTargetActiveSection(newActiveSectionId);
          }
        }
      }
    };

    window.addEventListener("scroll", handleScroll, true);

    // Unmount cleanup
    return () => {
      window.removeEventListener("scroll", handleScroll, true);
    };
  }, [globalLayoutWrapper, activeSectionId, targetActiveSectionId]);

  return (
    <>
      <CellCardsSidebarWrapper>
        <StickyWrapper>
          <SearchBarWrapper>
            <CellCardSearchBar />
          </SearchBarWrapper>

          <nav>
            <TableOfContents>
              {SECTIONS.map((section) => (
                <StyledJumpLink
                  isActive={section.id === targetActiveSectionId}
                  key={section.id}
                  onClick={() => {
                    setTargetActiveSection(section.id);
                    const sectionOffsets = calculateSectionOffsets(SECTIONS);
                    if (globalLayoutWrapper) {
                      // Jumping like this will also trigger the scroll listener
                      if (section.id === INTRO_SECTION_ID) {
                        // If scrolling to intro just scroll to top
                        globalLayoutWrapper.scrollTo({
                          top: 0,
                          behavior: "smooth",
                        });
                      } else {
                        // Jump to section but add some top padding
                        globalLayoutWrapper.scrollTo({
                          top:
                            sectionOffsets[section.id] -
                            TOP_PADDING_PX -
                            HEADER_HEIGHT_PX,
                          behavior: "smooth",
                        });
                      }
                    }
                  }}
                >
                  {section.title}
                </StyledJumpLink>
              ))}
            </TableOfContents>
          </nav>
        </StickyWrapper>
      </CellCardsSidebarWrapper>
    </>
  );
}
