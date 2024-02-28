import { DefaultDropdownMenuOption } from "@czi-sds/components";
import { throttle } from "lodash";
import { useRouter } from "next/router";
import { useState, useRef, useCallback, useMemo, useEffect } from "react";
import { ROUTES } from "src/common/constants/routes";
import {
  ALL_TISSUES,
  NO_ORGAN_ID,
  TISSUE_AGNOSTIC,
} from "src/views/CellGuide/components/CellGuideCard/components/MarkerGeneTables/constants";
import { useOrganAndOrganismFilterListForCellType } from "src/views/CellGuide/components/CellGuideCard/components/MarkerGeneTables/hooks/common";
import { SKINNY_MODE_BREAKPOINT_WIDTH } from "src/views/CellGuide/components/CellGuideCard/constants";
import { LEFT_RIGHT_PADDING_PX_XXL } from "src/views/CellGuide/components/CellGuideCard/style";
import { SDSOrgan } from "src/views/CellGuide/components/CellGuideCard/types";
import { CellType } from "src/views/CellGuide/components/common/OntologyDagView/common/types";
import { Gene } from "src/views/WheresMyGeneV2/common/types";

// Hardcode this for now
export const CELLGUIDE_ORGANISMS_LIST = ["Homo sapiens", "Mus musculus"];

export function useConnect() {
  const router = useRouter();

  const [pageNavIsOpen, setPageNavIsOpen] = useState(false);
  const [selectedGene, setSelectedGene] = useState<string | undefined>(
    undefined
  );

  const [skinnyMode, setSkinnyMode] = useState<boolean>(false);

  // Navigation
  const sectionRef0 = useRef<HTMLDivElement>(null);
  const sectionRef1 = useRef<HTMLDivElement>(null);
  const sectionRef2 = useRef<HTMLDivElement>(null);
  const sectionRef3 = useRef<HTMLDivElement>(null);

  const selectGene = useCallback(
    (gene: string) => {
      if (gene === selectedGene) {
        setSelectedGene(undefined);
      } else {
        setSelectedGene(gene);
        if (sectionRef1.current) {
          window.scrollTo({
            top:
              sectionRef1.current.getBoundingClientRect().top +
              window.scrollY -
              50,
            behavior: "smooth",
          });
        }
      }
    },
    [selectedGene]
  );

  // Set the mobile tooltip view content
  const [tooltipContent, setTooltipContent] = useState<{
    title: string;
    element: JSX.Element;
  } | null>(null);

  const { cellTypeId: queryCellTypeId, tissueId: queryTissueId } = router.query;
  const cellTypeId = (queryCellTypeId as string)?.replace("_", ":") ?? "";

  /**
   * (thuang): `queryTissueId` can be undefined when a user is on `/cellguide/:cellTypeId` route
   * instead of `/cellguide/tissues/:tissueId/cell-types/:cellTypeId` route.
   *
   * NOTE: `tissue` and `organ` in variable names are interchangeable here.
   */
  const tissueId = (queryTissueId as string)?.replace("_", ":") || NO_ORGAN_ID;

  const handleResize = useCallback(() => {
    setSkinnyMode(
      window.innerWidth <
        SKINNY_MODE_BREAKPOINT_WIDTH + 2 * LEFT_RIGHT_PADDING_PX_XXL
    );
  }, [setSkinnyMode]);

  const throttledHandleResize = useMemo(() => {
    return throttle(handleResize, 100);
  }, [handleResize]);

  useEffect(() => {
    throttledHandleResize();
    window.addEventListener("resize", throttledHandleResize);

    return () => window.removeEventListener("resize", throttledHandleResize);
  }, [throttledHandleResize]);

  const [geneInfoGene, setGeneInfoGene] = useState<Gene["name"] | null>(null);
  const [cellInfoCellType, setCellInfoCellType] = useState<CellType | null>(
    null
  );
  const sdsOrganismsList = CELLGUIDE_ORGANISMS_LIST.map((organism) => ({
    name: organism,
  }));
  const [selectedOrganism, setSelectedOrganism] =
    useState<DefaultDropdownMenuOption>(sdsOrganismsList[0]);

  const { organsMap, isSuccess } = useOrganAndOrganismFilterListForCellType(
    cellTypeId,
    selectedOrganism.name
  );

  const sdsOrgansList = useMemo<SDSOrgan[]>(
    () =>
      Array.from(organsMap).map(([name, id]) => ({
        name: name === ALL_TISSUES ? TISSUE_AGNOSTIC : name,
        id,
      })),
    [organsMap]
  );

  const selectedOrgan = sdsOrgansList.find((organ) => organ.id === tissueId);

  /**
   * (thuang): Push the user to tissue agnostic cell type page if the user is on
   * a tissue-specific cell type page and the tissue is not in the filter list.
   */
  useEffect(() => {
    if (isSuccess && tissueId !== NO_ORGAN_ID && !selectedOrgan) {
      router.replace(
        ROUTES.CELL_GUIDE_CELL_TYPE.replace(
          ":cellTypeId",
          queryCellTypeId as string
        )
      );
    }
  }, [queryCellTypeId, router, selectedOrgan, tissueId, isSuccess]);

  useEffect(() => {
    setSelectedGene(undefined);
  }, [selectedOrgan, selectedOrganism, setSelectedGene]);

  return {
    router,
    pageNavIsOpen,
    setPageNavIsOpen,
    selectedGene,
    sectionRef0,
    sectionRef1,
    sectionRef2,
    sectionRef3,
    skinnyMode,
    selectGene,
    tooltipContent,
    setTooltipContent,
    queryCellTypeId,
    cellTypeId,
    geneInfoGene,
    setGeneInfoGene,
    cellInfoCellType,
    setCellInfoCellType,
    organismsList: CELLGUIDE_ORGANISMS_LIST,
    organsMap,
    sdsOrganismsList,
    sdsOrgansList,
    selectedOrgan,
    selectedOrganId: tissueId,
    selectedOrganism,
    setSelectedOrganism,
  };
}
