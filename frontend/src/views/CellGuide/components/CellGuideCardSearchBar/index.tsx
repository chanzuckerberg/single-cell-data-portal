import { useMemo, useState } from "react";
import { StyledPopper, SectionItem } from "./style";
import { ROUTES } from "src/common/constants/routes";
import {
  useCellTypeMetadata,
  useTissueMetadata,
} from "src/common/queries/cellGuide";
import { useRouter } from "next/router";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  CELL_GUIDE_CARD_SEARCH_BAR,
  CELL_GUIDE_CARD_SEARCH_BAR_TEXT_INPUT,
  CELL_GUIDE_SEARCH_BAR_LABEL,
} from "src/views/CellGuide/components/CellGuideCardSearchBar/constants";
import { Autocomplete, SDSAutocompleteOptions } from "@czi-sds/components";
import { Paper } from "@mui/material";

export interface Entity {
  id: string;
  name: string;
  synonyms?: string[]; // only cell types (optionally) have synonyms
}

const TISSUE_PREFIX = "UBERON:";

export default function CellGuideCardSearchBar({
  skinnyMode,
}: {
  skinnyMode: boolean;
}): JSX.Element {
  const router = useRouter();
  const [open, setOpen] = useState(false);
  const { data: cellTypes } = useCellTypeMetadata();
  const { data: tissueData } = useTissueMetadata();
  const options: SDSAutocompleteOptions<Entity, false, false, false> =
    useMemo(() => {
      if (!tissueData || !cellTypes) return [];
      const entities1: Entity[] = [];
      for (const cellType in cellTypes) {
        entities1.push(cellTypes[cellType]);
      }
      const entities2: Entity[] = [];
      for (const tissue in tissueData) {
        entities2.push(tissueData[tissue]);
      }
      return [
        { name: "Cell Types", options: entities1 },
        { name: "Tissues", options: entities2 },
      ];
    }, [tissueData, cellTypes]);

  function changeEntity(entityId: string) {
    if (!entityId) return;

    const formattedEntityId = entityId.replace(":", "_");

    if (entityId) {
      if (entityId.startsWith("CL:")) {
        router.push(`${ROUTES.CELL_GUIDE}/${formattedEntityId}`);
      } else {
        router.push(`${ROUTES.CELL_GUIDE}/tissues/${formattedEntityId}`);
      }
      document.getElementById(CELL_GUIDE_CARD_SEARCH_BAR)?.blur();
    }
  }
  return (
    <div data-testid={CELL_GUIDE_CARD_SEARCH_BAR}>
      <Autocomplete<Entity, false, false, false>
        data-testid={CELL_GUIDE_CARD_SEARCH_BAR_TEXT_INPUT}
        label={CELL_GUIDE_SEARCH_BAR_LABEL}
        search
        InputBaseProps={{
          onFocus: () => {
            setOpen(true);
          },
          onBlur: () => {
            setOpen(false);
          },
        }}
        PopperComponent={(popperProps) => (
          <StyledPopper
            {...popperProps}
            disablePortal
            open={open}
            fullWidth={skinnyMode}
          />
        )}
        id={CELL_GUIDE_CARD_SEARCH_BAR}
        options={options}
        clearOnEscape
        renderOption={(props, option) => {
          return (
            <SectionItem
              {...props}
              key={option.id}
              onClick={() => {
                if (!option.id.startsWith(TISSUE_PREFIX)) {
                  track(EVENTS.CG_SEARCH_CT_TISSUE, {
                    cell_type: option.name,
                  });
                } else {
                  track(EVENTS.CG_SEARCH_CT_TISSUE, {
                    tissue: option.name,
                  });
                }

                changeEntity(option.id);
              }}
            >
              {option.name}
            </SectionItem>
          );
        }}
        filterOptions={(options, state) => {
          return options
            .filter((option) => {
              const entity = option as Entity;
              const searchTerm = state.inputValue.toLowerCase();
              // Determine if each item starts with any of the synonyms
              const synonyms = entity.synonyms ?? [];
              const synonymStartsWithSearch = synonyms.some((synonym) =>
                synonym.toLowerCase().startsWith(searchTerm)
              );

              return (
                entity.name &&
                (entity.name.toLowerCase().includes(searchTerm) ||
                  entity.id.toLowerCase().includes(searchTerm) ||
                  synonymStartsWithSearch)
              );
            })
            .sort((entity1, entity2) => {
              const entityA = entity1 as Entity;
              const entityB = entity2 as Entity;
              const aRaw = entityA.name;
              const bRaw = entityB.name;
              const a = aRaw.toLowerCase();
              const b = bRaw.toLowerCase();
              const searchTerm = state.inputValue.toLowerCase();

              // Determine if each item starts with the search term
              const aStartsWithSearch = a.startsWith(searchTerm);
              const bStartsWithSearch = b.startsWith(searchTerm);

              // First, sort by search term
              if (aStartsWithSearch && !bStartsWithSearch) {
                return -1;
              }
              if (!aStartsWithSearch && bStartsWithSearch) {
                return 1;
              }

              return a.localeCompare(b);
            });
        }}
      />
    </div>
  );
}
