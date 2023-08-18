import { useState, useMemo } from "react";
import { TextField } from "@mui/material";
import SearchIcon from "@mui/icons-material/Search";
import InputAdornment from "@mui/material/InputAdornment";
import { SectionItem, SectionTitle, StyledAutocomplete } from "./style";
import { ROUTES } from "src/common/constants/routes";
import { useCellGuide, useTissueCards } from "src/common/queries/cellGuide";
import { useRouter } from "next/router";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  CELL_GUIDE_CARD_SEARCH_BAR,
  CELL_GUIDE_CARD_SEARCH_BAR_TEXT_INPUT,
} from "src/views/CellGuide/components/CellGuideCardSearchBar/constants";

interface Entity {
  id: string;
  cell_types?: string[]; // only tissues have cell_types
  label: string;
  synonyms?: string[]; // only cell types (optionally) have synonyms
}

export default function CellGuideCardSearchBar(): JSX.Element {
  const router = useRouter();
  const { data: tissueData } = useTissueCards();
  const { data: cellTypes } = useCellGuide();

  const options: Entity[] = useMemo(() => {
    if (!tissueData || !cellTypes) return [];
    const entities: Entity[] = [];
    for (const cellType of cellTypes) {
      entities.push(cellType);
    }
    for (const tissue of tissueData) {
      entities.push(tissue);
    }
    return entities;
  }, [tissueData, cellTypes]);

  const [open, setOpen] = useState(false);
  const [inputValue, setValue] = useState("");

  // Used for keyboard navigation for cell type search
  const [highlightedEntity, setHighlightedEntity] = useState<Entity | null>(
    null
  );

  const handleFocus = () => {
    setOpen(true);
  };

  const handleBlur = () => {
    setOpen(false);
  };

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
      setOpen(false);
    }
  }
  return (
    <div data-testid={CELL_GUIDE_CARD_SEARCH_BAR}>
      <StyledAutocomplete
        open={open}
        value={inputValue}
        onChange={() => {
          // Clears the input after selection
          setValue("");
        }}
        onKeyDown={(event) => {
          if (highlightedEntity && event.key === "Enter") {
            changeEntity(highlightedEntity.id);
          }
        }}
        onHighlightChange={(_, value) => {
          const entity = value as Entity;
          setHighlightedEntity(entity);
        }}
        disablePortal
        id={CELL_GUIDE_CARD_SEARCH_BAR}
        options={options}
        groupBy={(option) => {
          const entity = option as Entity;
          return entity.id.split(":").at(0) ?? "";
        }}
        renderGroup={(params) => {
          const title = params.group === "CL" ? "Cell Type" : "Tissue";
          return (
            <div key={params.key}>
              <SectionTitle> {title} </SectionTitle>
              {params.children}
            </div>
          );
        }}
        renderInput={(params) => (
          <TextField
            {...params}
            data-testid={CELL_GUIDE_CARD_SEARCH_BAR_TEXT_INPUT}
            onFocus={handleFocus}
            onBlur={handleBlur}
            InputProps={{
              ...params.InputProps,
              endAdornment: (
                <InputAdornment position="end">
                  <SearchIcon />
                </InputAdornment>
              ),
            }}
            label="Search cell types or tissues"
          />
        )}
        renderOption={(props, option) => {
          const entity = option as Entity;
          return (
            <SectionItem
              {...props}
              key={entity.id}
              onClick={() => {
                // Only tissues have cell_types
                if (!entity.cell_types) {
                  track(EVENTS.CG_SEARCH_CT_TISSUE, {
                    cell_type: entity.label,
                  });
                } else {
                  track(EVENTS.CG_SEARCH_CT_TISSUE, {
                    tissue: entity.label,
                  });
                }

                changeEntity(entity.id);
              }}
            >
              {entity.label}
            </SectionItem>
          );
        }}
        autoComplete
        filterOptions={(options, state) => {
          return options
            .filter((option) => {
              const entity = option as Entity;
              const searchTerm = state.inputValue.toLowerCase();
              // Determine if each item starts with any of the snyonyms
              const synonyms = entity.synonyms ?? [];
              const synonymStartsWithSearch = synonyms.some((synonym) =>
                synonym.toLowerCase().startsWith(searchTerm)
              );

              return (
                entity.label &&
                (entity.label.toLowerCase().includes(searchTerm) ||
                  synonymStartsWithSearch)
              );
            })
            .sort((entity1, entity2) => {
              const entityA = entity1 as Entity;
              const entityB = entity2 as Entity;
              const aRaw = entityA.label;
              const bRaw = entityB.label;
              const a = aRaw.toLowerCase();
              const b = bRaw.toLowerCase();
              const searchTerm = state.inputValue.toLowerCase();

              // Determine if each item starts with the search term
              const aStartsWithSearch = a.startsWith(searchTerm);
              const bStartsWithSearch = b.startsWith(searchTerm);

              // Determine if each item starts with "CL:"
              const isA_CL = entityA.id.startsWith("CL:");
              const isB_CL = entityB.id.startsWith("CL:");

              // First, sort by search term
              if (aStartsWithSearch && !bStartsWithSearch) {
                return -1;
              }
              if (!aStartsWithSearch && bStartsWithSearch) {
                return 1;
              }

              // If neither or both start with the search term, then sort by "CL:" vs "UBERON:"
              if (isA_CL && !isB_CL) {
                return -1;
              }
              if (!isA_CL && isB_CL) {
                return 1;
              }

              // If they are both "CL:" or both "UBERON:", then sort alphabetically
              return a.localeCompare(b);
            });
        }}
      />
    </div>
  );
}
