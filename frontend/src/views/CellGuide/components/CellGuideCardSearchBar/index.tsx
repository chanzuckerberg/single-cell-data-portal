import { useState, useMemo } from "react";
import { TextField } from "@mui/material";
import SearchIcon from "@mui/icons-material/Search";
import InputAdornment from "@mui/material/InputAdornment";
import { SectionItem, SectionTitle, StyledAutocomplete } from "./style";
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
} from "src/views/CellGuide/components/CellGuideCardSearchBar/constants";
import { SKINNY_MODE_BREAKPOINT_WIDTH } from "../CellGuideCard/constants";

interface Entity {
  id: string;
  label: string;
  synonyms?: string[]; // only cell types (optionally) have synonyms
}

const TISSUE_PREFIX = "UBERON:";

export default function CellGuideCardSearchBar({
  autoFocus = false,
}: {
  autoFocus?: boolean;
}): JSX.Element {
  const router = useRouter();
  const { data: cellTypes } = useCellTypeMetadata();
  const { data: tissueData } = useTissueMetadata();

  const options: Entity[] = useMemo(() => {
    if (!tissueData || !cellTypes) return [];
    const entities: Entity[] = [];
    for (const cellType in cellTypes) {
      entities.push({
        label: cellTypes[cellType].name,
        ...cellTypes[cellType],
      });
    }
    for (const tissue in tissueData) {
      entities.push({
        label: tissueData[tissue].name,
        ...tissueData[tissue],
      });
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
        // This is used to style the autocomplete dropdown for mobile
        sx={{
          [`@media (max-width: ${SKINNY_MODE_BREAKPOINT_WIDTH}px)`]: {
            "& + .MuiAutocomplete-popper": {
              width: "100% !important",
            },
            "& + .MuiAutocomplete-popper .MuiPaper-root": {
              boxShadow: "0 4px 4px 0 rgba(0,0,0, 0.25)", // Hides top shadow for seamless look
            },
          },
        }}
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
            autoFocus={autoFocus}
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
                if (!entity.id.startsWith(TISSUE_PREFIX)) {
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
              // Determine if each item starts with any of the synonyms
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
