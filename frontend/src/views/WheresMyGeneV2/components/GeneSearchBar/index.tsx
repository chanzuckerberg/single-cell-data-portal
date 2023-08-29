import { Intent } from "@blueprintjs/core";
import { Autocomplete, DefaultAutocompleteOption } from "@czi-sds/components";
import React, {
  ReactChild,
  SyntheticEvent,
  createContext,
  useCallback,
  useContext,
  useMemo,
  useState,
} from "react";
import { EVENTS } from "src/common/analytics/events";
import { usePrimaryFilterDimensions } from "src/common/queries/wheresMyGene";
import Toast from "src/views/Collection/components/Toast";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import {
  deleteAllGenes,
  selectGenes,
} from "src/views/WheresMyGene/common/store/actions";
import { Gene } from "src/views/WheresMyGene/common/types";
import {
  ActionWrapper,
  AutocompleteWrapper,
  Container,
  StyledButtonWrapper,
  StyledClearButton,
} from "./style";

import { track } from "src/common/analytics";
import {
  AutocompleteCloseReason,
  AutocompleteInputChangeReason,
} from "@mui/base";
import { pull, uniq } from "lodash";
import { FixedSizeList, ListChildComponentProps } from "react-window";

const MAX_ITEMS_TO_SHOW = 9.5;
const LISTBOX_ITEM_HEIGHT_PX = 48;

const ListBoxContext = createContext({});

const OuterElementType = React.forwardRef<HTMLDivElement>(
  function OuterElementType(props, ref) {
    const outerProps = React.useContext(ListBoxContext);
    return <div ref={ref} {...props} {...outerProps} />;
  }
);

function rowRender(props: ListChildComponentProps) {
  const { data, index, style } = props;

  return <div style={style}>{data[index]}</div>;
}

interface ListboxProps {
  children: ReactChild;
}

const ListboxComponent = React.forwardRef<HTMLDivElement, ListboxProps>(
  function ListboxComponent(props: ListboxProps, ref) {
    const { children, ...other } = props;

    const itemData = React.Children.toArray(children);
    const itemCount = itemData.length;

    const height = Math.min(
      LISTBOX_ITEM_HEIGHT_PX * itemCount,
      LISTBOX_ITEM_HEIGHT_PX * MAX_ITEMS_TO_SHOW
    );

    return (
      <div ref={ref}>
        <ListBoxContext.Provider value={other}>
          <FixedSizeList
            height={height}
            itemCount={itemCount}
            outerElementType={OuterElementType}
            itemSize={LISTBOX_ITEM_HEIGHT_PX}
            width="100%"
            overscanCount={10}
            itemData={itemData}
          >
            {rowRender}
          </FixedSizeList>
        </ListBoxContext.Provider>
      </div>
    );
  }
);

export default function GeneSearchBar({
  className,
}: {
  className?: string;
}): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { selectedGenes, selectedOrganismId } = useContext(StateContext);

  const [hasComma, setHasComma] = useState(false);
  const [input, setInput] = useState("");
  const [open, setOpen] = useState(false);

  const { data } = usePrimaryFilterDimensions(2); //temp version 2

  const { genes: rawGenes } = data || {};

  const genes: Gene[] = useMemo(() => {
    if (!rawGenes) return [];

    return rawGenes[selectedOrganismId || ""] || [];
  }, [rawGenes, selectedOrganismId]);

  /**
   * NOTE: key is gene name in lowercase
   */
  const genesByName = useMemo(() => {
    return genes.reduce((acc, gene) => {
      return acc.set(gene.name.toLowerCase(), gene);
    }, new Map<Gene["name"], Gene>());
  }, [genes]);

  const selectedGeneOptions: Gene[] = useMemo(() => {
    return selectedGenes.map((gene: string) => {
      return genesByName.get(gene.toLowerCase()) as Gene;
    });
  }, [selectedGenes, genesByName]);

  const handleGeneNotFound = useCallback((geneName: string): void => {
    Toast.show({
      intent: Intent.DANGER,
      message: `Gene not found: ${geneName}`,
    });
  }, []);

  const handleInputChange = (
    _: SyntheticEvent<Element, Event>,
    value: string,
    reason: AutocompleteInputChangeReason
  ) => {
    if (reason === "reset") {
      return;
    }
    setInput(value);
    if (value.includes(",")) {
      if (hasComma) return;
      setHasComma(true);
    } else {
      if (!hasComma) return;
      setHasComma(false);
    }
  };

  const handleClose = (
    e: SyntheticEvent<Element, Event>,
    reason: AutocompleteCloseReason
  ) => {
    if (reason === "toggleInput") {
      return;
    }

    const { nativeEvent } = e;

    if (
      (nativeEvent instanceof FocusEvent &&
        nativeEvent.relatedTarget instanceof Element) ||
      (nativeEvent instanceof FocusEvent &&
        !(nativeEvent.relatedTarget instanceof Element)) ||
      !(nativeEvent instanceof FocusEvent)
    ) {
      setOpen((open) => !open);
    }

    setInput("");
  };

  const handleEnter = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === "Enter" && hasComma) {
      event.preventDefault();
      const newSelected = [...selectedGeneOptions];
      const parsedPaste = pull(uniq(input.split(/[ ,]+/)), "");

      parsedPaste.map((item) => {
        const newItem = genesByName.get(item.toLowerCase());
        if (!newItem) {
          handleGeneNotFound(item);
        } else if (!newSelected.includes(newItem)) newSelected.push(newItem);
      });
      setOpen((open) => !open);

      return handleSelectGenes(undefined, newSelected, undefined, undefined);
    }
  };

  return (
    <Container {...{ className }}>
      <ActionWrapper>
        <AutocompleteWrapper>
          <Autocomplete
            value={selectedGeneOptions}
            multiple
            fullWidth
            options={genes}
            inputValue={input}
            onChange={handleSelectGenes}
            onInputChange={handleInputChange}
            label="Add Genes"
            search
            onFocus={() => setOpen(true)}
            onClose={handleClose}
            onKeyDownCapture={handleEnter}
            open={open}
            keepSearchOnSelect={false}
            ListboxComponent={
              ListboxComponent as React.ComponentType<
                React.HTMLAttributes<HTMLElement>
              >
            }
            noOptionsText={
              hasComma
                ? "You can add multiple genes using a comma-separated list. Press enter to add."
                : "No options"
            }
          />
        </AutocompleteWrapper>
        {/* Clear Genes button */}
        {!selectedGenes.length || (
          <StyledButtonWrapper
            onClick={() => {
              if (dispatch) {
                track(EVENTS.WMG_CLEAR_GENES_CLICKED);

                dispatch(deleteAllGenes());
              }
            }}
          >
            <StyledClearButton
              data-testid="clear-genes-button"
              sdsType="primary"
              sdsStyle="minimal"
              isAllCaps={false}
            >
              Clear Genes
            </StyledClearButton>
          </StyledButtonWrapper>
        )}
      </ActionWrapper>
    </Container>
  );

  function handleSelectGenes(
    _: unknown,
    genes: DefaultAutocompleteOption[] | null,
    __: unknown,
    ___: unknown
  ) {
    if (!dispatch) return;

    dispatch(selectGenes(genes?.map((gene) => gene.name) || []));
  }
}
